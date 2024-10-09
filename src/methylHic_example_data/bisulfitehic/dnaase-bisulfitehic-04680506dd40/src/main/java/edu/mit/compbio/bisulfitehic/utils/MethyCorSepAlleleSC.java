/**
 * MethyCorAcrossHiccups.java
 * Sep 12, 2016
 * 1:03:54 PM
 * yaping    lyping1986@gmail.com
 */
package main.java.edu.mit.compbio.bisulfitehic.utils;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.IntervalTree.Node;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.exception.UnsortedFileException;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import main.java.edu.mit.compbio.ccinference.utils.CcInferenceUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;

import java.util.zip.GZIPOutputStream;

import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.util.Pair;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;


/**
 * accept BEDPE file. for each rows, locate reads from paired end bam file within the region, also locate SNP from VCF file within the region
 * for region within 20kb, skip it
 * iterate each bedpe region, first list all Cpg pairs that have interaction within two regions in any of the cells. by using bams.
 * for each cpg pairs, store a list,  i position in list , 1 represent i cells have interactions for these cpg pairs, 0 represent no interactions.
 *  only keep cpg pairs have at least two cells. then using these cpg_pairs at cell , extract methy value from related bw files. output each line as:
 *  cpg1_chr, cpg1_start, cpg1_end, cpg2_chr, cpg2_start, cpg2_end, region_name, cell_name, ob, exp, fdr, methy1, numC+numT1, methy2, numC+numT2
 * 
 * Input: bam file should be sorted, calmd, and read name for end1 and end2 should be the same
 * TODO: after get correct calmd bam file, need to change complementArray(bases) in negative strand case
 */

public class MethyCorSepAlleleSC {

	@Option(name="-minBaseQ",usage="minimum base quality score required to check. Default: 5")
	public int minBaseQ = 5;
	
	@Option(name="-bams",usage="input bams for each single cell. Default: null")
	public ArrayList<String> bams = null;
	
	@Option(name="-minCellNum",usage="minimum number of cell have methylation value there. Default: 2")
	public int minCellNum = 2;
	
	
	@Option(name="-useBadMate",usage="use not proper paired reads. Default: false")
	public boolean useBadMate = false;
	
	@Option(name="-minMapQ",usage="minimum mapping quality score required to check. Default: 20")
	public int minMapQ = 20;

	@Option(name="-turnOffBisulfiteFilter",usage="turn off bisulfite conversion reads filter. Default: false")
	public boolean turnOffBisulfiteFilter = false;


	
	@Option(name="-skipProcessFirstRow",usage="don't process first row, default: not enabled")
	public boolean skipProcessFirstRow = false;


	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "MethyCorSepAlleleSC [opts] outputPrefix snp_vcf.gz";
	
	private static Logger log = Logger.getLogger(MethyCorSepAlleleSC.class);

	private OutputStreamWriter methyCorWriter = null; 
	private OutputStreamWriter poissonBinomialWriter = null;
	
	private static long startTime = -1;

	
	/**
	 * @param args
	 */
	public static void main(String[] args)
		throws Exception {
		MethyCorSepAlleleSC mcah = new MethyCorSepAlleleSC();
			BasicConfigurator.configure();
			mcah.doMain(args);
		}
		

		public void doMain(String[] args)
				throws Exception {

						CmdLineParser parser = new CmdLineParser(this);
						//parser.setUsageWidth(80);
						try
						{
							if(help || args.length < 2) throw new CmdLineException(parser, USAGE, new Throwable());
							parser.parseArgument(args);
							
						
						}
						catch (CmdLineException e)
						{
							System.err.println(e.getMessage());
							// print the list of available options
							parser.printUsage(System.err);
							System.err.println();
							return;
						}
						if(bams==null  
								|| bams.isEmpty() ) 
							throw new CmdLineException(parser, USAGE, new Throwable());

						//read input bed file, for each row,
						//String intervalFile = arguments.get(0);
						String outputPrefix = arguments.get(0);
						String vcfFile = arguments.get(1);
						
						initiate(outputPrefix);
						
						
						log.info("Parsing input vcf file ...");
						
						File indexFile = Tribble.tabixIndexFile(new File(vcfFile));
						if(!indexFile.exists() || !indexFile.canRead()){
							throw new UnsortedFileException(vcfFile + " file's index " + indexFile.getName() + " does not exist, please use tabix to index it ...");
						}

						VCFFileReader vcfReader = new VCFFileReader(new File(vcfFile), true);
						SAMSequenceDictionary seqDict = VCFFileReader.getSequenceDictionary(new File(vcfFile));
						GenomeLocParser glp = new GenomeLocParser(seqDict);
						
						CloseableIterator<VariantContext> itSnp = vcfReader.iterator();
						
						HashMap<GenomeLoc, Pair<Byte, Byte>> genotypeMapList = new HashMap<GenomeLoc, Pair<Byte, Byte>>();
						while(itSnp.hasNext()){
							VariantContext vc = itSnp.next();
							if(vc.isFiltered() || vc.isIndel() || vc.isMixed() || !vc.isPointEvent())
								continue;
							
								Genotype gt = vc.getGenotype(0);
								
								if(!gt.isFiltered() && gt.isHet()){ //only phased Het SNP is considered
									GenomeLoc genomeLoc = glp.createGenomeLoc(vc.getContig(),vc.getStart(), vc.getEnd());
									Byte alleleA = gt.getAllele(0).getBases()[0];
									Byte alleleB = gt.getAllele(1).getBases()[0];
									Pair<Byte, Byte> tmp = new Pair<Byte, Byte>(alleleA, alleleB);
									genotypeMapList.put(genomeLoc, tmp);
								}
							
							
						}
						itSnp.close();
						vcfReader.close();
						log.info("Parsing input bam files ...");
						SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
						
						HashMap<GenomeLoc,ArrayList<Double>> cpgInterEachCell = new HashMap<GenomeLoc,ArrayList<Double>>();
						
						ArrayList<Double> cpgEachCellTemp = new ArrayList<Double>(bams.size()*2);
						for(int i = 0; i < bams.size(); i++) {
							cpgEachCellTemp.add(Double.NaN); //each allele
							cpgEachCellTemp.add(Double.NaN);
						}
						for(int i = 0; i < bams.size(); i++) {
							
							String bamFile = bams.get(i);
							log.info("Parsing input bam file " +  i + " : " + bamFile );
							SamReader wgsReader = factory.open(new File(bamFile));
							for(GenomeLoc genomeLoc : genotypeMapList.keySet()) {
								String chr = genomeLoc.getContig();
								int start = genomeLoc.getStart();
								int end = genomeLoc.getStop();
								Pair<Byte, Byte> alleles = genotypeMapList.get(genomeLoc);
								SAMRecordIterator samIt = wgsReader.queryOverlapping(chr,start,end);
								while(samIt.hasNext()){
									SAMRecord r1 = samIt.next();
									if(failFlagFilter(r1)){
										continue;
									}
									int readStart = r1.getAlignmentStart();
									boolean negStrand = r1.getReadNegativeStrandFlag();
									boolean secondPair = r1.getReadPairedFlag() && r1.getSecondOfPairFlag();

									
									byte[] readBases = r1.getReadBases();
									byte[] readBasesQ = r1.getBaseQualities();
									
									
									//define which allele the reads belonged to
									
									//in C/T or A/G SNP only consider bisulfite uncoverted strand
									if((SequenceUtil.basesEqual(alleles.getFirst(), SequenceUtil.C) && SequenceUtil.basesEqual(alleles.getSecond(), SequenceUtil.T))
											|| (SequenceUtil.basesEqual(alleles.getFirst(), SequenceUtil.T) && SequenceUtil.basesEqual(alleles.getSecond(), SequenceUtil.C))
											){
										if((!secondPair && !negStrand) || (secondPair && negStrand)){
											continue;
										}
									}else if((SequenceUtil.basesEqual(alleles.getFirst(), SequenceUtil.G) && SequenceUtil.basesEqual(alleles.getSecond(), SequenceUtil.A))
											|| (SequenceUtil.basesEqual(alleles.getFirst(), SequenceUtil.A) && SequenceUtil.basesEqual(alleles.getSecond(), SequenceUtil.G))){
										if((!secondPair && negStrand) || (secondPair && !negStrand)){
											continue;
										}
									}
									
									//int pos = loc.getStart() - d.getAlignmentStart();
									int pos = r1.getReadPositionAtReferencePosition(genomeLoc.getStart())-1;
									if(pos < 0) //no such a postion in the clipped reads
										continue;
										
									byte snpBase = CcInferenceUtils.toUpperCase(readBases[pos]);
									//if(negStrand){
										//snpBase = CcInferenceUtils.toUpperCase(SequenceUtil.complement(readBases[pos]));
									//}
									if(readBasesQ[pos] < minBaseQ){
										continue;
									}
									
									TreeMap<Integer, Integer> ctSummary = BisulfiteHicUtils.readsMethySummaryAtEachLoc(r1, minMapQ, minBaseQ, turnOffBisulfiteFilter, useBadMate);
						
									if(ctSummary == null || ctSummary.isEmpty()) {
										continue;
									}
									
									if(BisulfiteHicUtils.baseMatchInBisulfiteSpace(alleles.getFirst(), snpBase, negStrand, secondPair)){
										int cellNum = i*2 + 0;
										for(int cpgLoc : ctSummary.keySet()) {
											if(pos == cpgLoc) {
												continue;
											}
											int methy = ctSummary.get(cpgLoc);
											cpgLoc = cpgLoc + readStart;
											
											GenomeLoc gCpgLoc = glp.createGenomeLoc(chr,cpgLoc, cpgLoc);
											if(cpgInterEachCell.containsKey(gCpgLoc)) {
												ArrayList<Double> cpgEachCellTemp1 = cpgInterEachCell.get(gCpgLoc);
												cpgEachCellTemp1.set(cellNum, (double)methy);//need to go back and check it later for the overlapped one..
												cpgInterEachCell.put(gCpgLoc, cpgEachCellTemp1);
											}else {
												ArrayList<Double> cpgEachCellTemp1 = (ArrayList<Double>) cpgEachCellTemp.clone();
												//System.err.println(cpgEachCellTemp.size());
												//System.err.println(cpgEachCellTemp1.size());
												//System.err.println(cellNum);
												//System.err.println(methy);
												cpgEachCellTemp1.set(cellNum, (double)methy);
												cpgInterEachCell.put(gCpgLoc, cpgEachCellTemp1);
											}
											
										}
										
									}else if(BisulfiteHicUtils.baseMatchInBisulfiteSpace(alleles.getSecond(), snpBase, negStrand, secondPair)){
										int cellNum = i*2 + 1;
										for(int cpgLoc : ctSummary.keySet()) {
											if(pos == cpgLoc) {
												continue;
											}
											int methy = ctSummary.get(cpgLoc);
											cpgLoc = cpgLoc + readStart;
											
											GenomeLoc gCpgLoc = glp.createGenomeLoc(chr,cpgLoc, cpgLoc);
											if(cpgInterEachCell.containsKey(gCpgLoc)) {
												ArrayList<Double> cpgEachCellTemp1 = cpgInterEachCell.get(gCpgLoc);
												cpgEachCellTemp1.set(cellNum, (double)methy);//need to go back and check it later for the overlapped one..
												cpgInterEachCell.put(gCpgLoc, cpgEachCellTemp1);
											}else {
												ArrayList<Double> cpgEachCellTemp1 = (ArrayList<Double>) cpgEachCellTemp.clone();
												cpgEachCellTemp1.set(cellNum, (double)methy);
												cpgInterEachCell.put(gCpgLoc, cpgEachCellTemp1);
											}
											
										}
									}
									
								}
								samIt.close();
							}
							
							
							wgsReader.close();
							
						}
						log.info("calculating pearson correlation/poisson binomial score between cpg pairs ...");
						//calculate pearson correlations:
						ArrayList<GenomeLoc> cpgList = new ArrayList<GenomeLoc>();
						cpgList.addAll(cpgInterEachCell.keySet());
						//ArrayList<GenomeLoc> cpgList = (ArrayList<GenomeLoc>) cpgInterEachCell.keySet();
						//first row, write down each column's name:
						methyCorWriter.write("chr\tpos");
						poissonBinomialWriter.write("chr\tpos");
						
						for(int i = 0; i < cpgList.size(); i++) {
							methyCorWriter.write("\t"+ i + "_" + cpgList.get(i).getStop());
							poissonBinomialWriter.write("\t"+ i + "_" + cpgList.get(i).getStop());
						}
						methyCorWriter.write("\n");
						poissonBinomialWriter.write("\n");
						for(int i = 0; i < cpgList.size(); i++) {
							GenomeLoc l1 = cpgList.get(i);
							ArrayList<Double> cpgList1 = cpgInterEachCell.get(l1);
							double[] cpg1 = ArrayUtils.toPrimitive(cpgList1.toArray(new Double[cpgList1.size()]));
							methyCorWriter.write(l1.getContig() + "\t" + l1.getStop());
							poissonBinomialWriter.write(l1.getContig() + "\t" + l1.getStop());
							for(int j = 0; j < i+1; j++) {
								methyCorWriter.write("\t" + Double.NaN);
								poissonBinomialWriter.write("\t" + Double.NaN);
							}
							for(int j = i+1; j < cpgList.size(); j++) {
								GenomeLoc l2 = cpgList.get(j);
								ArrayList<Double> cpgList2 = cpgInterEachCell.get(l2);
								double[] cpg2 = ArrayUtils.toPrimitive(cpgList2.toArray(new Double[cpgList2.size()]));
								PearsonsCorrelation p = new PearsonsCorrelation();
								
								methyCorWriter.write("\t" + p.correlation(cpg1, cpg2));
								PoissonBinomial2 pb = new PoissonBinomial2(cpg1);
								int methyEvenet2 = 0;
								for(double m : cpg2) {
									methyEvenet2 += (int)m;
								}
								if(methyEvenet2 > cpg1.length) {
									poissonBinomialWriter.write("\t" + Double.NaN);
								}else {
									poissonBinomialWriter.write("\t" + pb.get_cdf(methyEvenet2));
								}
								//poissonBinomialWriter.write("\t" + p.correlation(cpg1, cpg2));
							}
							methyCorWriter.write("\n");
							poissonBinomialWriter.write("\n");
						}
						
						finish();
									
		}
		
	
		
		private boolean failFlagFilter(SAMRecord r){
			return r.getReadUnmappedFlag() || r.getNotPrimaryAlignmentFlag() || r.getMappingQuality() < minMapQ
					|| r.getReadFailsVendorQualityCheckFlag() || r.getDuplicateReadFlag() || !r.getReadPairedFlag() || (useBadMate ? false : !r.getProperPairFlag());
		}

		
		
		
		private void initiate(String outputPrefix) throws Exception{
			startTime = System.currentTimeMillis();
			methyCorWriter = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputPrefix + ".methyCor.txt.gz")), "UTF-8");
			poissonBinomialWriter = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputPrefix + ".poissonBinomial.txt.gz")), "UTF-8");
		}
		

		private void finish() throws Exception{
			
			
				methyCorWriter.close();
				poissonBinomialWriter.close();
			
			
			
			long endTime   = System.currentTimeMillis();
			double totalTime = endTime - startTime;
			totalTime /= 1000;
			double totalTimeMins = totalTime/60;
			double totalTimeHours = totalTime/3600;
			
			log.info("MethyCorAcrossHiccupsSC's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
		}
		

}
