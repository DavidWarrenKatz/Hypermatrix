/**
 * SepAlleleMethySC.java
 * Feb 11, 2018
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
 * seperate the alelle's methylation in single cell level
 */

public class SepAlleleMethySC {

	@Option(name="-minBaseQ",usage="minimum base quality score required to check. Default: 5")
	public int minBaseQ = 5;	
	
	@Option(name="-useBadMate",usage="use not proper paired reads. Default: false")
	public boolean useBadMate = false;
	
	@Option(name="-minMapQ",usage="minimum mapping quality score required to check. Default: 30")
	public int minMapQ = 30;

	@Option(name="-turnOffBisulfiteFilter",usage="turn off bisulfite conversion reads filter. Default: false")
	public boolean turnOffBisulfiteFilter = false;

	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "SepAlleleMethySC [opts] outputPrefix snp_vcf.gz input.bam";
	
	private static Logger log = Logger.getLogger(SepAlleleMethySC.class);

	private OutputStreamWriter alleleAWriter = null; 
	private OutputStreamWriter alleleBWriter = null;
	
	private static long startTime = -1;

	
	/**
	 * @param args
	 */
	public static void main(String[] args)
		throws Exception {
		SepAlleleMethySC sams = new SepAlleleMethySC();
			BasicConfigurator.configure();
			sams.doMain(args);
		}
		

		public void doMain(String[] args)
				throws Exception {

						CmdLineParser parser = new CmdLineParser(this);
						//parser.setUsageWidth(80);
						try
						{
							if(help || args.length < 3) throw new CmdLineException(parser, USAGE, new Throwable());
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
						

						//read input bed file, for each row,
						//String intervalFile = arguments.get(0);
						String outputPrefix = arguments.get(0);
						String vcfFile = arguments.get(1);
						String bamFile = arguments.get(2);
						
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
						
							log.info("Parsing input bam file " );
							SamReader wgsReader = factory.open(new File(bamFile));
							for(GenomeLoc genomeLoc : genotypeMapList.keySet()) {
								String chr = genomeLoc.getContig();
								int start = genomeLoc.getStart();
								int end = genomeLoc.getStop();
								Pair<Byte, Byte> alleles = genotypeMapList.get(genomeLoc);
								//System.err.println(chr + "\t" + start + "\t" + end + "\n");
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
										for(int cpgLoc : ctSummary.keySet()) {
											if(pos == cpgLoc) {
												continue;
											}
											int methy = ctSummary.get(cpgLoc);
											int cpg_end = cpgLoc + readStart;
											int cpg_start = cpg_end - 1;
											alleleAWriter.write(chr + "\t" + cpg_start + "\t" + cpg_end + "\t" + methy + "\n");
										}
										
									}else if(BisulfiteHicUtils.baseMatchInBisulfiteSpace(alleles.getSecond(), snpBase, negStrand, secondPair)){
										for(int cpgLoc : ctSummary.keySet()) {
											if(pos == cpgLoc) {
												continue;
											}
											int methy = ctSummary.get(cpgLoc);
											int cpg_end = cpgLoc + readStart;
											int cpg_start = cpg_end - 1;
											alleleBWriter.write(chr + "\t" + cpg_start + "\t" + cpg_end + "\t" + methy + "\n");
											
										}
									}
									
								}
								samIt.close();
							}
							wgsReader.close();
							
						
						
						
						finish();
									
		}
		
	
		
		private boolean failFlagFilter(SAMRecord r){
			return r.getReadUnmappedFlag() || r.getNotPrimaryAlignmentFlag() || r.getMappingQuality() < minMapQ
					|| r.getReadFailsVendorQualityCheckFlag() || r.getDuplicateReadFlag() || !r.getReadPairedFlag() || (useBadMate ? false : !r.getProperPairFlag());
		}

		
		
		
		private void initiate(String outputPrefix) throws Exception{
			startTime = System.currentTimeMillis();
			alleleAWriter = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputPrefix + ".alleleA_methy.bedgraph.gz")), "UTF-8");
			alleleBWriter = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputPrefix + ".alleleB_methy.bedgraph.gz")), "UTF-8");
		}
		

		private void finish() throws Exception{
			
			
			alleleAWriter.close();
			alleleBWriter.close();
			
			
			
			long endTime   = System.currentTimeMillis();
			double totalTime = endTime - startTime;
			totalTime /= 1000;
			double totalTimeMins = totalTime/60;
			double totalTimeHours = totalTime/3600;
			
			log.info("SepAlleleMethySC's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
		}
		

}
