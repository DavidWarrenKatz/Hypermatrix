/**
 * SepAlleleMethySC.java
 * Feb 11, 2018
 * 1:03:54 PM
 * yaping    lyping1986@gmail.com
 */
package main.java.edu.mit.compbio.bisulfitehic.utils;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.exception.UnsortedFileException;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import main.java.edu.mit.compbio.ccinference.utils.CcInferenceUtils;
import org.apache.commons.math3.util.Pair;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.TreeMap;
import java.util.zip.GZIPOutputStream;


/**
 * seperate the alelle's methylation in single cell level
 * bam file is sorted by readname
 */

public class SepAlleleMethySC2 {

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

	final private static String USAGE = "SepAlleleMethySC2 [opts] outputPrefix snp_vcf.gz input.sortByReadName.bam";
	
	private static Logger log = Logger.getLogger(SepAlleleMethySC2.class);

	private OutputStreamWriter alleleAWriter = null; 
	private OutputStreamWriter alleleBWriter = null;
	
	private static long startTime = -1;

	
	/**
	 * @param args
	 */
	public static void main(String[] args)
		throws Exception {
		SepAlleleMethySC2 sams = new SepAlleleMethySC2();
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
						

						HashMap<String,IntervalTree<Pair<Byte, Byte>>> genotypeMapList = new HashMap<String,IntervalTree<Pair<Byte, Byte>>>();
						while(itSnp.hasNext()){
							VariantContext vc = itSnp.next();
							if(vc.isFiltered() || vc.isIndel() || vc.isMixed() || !vc.isPointEvent())
								continue;
							
								Genotype gt = vc.getGenotype(0);
								
								if(!gt.isFiltered() && gt.isHet()){ //only phased Het SNP is considered
									IntervalTree<Pair<Byte, Byte>> tree;
									String chr = vc.getContig();
									if(genotypeMapList.containsKey(chr)){
										tree = genotypeMapList.get(chr);
									}else{
										tree = new IntervalTree<Pair<Byte, Byte>>();
									}

									Byte alleleA = gt.getAllele(0).getBases()[0];
									Byte alleleB = gt.getAllele(1).getBases()[0];
									Pair<Byte, Byte> tmp = new Pair<Byte, Byte>(alleleA, alleleB);
									tree.put(vc.getStart(), vc.getEnd(), tmp);
									genotypeMapList.put(chr, tree);
								}
							
							
						}
						itSnp.close();
						vcfReader.close();
						log.info("Parsing input bam files ...");
							SamReader wgsReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamFile));
								//System.err.println(chr + "\t" + start + "\t" + end + "\n");
								SAMRecordIterator samIt = wgsReader.iterator();
								while(samIt.hasNext()){
									SAMRecord r1 = samIt.next();
									SAMRecord r2 = samIt.next();
									if(failFlagFilter(r1) || failFlagFilter(r2)){
										continue;
									}
									String chr1 = r1.getContig();
									String chr2 = r2.getContig();
									if(!chr1.equalsIgnoreCase(chr2)) {
										continue;
									}
									int start1 = r1.getAlignmentStart();
									int end1 = r1.getAlignmentEnd();
									int start2 = r2.getAlignmentStart();
									int end2 = r2.getAlignmentEnd();
									Pair<Byte, Byte> alleles = null;
									int snpStart=-1;
									//int snpEnd;
									boolean snpInRead2 = false;
									if(genotypeMapList.containsKey(chr1)){
										IntervalTree<Pair<Byte, Byte>> cg = genotypeMapList.get(chr1);

										if(cg.minOverlapper(start1, end1)==null && cg.minOverlapper(start2, end2)==null){
											continue;
										}else{
											if(cg.minOverlapper(start1, end1)!=null){
												Node<Pair<Byte, Byte>> node = cg.minOverlapper(start1, end1);
												snpStart = node.getStart();
												//snpEnd = node.getEnd();
												alleles = node.getValue();
												snpInRead2 = false;
											}else if(cg.minOverlapper(start2, end2)!=null){
												Node<Pair<Byte, Byte>> node = cg.minOverlapper(start2, end2);
												snpStart = node.getStart();
												//snpEnd = node.getEnd();
												alleles = node.getValue();
												snpInRead2 = true;
											}

										}
									}else{
										continue;
									}

									int read1Start = r1.getAlignmentStart();
									int read2Start = r2.getAlignmentStart();

									boolean negStrand = snpInRead2 ? r2.getReadNegativeStrandFlag() : r1.getReadNegativeStrandFlag();
									boolean secondPair = snpInRead2 ? (r2.getReadPairedFlag() && r2.getSecondOfPairFlag()) : (r1.getReadPairedFlag() && r1.getSecondOfPairFlag());


									byte[] readBases = snpInRead2 ? r2.getReadBases() : r1.getReadBases();
									byte[] readBasesQ = snpInRead2 ? r2.getBaseQualities() : r1.getBaseQualities();


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
									int pos = snpInRead2 ? (r2.getReadPositionAtReferencePosition(snpStart)-1) : (r1.getReadPositionAtReferencePosition(snpStart)-1);
									if(pos < 0) //no such a postion in the clipped reads
										continue;

									byte snpBase = CcInferenceUtils.toUpperCase(readBases[pos]);
									//if(negStrand){
										//snpBase = CcInferenceUtils.toUpperCase(SequenceUtil.complement(readBases[pos]));
									//}
									if(readBasesQ[pos] < minBaseQ){
										continue;
									}

									TreeMap<Integer, Integer> ctSummary1 = BisulfiteHicUtils.readsMethySummaryAtEachLoc(r1, minMapQ, minBaseQ, turnOffBisulfiteFilter, useBadMate);
									TreeMap<Integer, Integer> ctSummary2 = BisulfiteHicUtils.readsMethySummaryAtEachLoc(r2, minMapQ, minBaseQ, turnOffBisulfiteFilter, useBadMate);

									if((ctSummary1 == null || ctSummary1.isEmpty()) && (ctSummary2 == null || ctSummary2.isEmpty())) {
										continue;
									}

									if(BisulfiteHicUtils.baseMatchInBisulfiteSpace(alleles.getFirst(), snpBase, negStrand, secondPair)){
										for(int cpgLoc : ctSummary1.keySet()) {
											if(pos == cpgLoc) {
												continue;
											}
											int methy = ctSummary1.get(cpgLoc);
											int cpg_end = cpgLoc + read1Start;
											int cpg_start = cpg_end - 1;
											alleleAWriter.write(chr1 + "\t" + cpg_start + "\t" + cpg_end + "\t" + methy + "\n");
										}
										for(int cpgLoc : ctSummary2.keySet()) {
											if(pos == cpgLoc) {
												continue;
											}
											int methy = ctSummary2.get(cpgLoc);
											int cpg_end = cpgLoc + read2Start;
											int cpg_start = cpg_end - 1;
											alleleAWriter.write(chr1 + "\t" + cpg_start + "\t" + cpg_end + "\t" + methy + "\n");
										}
									}else if(BisulfiteHicUtils.baseMatchInBisulfiteSpace(alleles.getSecond(), snpBase, negStrand, secondPair)){
										for(int cpgLoc : ctSummary1.keySet()) {
											if(pos == cpgLoc) {
												continue;
											}
											int methy = ctSummary1.get(cpgLoc);
											int cpg_end = cpgLoc + read1Start;
											int cpg_start = cpg_end - 1;
											alleleBWriter.write(chr1 + "\t" + cpg_start + "\t" + cpg_end + "\t" + methy + "\n");
										}
										for(int cpgLoc : ctSummary2.keySet()) {
											if(pos == cpgLoc) {
												continue;
											}
											int methy = ctSummary2.get(cpgLoc);
											int cpg_end = cpgLoc + read2Start;
											int cpg_start = cpg_end - 1;
											alleleBWriter.write(chr1 + "\t" + cpg_start + "\t" + cpg_end + "\t" + methy + "\n");
										}
									}
									
								}
								samIt.close();
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
			
			log.info("SepAlleleMethySC2's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
		}
		

}
