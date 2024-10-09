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
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.exception.UnsortedFileException;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.*;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import main.java.edu.mit.compbio.ccinference.utils.CcInferenceUtils;

import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.stat.inference.ChiSquareTest;
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
 * iterate each reads that two ends are located at two regions, each fragment is a datapoint in the methy_cor output txt.gz file. (chr1,read.start1,read.end1,chr2,read.start2,read.end2,region_name, id, hiccups_ob, hiccups_exp, hiccups_fdr, methy_end1, cpg_num_end1, methy_end2, cpg_num_end2) 
 * iterate each reads that overlapped with heaterozygous SNPs in VCF file (also have another end mapped to the other part of hiccups.) summaryize the methy ratio on two different allele in two end
 * SNP.chr, SNP.start,SNP.end,region_name, id,hiccups_ob, hiccups_exp, hiccups_fdr, SNP_end.alleleA_all.methy_count,SNP_end.alleleA_all.unmethy_count,SNP_end.alleleB_all.methy_count,SNP_end.alleleB_all.unmethy_count,methyDiff(A-B),fisher_test_p_value, FDR, 
 * SNP_end.alleleA.methy_count,SNP_end.alleleA.unmethy_count,SNP_end.alleleB.methy_count,SNP_end.alleleB.unmethy_count,methyDiff(A-B),fisher_test_p_value, FDR,
 * nonSNP_end.alleleA.methy_count,nonSNP_end.alleleA.unmethy_count,nonSNP_end.alleleB.methy_count,nonSNP_end.alleleB.unmethy_count,methyDiff(A-B), fisher_test_p_value, FDR
 * (maybe further consider the other end have SNP also or not...)
 * then overlapped with GWAS other SNP loci...
 * currently only work for hiccups within the same chromosomes
 * 
 * Input: bam file should be sorted, calmd, and read name for end1 and end2 should be the same
 * TODO: after get correct calmd bam file, need to change complementArray(bases) in negative strand case
 */

public class MethyCorAcrossHiccups {

	@Option(name="-minBaseQ",usage="minimum base quality score required to check. Default: 5")
	public int minBaseQ = 5;
	
	@Option(name="-turnOffBisulfiteFilter",usage="turn off bisulfite conversion reads filter. Default: false")
	public boolean turnOffBisulfiteFilter = false;
	
	@Option(name="-useBadMate",usage="use not proper paired reads. Default: false")
	public boolean useBadMate = false;

	@Option(name="-useUnpaired",usage="reads are mapped in single end mode, have not joined into pairs yet. Default: false")
	public boolean useUnpaired = false;

	@Option(name="-removeReadnameSuffix",usage="remove the last two characters in read name, such as /1, /2 .... Default: false")
	public boolean removeReadnameSuffix = false;
	
	@Option(name="-minMapQ",usage="minimum mapping quality score required to check. Default: 30")
	public int minMapQ = 30;

	@Option(name="-minDistEnds",usage="minimum distance between two end to be considered as HiC links. Default: 20000")
	public int minDistEnds = 20000;
	
	@Option(name="-vcf",usage="vcf file. Default: null")
	public String vcf = null;
	
	@Option(name="-coverageRef",usage="For ASM mode, specify the minimum  coverage in SNP's allele Ref, default: 2")
	public int coverageRef = 2;

	@Option(name="-coverageAlt",usage="For ASM mode, specify the minimum coverage in SNP's allele Alt, default: 2")
	public int coverageAlt = 2;
	
	@Option(name="-minCgNum",usage="For ASM mode, specify the minimum  number of CpG to call ASM. Default: 1")
	public int minCgNum = 1;
	
	@Option(name="-sample",usage="specify the sample number in the VCF file if contain multiple sample genotype. number started at 0. default: 0")
	public int sample = 0;

	@Option(name="-randomHiccupsInput",usage="when the input file is random hiccups links for the control, no need to restrict each end of fragment at the same hiccups links.., default: not enabled")
	public boolean randomHiccupsInput = false;

	@Option(name="-randomNotLinked",usage="randomly sample the reads from the regions specified, the random reads do not need to be linked with each other in long-range, good for low-coverage data. default: not enabled")
	public boolean randomNotLinked = false;

	@Option(name="-shortLowerBound",usage="only use the paired end within certain short distance. this is the lower bound of distance, use it together with -onlyShort. default: -1")
	public int shortLowerBound = -1;

	@Option(name="-onlyShort",usage="only use the paired end within very short region. when enable it with positive number, it will automate disable -minDistEnds option, default: -1")
	public int onlyShort = -1;

	@Option(name="-randomSampleReads",usage="only use part of the reads in each interval. good for -randomHiccupsInput mode or -onlyShort mode, default: -1. not enabled")
	public int randomSampleReads = -1;

	@Option(name="-minReads",usage="minimum pairs of reads within the regions. default: 0")
	public int minReads = 0;
	
	@Option(name="-skipProcessFirstRow",usage="don't process first row, default: not enabled")
	public boolean skipProcessFirstRow = false;

	@Option(name="-diffChr",usage="don't skip the region from two different chromosomes, default: not enabled")
	public boolean diffChr = false;

	@Option(name="-useGeneralBedPe",usage="use general bedpe file rather than HiCUUPS link file, default: not enabled")
	public boolean useGeneralBedPe = false;

	@Option(name="-outputCount",usage="output the raw read count frequency between two end when not in ASM mode, default: not enabled")
	public boolean outputCount = false;

	@Option(name="-methyPatternSearch",usage="methylation pattern to search, allow IUPAC match code, such as HCG or GCH for NOMe-seq, default: CG")
	public String methyPatternSearch = "CG";

	@Option(name="-methyPatternCytPos",usage="cytosine position in -methyPatternSearch. for example, 1 for CG and 2 for HCG, default: 1")
	public int methyPatternCytPos = 1;

	@Option(name="-bsConvPatternSearch",usage="bisulfite conversion pattern to search, allow IUPAC match code, such as WCH for NOMe-seq, default: CH")
	public String bsConvPatternSearch = "CH";

	@Option(name="-bsConvPatternCytPos",usage="cytosine position in -bsConvPatternSearch. for example, 1 for CH and 2 for WCH, default: 1")
	public int bsConvPatternCytPos = 1;

	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "MethyCorAcrossHiccups [opts] outputPrefix hiccups.bedpe.bed[.gz] input.bam";
	
	private static Logger log = Logger.getLogger(MethyCorAcrossHiccups.class);

	private OutputStreamWriter methyCorWriter = null; 
	private OutputStreamWriter asmWriter = null; 
	private MersenneTwister generator = null;
	
	private VCFFileReader vcfReader = null;
	private GenomeLocParser glp = null;
	
	private static long startTime = -1;
	private static long lineNum=0;
	private static long queryNum=0;
	private static long totalReadPairNum=0;
	/**
	 * @param args
	 */
	public static void main(String[] args)
		throws Exception {
		MethyCorAcrossHiccups mcah = new MethyCorAcrossHiccups();
			BasicConfigurator.configure();
			mcah.doMain(args);
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
						//String refFile = arguments.get(1);
						String hiccupsFile = arguments.get(1);
						String bamFile = arguments.get(2);

						initiate(outputPrefix);
						//SamReader wgsReader = SamReaderFactory.makeDefault().setOption(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES, true).validationStringency(ValidationStringency.SILENT).open(new File(bamFile));
						SamReader wgsReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamFile));
						log.info("Parsing input hiccups links BEDPE format bed file ...");
						GZIPInputStream gzipInputStream = null;
						BufferedReader br;
						if(hiccupsFile.endsWith(".gz")){
							gzipInputStream = new GZIPInputStream(new FileInputStream(hiccupsFile));
							br = new BufferedReader(new InputStreamReader(gzipInputStream));
							
						}else{
							br = new BufferedReader(new FileReader(hiccupsFile));
						}
						String line;
						
						while( (line = br.readLine()) != null){
							if((skipProcessFirstRow && lineNum==0) || line.startsWith("#")){
								lineNum++;
								continue;
							}else{
								String[] splitin = line.split("\t");
								String chr1 = splitin[0];
								int start1 = Integer.parseInt(splitin[1]);
								int end1 = Integer.parseInt(splitin[2]);
								String chr2 = splitin[3];
								int start2 = Integer.parseInt(splitin[4]);
								int end2 = Integer.parseInt(splitin[5]);
								if((!diffChr && !chr1.equalsIgnoreCase(chr2)) || (onlyShort<0 && Math.abs(start1-start2)<minDistEnds) || (onlyShort<0 && Math.abs(end1-end2)<minDistEnds)){
									continue;
								}
								double ob = Double.NaN;
								double exp = Double.NaN;
								double fdr = Double.NaN;
								if(!useGeneralBedPe){
									ob = Double.parseDouble(splitin[7]);
									exp = Double.parseDouble(splitin[11]);
									fdr = Double.parseDouble(splitin[15]);
								}
								
								if(vcf != null){
									asmCall(chr1, start1, end1, start2, end2, wgsReader, ob, exp, fdr);
								}else{
									if(randomHiccupsInput){
										//log.info("randomHiccupsInput mode: ");
										methyCorrelationRandom(chr1, start1, end1, chr2, start2, end2, wgsReader, ob, exp, fdr, bamFile);
									}else{
										if(onlyShort>0){
											//log.info("only use paired end reads within short distance: ");
											methyCorrelationShort(chr1, start1, end1, start2, end2, wgsReader, ob, exp, fdr, bamFile);
										}else{
											methyCorrelation(chr1, start1, end1, chr2, start2, end2, wgsReader, ob, exp, fdr, bamFile);
										}
										
									}
								}
								
								
								
							}
								
							
							lineNum++;
							if(lineNum % 1000 == 0){
								log.info("Processing line: " + lineNum);
								if(vcf != null){
									asmWriter.flush();
								}else{
									methyCorWriter.flush();
								}
							}
							
						}
						if(hiccupsFile.endsWith(".gz")){
							gzipInputStream.close();
						}
						br.close();
						
						finish();
		}
		
		//currently, not consider condition that each end overlapped a SNP
		//TODO
		private void asmCall(String chr, int start1, int end1, int start2, int end2, SamReader wgsReader, double ob, double exp, double fdr) throws Exception{ 
			String regionName = chr + ":" + start1 + ":" + end1 + ":" + start2 + ":" + end2;
			
			//read SNP list
			//region1
			CloseableIterator<VariantContext> itSnp = vcfReader.query(chr, start1, end1);
			
			HashMap<GenomeLoc, Pair<Byte, Byte>> genotypeMapList1 = new HashMap<GenomeLoc, Pair<Byte, Byte>>();
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
						genotypeMapList1.put(genomeLoc, tmp);
					}
				
				
			}
			itSnp.close();
			
			//region2
			itSnp = vcfReader.query(chr, start2, end2);
			HashMap<GenomeLoc, Pair<Byte, Byte>> genotypeMapList2 = new HashMap<GenomeLoc, Pair<Byte, Byte>>();
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
						genotypeMapList2.put(genomeLoc, tmp);
					}
				
				
			}
			itSnp.close();
			
			//read all reads within the region
			SAMRecordIterator wgsIt = wgsReader.queryOverlapping(chr,start1,end1);
			HashMap<String, SAMRecord> countedReadsEnd1 = new HashMap<String, SAMRecord>();	
			while(wgsIt.hasNext()){
				SAMRecord r = wgsIt.next();
				
				if(failFlagFilter(r)){
					continue;
				}else{
					
					if(countedReadsEnd1.containsKey(r.getReadName())){
						countedReadsEnd1.remove(r.getReadName()); //if two ends are within the same half end of hiccup links
					}else{
						countedReadsEnd1.put(r.getReadName(), r);
					}
				}
				
			}
			wgsIt.close();
			
			HashMap<String, SAMRecord> countedReadsEnd2 = new HashMap<String, SAMRecord>();
			wgsIt = wgsReader.queryOverlapping(chr,start2,end2);
			while(wgsIt.hasNext()){
				SAMRecord r = wgsIt.next();
				if(failFlagFilter(r)){
					continue;
				}else{
					if(countedReadsEnd2.containsKey(r.getReadName())){
						countedReadsEnd2.remove(r.getReadName()); //if two ends are within the same half end of hiccup links
					}else{
						countedReadsEnd2.put(r.getReadName(), r);
					}
				}
			}
			wgsIt.close();
			//end1's SNP
			processBam(genotypeMapList1, countedReadsEnd1, countedReadsEnd2, wgsReader, regionName, ob, exp, fdr);
			
			//end2's SNP
			processBam(genotypeMapList2, countedReadsEnd2, countedReadsEnd1, wgsReader, regionName, ob, exp, fdr);
			
		}
		
		private void processBam(HashMap<GenomeLoc, Pair<Byte, Byte>> genotypeMapList1,  HashMap<String, SAMRecord> countedReadsEnd1, 
				HashMap<String, SAMRecord> countedReadsEnd2, SamReader wgsReader, String regionName, double ob, double exp, double fdr) throws Exception{
			//end1
			int snpNum = 0;
			for(GenomeLoc loc : genotypeMapList1.keySet()){
				Pair<Byte, Byte> alleles = genotypeMapList1.get(loc);
				SAMRecordIterator samIt = wgsReader.queryOverlapping(loc.getContig(), loc.getStart(), loc.getStop());

				Pair<Integer, Integer> ctAlleleRefAll = new Pair<Integer, Integer>(0,0);
				Pair<Integer, Integer> ctAlleleAltAll = new Pair<Integer, Integer>(0,0);
				Pair<Integer, Integer> ctAlleleRefEnd1 = new Pair<Integer, Integer>(0,0);
				Pair<Integer, Integer> ctAlleleAltEnd1 = new Pair<Integer, Integer>(0,0);
				Pair<Integer, Integer> ctAlleleRefEnd2 = new Pair<Integer, Integer>(0,0);
				Pair<Integer, Integer> ctAlleleAltEnd2 = new Pair<Integer, Integer>(0,0);
				
				int refAlleleCount = 0;
				int altAlleleCount = 0;


				while(samIt.hasNext()){
					SAMRecord r1 = samIt.next();
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
					int pos = r1.getReadPositionAtReferencePosition(loc.getStart())-1;
					if(pos < 0) //no such a postion in the clipped reads
						continue;
						
					byte snpBase = CcInferenceUtils.toUpperCase(readBases[pos]);
					//if(negStrand){
						//snpBase = CcInferenceUtils.toUpperCase(SequenceUtil.complement(readBases[pos]));
					//}
					if(readBasesQ[pos] < minBaseQ){
						continue;
					}

					//System.err.println((char)snpBase + "\t" + pos + "\t" + (char)alleles.getFirst().byteValue() + "\t" + (char)alleles.getSecond().byteValue()+ "\t" + negStrand + "\t" + secondPair + "\t" +
					//BisulfiteHicUtils.baseMatchInBisulfiteSpace(alleles.getFirst(), snpBase, negStrand, secondPair) + "\t" + 
					//new String(CcInferenceUtils.toUpperCase(readBases)));
					
					
					Triple<Integer,Integer, Integer> ctSummary = BisulfiteHicUtils.readsMethySummaryGeneral(r1, methyPatternSearch, methyPatternCytPos, bsConvPatternSearch, bsConvPatternCytPos, minMapQ, minBaseQ, turnOffBisulfiteFilter, useBadMate);
					int numC = ctSummary.getLeft();
					int numT = ctSummary.getMiddle();
					int numCg = ctSummary.getRight();
					if((numC + numT) == 0 || numCg < minCgNum)
						continue;
					
					Triple<Integer,Integer, Integer> ctSummary2 = Triple.of(0,0,0);
					boolean matchedEndFound = false;
					if(countedReadsEnd2.containsKey(r1.getReadName())){
						matchedEndFound = true;
						ctSummary2 = BisulfiteHicUtils.readsMethySummaryGeneral(countedReadsEnd2.get(r1.getReadName()), methyPatternSearch, methyPatternCytPos, bsConvPatternSearch, bsConvPatternCytPos, minMapQ, minBaseQ, turnOffBisulfiteFilter, useBadMate);
					}else{
						
					}
							
					if(BisulfiteHicUtils.baseMatchInBisulfiteSpace(alleles.getFirst(), snpBase, negStrand, secondPair)){
						int numCRef = ctAlleleRefAll.getFirst() + numC;
						int numTRef = ctAlleleRefAll.getSecond() + numT;
						ctAlleleRefAll = new Pair<Integer, Integer>(numCRef,numTRef);
						refAlleleCount++;
						
						if(matchedEndFound){
							numCRef = ctAlleleRefEnd1.getFirst() + ctSummary.getLeft();
							numTRef = ctAlleleRefEnd1.getSecond() + ctSummary.getMiddle();
							ctAlleleRefEnd1 = new Pair<Integer, Integer>(numCRef,numTRef);
							
							numCRef = ctAlleleRefEnd2.getFirst() + ctSummary2.getLeft();
							numTRef = ctAlleleRefEnd2.getSecond() + ctSummary2.getMiddle();
							ctAlleleRefEnd2 = new Pair<Integer, Integer>(numCRef,numTRef);
						}
						
					}else if(BisulfiteHicUtils.baseMatchInBisulfiteSpace(alleles.getSecond(), snpBase, negStrand, secondPair)){
						int numCAlt = ctAlleleAltAll.getFirst() + numC;
						int numTAlt = ctAlleleAltAll.getSecond() + numT;
						ctAlleleAltAll = new Pair<Integer, Integer>(numCAlt,numTAlt);
						altAlleleCount++;
						
						if(matchedEndFound){
							numCAlt = ctAlleleAltEnd1.getFirst() + ctSummary.getLeft();
							numTAlt = ctAlleleAltEnd1.getSecond() + ctSummary.getMiddle();
							ctAlleleAltEnd1 = new Pair<Integer, Integer>(numCAlt,numTAlt);
							
							numCAlt = ctAlleleAltEnd2.getFirst() + ctSummary2.getLeft();
							numTAlt = ctAlleleAltEnd2.getSecond() + ctSummary2.getMiddle();
							ctAlleleAltEnd2 = new Pair<Integer, Integer>(numCAlt,numTAlt);
						}
					}
					
					
					
				}
				samIt.close();
				if(refAlleleCount < coverageRef || altAlleleCount< coverageAlt || (ctAlleleRefAll.getFirst() + ctAlleleRefAll.getSecond() ==0 && ctAlleleAltAll.getFirst() + ctAlleleAltAll.getSecond() ==0)
						|| (ctAlleleRefEnd1.getFirst() + ctAlleleRefEnd1.getSecond() ==0 && ctAlleleAltEnd1.getFirst() + ctAlleleAltEnd1.getSecond() ==0) ||
						(ctAlleleRefEnd2.getFirst() + ctAlleleRefEnd2.getSecond() ==0 && ctAlleleAltEnd2.getFirst() + ctAlleleAltEnd2.getSecond() ==0))
					continue;
				
				/* SNP.chr, SNP.start,SNP.end,SNP basesA, SNP baseB, region_name, id,hiccups_ob, hiccups_exp, hiccups_fdr, SNP_end.alleleA_all.methy_count,SNP_end.alleleA_all.unmethy_count,SNP_end.alleleB_all.methy_count,SNP_end.alleleB_all.unmethy_count,methyDiff(B-A),fisher_test_p_value, 
				 * SNP_end.alleleA.methy_count,SNP_end.alleleA.unmethy_count,SNP_end.alleleB.methy_count,SNP_end.alleleB.unmethy_count,methyDiff(B-A),fisher_test_p_value,
				 * nonSNP_end.alleleA.methy_count,nonSNP_end.alleleA.unmethy_count,nonSNP_end.alleleB.methy_count,nonSNP_end.alleleB.unmethy_count,methyDiff(B-A), fisher_test_p_value
				 */
				
				
				Pair<Double, Double> summaryStatAll = getPvalueAndDiff(ctAlleleRefAll.getFirst(), ctAlleleRefAll.getSecond(), ctAlleleAltAll.getFirst(), ctAlleleAltAll.getSecond());
				Pair<Double, Double> summaryStatEnd1 = getPvalueAndDiff(ctAlleleRefEnd1.getFirst(), ctAlleleRefEnd1.getSecond(), ctAlleleAltEnd1.getFirst(), ctAlleleAltEnd1.getSecond());
				Pair<Double, Double> summaryStatEnd2 = getPvalueAndDiff(ctAlleleRefEnd2.getFirst(), ctAlleleRefEnd2.getSecond(), ctAlleleAltEnd2.getFirst(), ctAlleleAltEnd2.getSecond());
				
				
				asmWriter.write(loc.getContig() + "\t" + (loc.getStart()-1)  + "\t" + loc.getStop()  + "\t" + String.format("%c", alleles.getFirst()) + "\t" + String.format("%c", alleles.getSecond()) + "\t"
						+ regionName + "\t" + snpNum + "\t" + ob + "\t" + exp + "\t" + fdr  + "\t" + ctAlleleRefAll.getFirst() + "\t" + ctAlleleRefAll.getSecond()  + "\t" + ctAlleleAltAll.getFirst() + "\t" + ctAlleleAltAll.getSecond() + "\t" + summaryStatAll.getFirst() + "\t" + summaryStatAll.getSecond()
						+ "\t" + ctAlleleRefEnd1.getFirst() + "\t" + ctAlleleRefEnd1.getSecond()+ "\t" + ctAlleleAltEnd1.getFirst() + "\t" + ctAlleleAltEnd1.getSecond() + "\t" + summaryStatEnd1.getFirst() + "\t" + summaryStatEnd1.getSecond()
						+ "\t" + ctAlleleRefEnd2.getFirst() + "\t" + ctAlleleRefEnd2.getSecond()+ "\t" + ctAlleleAltEnd2.getFirst() + "\t" + ctAlleleAltEnd2.getSecond() + "\t" + summaryStatEnd2.getFirst() + "\t" + summaryStatEnd2.getSecond() + "\n");
				snpNum++;
			}
			
		}
		
		//first is snpB_methy_level-snpA_methy_level , second is p value
		private Pair<Double, Double> getPvalueAndDiff(int snpA_methy, int snpA_unmethy, int snpB_methy, int snpB_unmethy){
			double methyDiff = Double.NaN;
			if((snpA_methy + snpA_unmethy > 0) && (snpB_methy + snpB_unmethy > 0) ){
				methyDiff = (double)snpB_methy/(double)(snpB_methy + snpB_unmethy) - (double)snpA_methy/(double)(snpA_methy + snpA_unmethy);
				methyDiff*=100;
			}
			double pvalue = Double.NaN;
			if(snpA_methy + snpA_unmethy + snpB_methy + snpB_unmethy > 10000){ //use ChiSquareTest
				ChiSquareTest test = new ChiSquareTest();
				long[][] input = new long[][]{{snpA_methy, snpA_unmethy},{snpB_methy, snpB_unmethy}};
				pvalue = test.chiSquareTest(input);
			}else{
				FisherExactTest test = new FisherExactTest(10000);
				pvalue = test.getTwoTailedP(snpA_methy, snpA_unmethy, snpB_methy, snpB_unmethy);
			}
			return new Pair<Double, Double>(methyDiff, pvalue);
		}
		
		private void methyCorrelation(String chr, int start1, int end1,String chr2, int start2, int end2, SamReader wgsReader, double ob, double exp, double fdr, String bamFile) throws Exception{
			
			if(outputCount) {
				SAMRecordIterator wgsIt = wgsReader.iterator();
				while(wgsIt.hasNext()){
					SAMRecord r = wgsIt.next();
					if(failFlagFilter(r)){
						continue;
					}else{
						totalReadPairNum++;
					}
				}
				wgsIt.close();
				totalReadPairNum/=2;
			}
			
			String regionName = chr + ":" + start1 + ":" + end1 + ":" + chr2 + ":" + start2 + ":" + end2;
			
			SAMRecordIterator wgsIt = wgsReader.queryOverlapping(chr,start1,end1);
			HashMap<String, SAMRecord> countedReadsEnd1 = new HashMap<String, SAMRecord>();	
			boolean reopenFlag=false;
			while(wgsIt.hasNext()){
				SAMRecord r = wgsIt.next();
				queryNum++;
				if(queryNum % 1000000 == 0){
					reopenFlag=true;
				}
			//	if(r.getReadName().equalsIgnoreCase("K00168:28:H3NKLBBXX:7:1214:10064:40385")){
					//refParser.setCurrentSequence(chr);
					//byte[] refBaseEnd1 = refParser.loadFragment(r.getAlignmentStart()-1, r.getAlignmentEnd()-r.getAlignmentStart()+1).getBytes();
					//if(r.getReadNegativeStrandFlag()){
					//	refBaseEnd1=BisulfiteHicUtils.complementArray(refBaseEnd1);
					//}
					//byte[] readBases = BisulfiteHicUtils.getClippedReadsBase(r);
					//int pos = CcInferenceUtils.bisulfiteIncompleteReads(readBases, refBaseEnd1, r.getReadNegativeStrandFlag(), r.getReadPairedFlag() & r.getSecondOfPairFlag(), false, true, false, 0.1, 0.4, 0, "CH", 1);
					//int truncateStart = r.getReadPositionAtReferencePosition(r.getAlignmentStart());
					//int truncateEnd = r.getReadPositionAtReferencePosition(r.getAlignmentEnd());
					
					//System.err.println(truncateStart + "\t" + truncateEnd );
					//System.err.println(new String(refBaseEnd1));
					//System.err.println(new String(readBases));
					//System.err.println(r.getReadName() + "\t" + countedReadsEnd1.containsKey(r.getReadName()) + "\t" + r + "\t" + pos);
				//}
				if(failFlagFilter(r)){
					continue;
				}else{
					//int bisulfiteStartPos = CcInferenceUtils.bisulfiteIncompleteReads(r);
					//if(bisulfiteStartPos<0){
					//	continue;
					//}
					//if(r.getReadName().equalsIgnoreCase("K00168:28:H3NKLBBXX:7:1214:10064:40385")){
					//	refParser.setCurrentSequence(chr);
					//	byte[] refBaseEnd1 = refParser.loadFragment(r.getAlignmentStart()-1, r.getAlignmentEnd()-r.getAlignmentStart()+1).getBytes();
					//	int pos = CcInferenceUtils.bisulfiteIncompleteReads(r.getReadBases(), refBaseEnd1, r.getReadNegativeStrandFlag(), r.getReadPairedFlag() & r.getSecondOfPairFlag(), false, true, false, 0.1, 0.4, 0, "CH", 1);
						
					//	System.err.println(new String(refBaseEnd1));
					//	System.err.println(new String(r.getReadBases()));
					//	System.err.println(r.getReadName() + "\t" + countedReadsEnd1.containsKey(r.getReadName()) + "\t" + r + "\t" + pos);
					//}
					String name = r.getReadName();
					if(removeReadnameSuffix){
						name = name.substring(0, name.length() - 2);
					}
					if(countedReadsEnd1.containsKey(name)){
						countedReadsEnd1.remove(name); //if two ends are within the same half end of hiccup links
					}else{
						countedReadsEnd1.put(name, r);
					}
				}
				
			}
			wgsIt.close();
			if(reopenFlag){
				//wgsReader.close();
				//wgsReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamFile));
				reopenFlag=false;
			}
			HashMap<String, SAMRecord> countedReadsEnd2 = new HashMap<String, SAMRecord>();
			wgsIt = wgsReader.queryOverlapping(chr2,start2,end2);
			while(wgsIt.hasNext()){
				SAMRecord r = wgsIt.next();
				queryNum++;
				if(queryNum % 1000000 == 0){
					reopenFlag=true;
				}
				//if(r.getReadName().equalsIgnoreCase("K00168:28:H3NKLBBXX:7:1214:10064:40385")){
					//byte[] refBaseEnd1 = refParser.loadFragment(r.getAlignmentStart()-1, r.getAlignmentEnd()-r.getAlignmentStart()+1).getBytes();
					//if(r.getReadNegativeStrandFlag()){
					//	refBaseEnd1=BisulfiteHicUtils.complementArray(refBaseEnd1);
					//}
					//byte[] readBases = BisulfiteHicUtils.getClippedReadsBase(r);
					//int truncateStart = r.getReadPositionAtReferencePosition(r.getAlignmentStart());
					//int truncateEnd = r.getReadPositionAtReferencePosition(r.getAlignmentEnd());
					
					//System.err.println(truncateStart + "\t" + truncateEnd );
					//System.err.println(new String(refBaseEnd1));
					//System.err.println(new String(readBases));
					
					//int pos = CcInferenceUtils.bisulfiteIncompleteReads(readBases, refBaseEnd1, r.getReadNegativeStrandFlag(), r.getReadPairedFlag() & r.getSecondOfPairFlag(), false, true, false, 0.1, 0.4, 0, "CH", 1);
					//System.err.println(r.getReadName() + "\t" + countedReadsEnd1.containsKey(r.getReadName()) + "\t" + r + "\t" + pos);
					//System.err.println(new String(refBaseEnd1));
					//System.err.println(new String(readBases));
				//}
				if(failFlagFilter(r)){
					continue;
				}else{
					//int bisulfiteStartPos = CcInferenceUtils.bisulfiteIncompleteReads(r);
					//if(bisulfiteStartPos<0){
					//	continue;
					//}
					//if(r.getReadName().equalsIgnoreCase("K00168:28:H3NKLBBXX:7:1214:10064:40385")){
					//	System.err.println(r.getReadName() + "\t" + countedReadsEnd2.containsKey(r.getReadName()) + "\t" + r);
					//}
					String name = r.getReadName();
					if(removeReadnameSuffix){
						name = name.substring(0, name.length() - 2);
					}
					if(countedReadsEnd2.containsKey(name)){
						countedReadsEnd2.remove(name); //if two ends are within the same half end of hiccup links
					}else{
						countedReadsEnd2.put(name, r);
					}
				}
			}
			wgsIt.close();
			if(reopenFlag){
				//wgsReader.close();
				//wgsReader = SamReaderFactory.makeDefault().setOption(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES, true).validationStringency(ValidationStringency.SILENT).open(new File(bamFile));
				log.info("reopen bam file to avoid too many connections...");
				methyCorWriter.flush();
				reopenFlag=false;
			}
			//System.err.println(countedReadsEnd1.size() + "\t" + countedReadsEnd2.size());

			int totalPairs = 0;
			for(String readName : countedReadsEnd1.keySet()) {
				if (countedReadsEnd2.containsKey(readName)) {
					totalPairs++;
				}
			}
			if(totalPairs<minReads){
				return; //too few reads, not wirte the pairs...
			}
			int readNum = 1;
			for(String readName : countedReadsEnd1.keySet()){
				if(countedReadsEnd2.containsKey(readName)){
					
					SAMRecord r1 = countedReadsEnd1.get(readName);
					SAMRecord r2 = countedReadsEnd2.get(readName);
					
					Triple<Integer,Integer, Integer> ctSummary1 = BisulfiteHicUtils.readsMethySummaryGeneral(r1, methyPatternSearch, methyPatternCytPos, bsConvPatternSearch, bsConvPatternCytPos, minMapQ, minBaseQ, turnOffBisulfiteFilter, useBadMate);
					Triple<Integer,Integer, Integer> ctSummary2 = BisulfiteHicUtils.readsMethySummaryGeneral(r2, methyPatternSearch, methyPatternCytPos, bsConvPatternSearch, bsConvPatternCytPos, minMapQ, minBaseQ, turnOffBisulfiteFilter, useBadMate);
					int numC1 = ctSummary1.getLeft();
					int numT1 = ctSummary1.getMiddle();
					int numCg1 = ctSummary1.getRight();
					int numC2 = ctSummary2.getLeft();
					int numT2 = ctSummary2.getMiddle();
					int numCg2 = ctSummary2.getRight();
					//System.err.println(chr + "\t" + start1 + "\t" + end1 + "\t" + chr2 + "\t" + start2 + "\t" + end2 + "\t" + (numC1 + numT1) + "\t" + (numC2 + numT2) + "\t" + numCg1 + "\t" + numCg2);
					if(outputCount) {
						
						methyCorWriter.write(chr + "\t" + start1 + "\t" + end1 + "\t" + chr2 + "\t" + start2 + "\t" + end2 + "\t" +
								regionName + "\t" + readNum + "\t" + ob+ "\t" + exp+ "\t" + fdr + "\t" + ((numC1 + numT1) == 0 ? Double.NaN : (double)numC1/(double)(numC1 + numT1)) + "\t" + (numC1 + numT1)
								 + "\t" + ((numC2 + numT2) == 0 ? Double.NaN : (double)numC2/(double)(numC2 + numT2)) + "\t" + (numC2 + numT2) + "\t" + totalReadPairNum + "\n");
								readNum++;
					}else {
						if((numC1 + numT1) == 0 || numCg1 < minCgNum || (numC2 + numT2) == 0 || numCg2 < minCgNum)
							continue;
						
						
						//chr1,read.start1,read.end1,chr2,read.start2,read.end2,region_name, id, hiccups_ob, hiccups_exp, hiccups_fdr, methy_end1, cpg_num_end1, methy_end2, cpg_num_end2
						
						methyCorWriter.write(chr + "\t" + start1 + "\t" + end1 + "\t" + chr2 + "\t" + start2 + "\t" + end2 + "\t" +
									regionName + "\t" + readNum + "\t" + ob+ "\t" + exp+ "\t" + fdr + "\t" + (double)numC1/(double)(numC1 + numT1) + "\t" + (numC1 + numT1)
									 + "\t" + (double)numC2/(double)(numC2 + numT2) + "\t" + (numC2 + numT2) + "\n");
									readNum++;
					}
					
					
					
				}
				
				
			}
			//System.err.println(readNum);
		}
		
		private void methyCorrelationRandom(String chr, int start1, int end1, String chr2, int start2, int end2, SamReader wgsReader, double ob, double exp, double fdr, String bamFile) throws Exception{
			
			String regionName = chr + ":" + start1 + ":" + end1 + ":" + chr2 + ":" + start2 + ":" + end2;
			
			SAMRecordIterator wgsIt = wgsReader.queryOverlapping(chr,start1,end1);
			HashMap<String, SAMRecord> countedReadsEnd1 = new HashMap<String, SAMRecord>();	
			boolean reopenFlag=false;
			while(wgsIt.hasNext()){
				SAMRecord r = wgsIt.next();
				queryNum++;
				if(queryNum % 1000000 == 0){
					reopenFlag=true;
				}
			
				if(failFlagFilter(r)){
					continue;
				}else{
					String name = r.getReadName();
					if(removeReadnameSuffix){
						name = name.substring(0, name.length() - 2);
					}
					if(countedReadsEnd1.containsKey(name)){
						if(!randomNotLinked){//when randomNotLinked, just keep the first one
							countedReadsEnd1.remove(name); //if two ends are within the same half end of hiccup links
						}

					}else{
						countedReadsEnd1.put(name, r);
					}
				}
				
			}
			wgsIt.close();
			
			HashMap<String, SAMRecord> countedReadsEnd2 = new HashMap<String, SAMRecord>();
			wgsIt = wgsReader.queryOverlapping(chr2,start2,end2);
			while(wgsIt.hasNext()){
				SAMRecord r = wgsIt.next();
				queryNum++;
				if(queryNum % 1000000 == 0){
					reopenFlag=true;
				}
				
				if(failFlagFilter(r)){
					continue;
				}else{
					String name = r.getReadName();
					if(removeReadnameSuffix){
						name = name.substring(0, name.length() - 2);
					}
					if(countedReadsEnd2.containsKey(name)){
						if(!randomNotLinked){ //when randomNotLinked, just keep the first one
							countedReadsEnd2.remove(name); //if two ends are within the same half end of hiccup links
						}
					}else{
						countedReadsEnd2.put(name, r);
					}
				}
			}
			wgsIt.close();
			
			//System.err.println(countedReadsEnd1.size() + "\t" + countedReadsEnd2.size());
			int readNum = 0;
			if(randomNotLinked){
				ArrayList<String> readNamesEnd1 = new ArrayList<String>(countedReadsEnd1.keySet());
				ArrayList<String> readNamesEnd2 = new ArrayList<String>(countedReadsEnd2.keySet());
				int minLen = readNamesEnd1.size() < readNamesEnd2.size() ? readNamesEnd1.size() : readNamesEnd2.size();
				if(minLen > minReads) {

						Collections.shuffle(readNamesEnd1);
						for (int j = 0; j < minLen; j++) {
							String readName1 = readNamesEnd1.get(j);
							String readName2 = readNamesEnd2.get(j);
							SAMRecord r1 = countedReadsEnd1.get(readName1);

							SAMRecord r2 = countedReadsEnd2.get(readName2);


							Triple<Integer, Integer, Integer> ctSummary1 = BisulfiteHicUtils.readsMethySummaryGeneral(r1, methyPatternSearch, methyPatternCytPos, bsConvPatternSearch, bsConvPatternCytPos, minMapQ, minBaseQ, turnOffBisulfiteFilter, useBadMate);
							Triple<Integer, Integer, Integer> ctSummary2 = BisulfiteHicUtils.readsMethySummaryGeneral(r2, methyPatternSearch, methyPatternCytPos, bsConvPatternSearch, bsConvPatternCytPos, minMapQ, minBaseQ, turnOffBisulfiteFilter, useBadMate);
							int numC1 = ctSummary1.getLeft();
							int numT1 = ctSummary1.getMiddle();
							int numCg1 = ctSummary1.getRight();
							int numC2 = ctSummary2.getLeft();
							int numT2 = ctSummary2.getMiddle();
							int numCg2 = ctSummary2.getRight();
							if ((numC1 + numT1) == 0 || numCg1 < minCgNum || (numC2 + numT2) == 0 || numCg2 < minCgNum)
								continue;


							//chr1,read.start1,read.end1,chr2,read.start2,read.end2,region_name, id, hiccups_ob, hiccups_exp, hiccups_fdr, methy_end1, cpg_num_end1, methy_end2, cpg_num_end2

							methyCorWriter.write(chr + "\t" + start1 + "\t" + end1 + "\t" + chr2 + "\t" + start2 + "\t" + end2 + "\t" +
									regionName + "\t" + readNum + "\t" + ob + "\t" + exp + "\t" + fdr + "\t" + (double) numC1 / (double) (numC1 + numT1) + "\t" + (numC1 + numT1)
									+ "\t" + (double) numC2 / (double) (numC2 + numT2) + "\t" + (numC2 + numT2) + "\n");


							readNum++;



						}

				}
			}else{
				ArrayList<String> sharedReadNames = new ArrayList<String>();
				for(String readName : countedReadsEnd1.keySet()) {
					if (countedReadsEnd2.containsKey(readName)) {
						sharedReadNames.add(readName);
					}
				}

				String[] readNames = sharedReadNames.toArray(new String[sharedReadNames.size()]);
				//if(readNames1.length>10 && readNames2.length>10){
				if(readNames.length > minReads) {
					for (int i = 0; i < randomSampleReads; i++) {
						Collections.shuffle(sharedReadNames);
						for (int j = 0; j < readNames.length; j++) {
							String readName1 = readNames[j];
							String readName2 = sharedReadNames.get(j);
							SAMRecord r1 = countedReadsEnd1.get(readName1);

							SAMRecord r2 = countedReadsEnd2.get(readName2);


							Triple<Integer, Integer, Integer> ctSummary1 = BisulfiteHicUtils.readsMethySummaryGeneral(r1, methyPatternSearch, methyPatternCytPos, bsConvPatternSearch, bsConvPatternCytPos, minMapQ, minBaseQ, turnOffBisulfiteFilter, useBadMate);
							Triple<Integer, Integer, Integer> ctSummary2 = BisulfiteHicUtils.readsMethySummaryGeneral(r2, methyPatternSearch, methyPatternCytPos, bsConvPatternSearch, bsConvPatternCytPos, minMapQ, minBaseQ, turnOffBisulfiteFilter, useBadMate);
							int numC1 = ctSummary1.getLeft();
							int numT1 = ctSummary1.getMiddle();
							int numCg1 = ctSummary1.getRight();
							int numC2 = ctSummary2.getLeft();
							int numT2 = ctSummary2.getMiddle();
							int numCg2 = ctSummary2.getRight();
							if ((numC1 + numT1) == 0 || numCg1 < minCgNum || (numC2 + numT2) == 0 || numCg2 < minCgNum)
								continue;


							//chr1,read.start1,read.end1,chr2,read.start2,read.end2,region_name, id, hiccups_ob, hiccups_exp, hiccups_fdr, methy_end1, cpg_num_end1, methy_end2, cpg_num_end2

							methyCorWriter.write(chr + "\t" + start1 + "\t" + end1 + "\t" + chr2 + "\t" + start2 + "\t" + end2 + "\t" +
									regionName + "\t" + readNum + "\t" + ob + "\t" + exp + "\t" + fdr + "\t" + (double) numC1 / (double) (numC1 + numT1) + "\t" + (numC1 + numT1)
									+ "\t" + (double) numC2 / (double) (numC2 + numT2) + "\t" + (numC2 + numT2) + "\n");


							readNum++;



						}
					}
				}
			}

			

			//System.err.println(readNum);
		}
		
		private void methyCorrelationShort(String chr, int start1, int end1, int start2, int end2, SamReader wgsReader, double ob, double exp, double fdr, String bamFile) throws Exception{
			
			String regionName = chr + ":" + start1 + ":" + end1 + ":" + start2 + ":" + end2;
			
			SAMRecordIterator wgsIt = wgsReader.queryOverlapping(chr,start1,end1);
			HashMap<String, SAMRecord> countedReadsEnd1 = new HashMap<String, SAMRecord>();	
			HashMap<String, SAMRecord> countedReadsEnd2 = new HashMap<String, SAMRecord>();
			

			while(wgsIt.hasNext()){
				SAMRecord r = wgsIt.next();
				
				if(failFlagFilter(r) || Math.abs(r.getInferredInsertSize())>onlyShort || Math.abs(r.getInferredInsertSize())<shortLowerBound){
					continue;
				}else{
					String name = r.getReadName();
					if(removeReadnameSuffix){
						name = name.substring(0, name.length() - 2);
					}
					if(countedReadsEnd1.containsKey(name)){
						countedReadsEnd2.put(name, r); //if two ends are within the same half end of hiccup links
					}else{
						countedReadsEnd1.put(name, r);
					}
				}
				
			}
			wgsIt.close();

			wgsIt = wgsReader.queryOverlapping(chr,start2,end2);
			while(wgsIt.hasNext()){
				SAMRecord r = wgsIt.next();
				if(failFlagFilter(r) || Math.abs(r.getInferredInsertSize())>onlyShort || Math.abs(r.getInferredInsertSize())<shortLowerBound){
					continue;
				}else{
					String name = r.getReadName();
					if(removeReadnameSuffix){
						name = name.substring(0, name.length() - 2);
					}
					if(countedReadsEnd1.containsKey(name)){
						countedReadsEnd2.put(name, r); //if two ends are within the same half end of hiccup links
					}else{
						countedReadsEnd1.put(name, r);
					}
				}
			}
			wgsIt.close();
			int totalPairs = 0;
			for(String readName : countedReadsEnd1.keySet()) {
				if (countedReadsEnd2.containsKey(readName)) {
					totalPairs++;
				}
			}
			if(totalPairs<minReads){
				return; //too few reads, not wirte the pairs...
			}
			//System.err.println(countedReadsEnd1.size() + "\t" + countedReadsEnd2.size());
			int readNum = 0;
			for(String readName : countedReadsEnd1.keySet()){
				if(countedReadsEnd2.containsKey(readName)){
					
					SAMRecord r1 = countedReadsEnd1.get(readName);
					SAMRecord r2 = countedReadsEnd2.get(readName);
					
					
					Triple<Integer,Integer, Integer> ctSummary1 = BisulfiteHicUtils.readsMethySummaryGeneral(r1, methyPatternSearch, methyPatternCytPos, bsConvPatternSearch, bsConvPatternCytPos, minMapQ, minBaseQ, turnOffBisulfiteFilter, useBadMate);
					Triple<Integer,Integer, Integer> ctSummary2 = BisulfiteHicUtils.readsMethySummaryGeneral(r2, methyPatternSearch, methyPatternCytPos, bsConvPatternSearch, bsConvPatternCytPos, minMapQ, minBaseQ, turnOffBisulfiteFilter, useBadMate);
					int numC1 = ctSummary1.getLeft();
					int numT1 = ctSummary1.getMiddle();
					int numCg1 = ctSummary1.getRight();
					int numC2 = ctSummary2.getLeft();
					int numT2 = ctSummary2.getMiddle();
					int numCg2 = ctSummary2.getRight();
					if((numC1 + numT1) == 0 || numCg1 < minCgNum || (numC2 + numT2) == 0 || numCg2 < minCgNum)
						continue;
					
					
					//chr1,read.start1,read.end1,chr2,read.start2,read.end2,region_name, id, hiccups_ob, hiccups_exp, hiccups_fdr, methy_end1, cpg_num_end1, methy_end2, cpg_num_end2
					
					methyCorWriter.write(chr + "\t" + start1 + "\t" + end1 + "\t" + chr + "\t" + start2 + "\t" + end2 + "\t" + 
								regionName + "\t" + readNum + "\t" + ob+ "\t" + exp+ "\t" + fdr + "\t" + (double)numC1/(double)(numC1 + numT1) + "\t" + (numC1 + numT1)
								 + "\t" + (double)numC2/(double)(numC2 + numT2) + "\t" + (numC2 + numT2) + "\n");
					readNum++;
					
					if(randomSampleReads>0 && readNum > randomSampleReads){
						break;
					}
					
				}
				
				
			}
			//System.err.println(readNum);
		}
		
		private boolean failFlagFilter(SAMRecord r){
			return r.getReadUnmappedFlag() || r.getNotPrimaryAlignmentFlag() || r.getMappingQuality() < minMapQ
					|| r.getReadFailsVendorQualityCheckFlag() || r.getDuplicateReadFlag() || (useUnpaired ? false : !r.getReadPairedFlag()) || (useBadMate ? false : !r.getProperPairFlag());
		}

		
		
		
		private void initiate(String outputPrefix) throws Exception{
			startTime = System.currentTimeMillis();

			
			if(vcf != null){
				asmWriter = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputPrefix + ".asm.txt.gz")), "UTF-8");
				File indexFile = Tribble.tabixIndexFile(new File(vcf));
				if(!indexFile.exists() || !indexFile.canRead()){
					throw new UnsortedFileException(vcf + " file's index " + indexFile.getName() + " does not exist, please use tabix to index it ...");
				}

				vcfReader = new VCFFileReader(new File(vcf), true);
				SAMSequenceDictionary seqDict = VCFFileReader.getSequenceDictionary(new File(vcf));
				glp = new GenomeLocParser(seqDict);
			}else{
				methyCorWriter = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputPrefix + ".methyCor.txt.gz")), "UTF-8");
			}
			
			generator = new MersenneTwister();
		}

		private void finish() throws Exception{
			
			
			
			if(vcf != null){
				vcfReader.close();
				asmWriter.close();
			}else{
				methyCorWriter.close();
			}
			
			
			long endTime   = System.currentTimeMillis();
			double totalTime = endTime - startTime;
			totalTime /= 1000;
			double totalTimeMins = totalTime/60;
			double totalTimeHours = totalTime/3600;
			
			log.info("MethyCorAcrossHiccups's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
		}
		

}
