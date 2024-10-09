/**
 * MethyCorAcrossHiccups.java
 * Sep 12, 2016
 * 1:03:54 PM
 * yaping    lyping1986@gmail.com
 */
package main.java.edu.mit.compbio.bisulfitehic.utils;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.exception.UnsortedFileException;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import main.java.edu.mit.compbio.ccinference.utils.CcInferenceUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
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
import scala.tools.cmd.gen.AnyVals;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

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

public class MethyCorGlobal {

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

	
	@Option(name="-minCgNum",usage="For ASM mode, specify the minimum  number of CpG to call ASM. Default: 1")
	public int minCgNum = 1;


	@Option(name="-randomHiccupsInput",usage="when the input file is random hiccups links for the control, no need to restrict each end of fragment at the same hiccups links.., default: not enabled")
	public boolean randomHiccupsInput = false;


	@Option(name="-randomSampleReads",usage="only use part of the reads in each interval. good for -randomHiccupsInput mode or -onlyShort mode, default: -1. not enabled")
	public int randomSampleReads = -1;


	@Option(name="-aveWindow",usage="the window size to average the test statistics. default value do not use window to average. Default: 1000000")
	public int aveWindow = 1000000;


	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "MethyCorGlobal [opts] outputPrefix mm9.chroms.size chr input.bam";
	
	private static Logger log = Logger.getLogger(MethyCorGlobal.class);

	private OutputStreamWriter methyCorWriter = null;
	private MersenneTwister generator = null;

	
	private static long startTime = -1;
	private static long queryNum=0;

	/**
	 * @param args
	 */
	public static void main(String[] args)
		throws Exception {
		MethyCorGlobal mcah = new MethyCorGlobal();
			BasicConfigurator.configure();
			mcah.doMain(args);
		}
		

		public void doMain(String[] args)
				throws Exception {

						CmdLineParser parser = new CmdLineParser(this);
						//parser.setUsageWidth(80);
						try
						{
							if(help || args.length < 4) throw new CmdLineException(parser, USAGE, new Throwable());
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
			String chromSizeFile = arguments.get(1);
			String chr = arguments.get(2);
			String bamFile = arguments.get(3);

			initiate(outputPrefix);
			SamReader wgsReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamFile));

			PearsonsCorrelation pearson = new PearsonsCorrelation();
			BufferedReader br = new BufferedReader(new FileReader(chromSizeFile));
			String line;
			while( (line = br.readLine()) != null){
				String[] splitin = line.split("\t");
				if(splitin[0].equalsIgnoreCase(chr)){
					int size = Integer.parseInt(splitin[1]);
					for(int start1 = 0; start1 < size-aveWindow; start1 += aveWindow){
						int end1 = start1 + aveWindow;
						for(int start2 = 0; start2 < size-aveWindow; start2 += aveWindow){
							int end2 = start2 + aveWindow;
							double cor = Double.NaN;
							if(randomHiccupsInput) {
								cor = methyCorrelationRandom(chr, start1, end1, start2, end2, pearson, wgsReader);

							}else{
								cor = methyCorrelation(chr, start1, end1, start2, end2, pearson, wgsReader);
							}


							methyCorWriter.write(cor + "\t");
						}
						methyCorWriter.write("\n");
						methyCorWriter.flush();
					}
				}
			}

			wgsReader.close();
			br.close();

			finish();
		}


		
		private double methyCorrelation(String chr, int start1, int end1, int start2, int end2, PearsonsCorrelation pearson, SamReader wgsReader) throws Exception{
			
			SAMRecordIterator wgsIt = wgsReader.queryOverlapping(chr,start1,end1);
			HashMap<String, SAMRecord> countedReadsEnd1 = new HashMap<String, SAMRecord>();
			while(wgsIt.hasNext()){
				SAMRecord r = wgsIt.next();
				queryNum++;


				if(failFlagFilter(r)){
					continue;
				}else{

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

			HashMap<String, SAMRecord> countedReadsEnd2 = new HashMap<String, SAMRecord>();
			wgsIt = wgsReader.queryOverlapping(chr,start2,end2);
			while(wgsIt.hasNext()){
				SAMRecord r = wgsIt.next();
				queryNum++;

				if(failFlagFilter(r)){
					continue;
				}else{

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

			//System.err.println(countedReadsEnd1.size() + "\t" + countedReadsEnd2.size());
			ArrayList<Double> regionArray1 = new ArrayList<Double>();
			ArrayList<Double> regionArray2 = new ArrayList<Double>();
			for(String readName : countedReadsEnd1.keySet()){
				if(countedReadsEnd2.containsKey(readName)){
					
					SAMRecord r1 = countedReadsEnd1.get(readName);
					SAMRecord r2 = countedReadsEnd2.get(readName);
					
					Triple<Integer,Integer, Integer> ctSummary1 = BisulfiteHicUtils.readsMethySummary(r1, minMapQ, minBaseQ, turnOffBisulfiteFilter);
					Triple<Integer,Integer, Integer> ctSummary2 = BisulfiteHicUtils.readsMethySummary(r2, minMapQ, minBaseQ, turnOffBisulfiteFilter);
					int numC1 = ctSummary1.getLeft();
					int numT1 = ctSummary1.getMiddle();
					int numCg1 = ctSummary1.getRight();
					int numC2 = ctSummary2.getLeft();
					int numT2 = ctSummary2.getMiddle();
					int numCg2 = ctSummary2.getRight();
					if((numC1 + numT1) == 0 || numCg1 < minCgNum || (numC2 + numT2) == 0 || numCg2 < minCgNum)
						continue;
					regionArray1.add((double)numC1/(double)(numC1+numT1));
					regionArray2.add((double)numC2/(double)(numC2+numT2));
				}
					
					
					
			}
			double[] region1 = ArrayUtils.toPrimitive(regionArray1.toArray(new Double[regionArray1.size()]));
			double[] region2 = ArrayUtils.toPrimitive(regionArray2.toArray(new Double[regionArray2.size()]));
			if(regionArray1.size()>=2 && regionArray2.size()>=2){
				return pearson.correlation(region1,region2);
			}else{
				return Double.NaN;
			}


		}
		
		private double methyCorrelationRandom(String chr, int start1, int end1, int start2, int end2, PearsonsCorrelation pearson, SamReader wgsReader) throws Exception{

			
			SAMRecordIterator wgsIt = wgsReader.queryOverlapping(chr,start1,end1);
			HashMap<String, SAMRecord> countedReadsEnd1 = new HashMap<String, SAMRecord>();
			while(wgsIt.hasNext()){
				SAMRecord r = wgsIt.next();
				queryNum++;

			
				if(failFlagFilter(r)){
					continue;
				}else{
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
			
			HashMap<String, SAMRecord> countedReadsEnd2 = new HashMap<String, SAMRecord>();
			wgsIt = wgsReader.queryOverlapping(chr,start2,end2);
			while(wgsIt.hasNext()){
				SAMRecord r = wgsIt.next();
				queryNum++;

				
				if(failFlagFilter(r)){
					continue;
				}else{
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
			
			//System.err.println(countedReadsEnd1.size() + "\t" + countedReadsEnd2.size());
			int readNum = 0;
			String[] readNames1 = countedReadsEnd1.keySet().toArray(new String[countedReadsEnd1.keySet().size()]);
			String[] readNames2 = countedReadsEnd2.keySet().toArray(new String[countedReadsEnd2.keySet().size()]);
			ArrayList<Double> regionArray1 = new ArrayList<Double>();
			ArrayList<Double> regionArray2 = new ArrayList<Double>();
			if(readNames1.length>10 && readNames2.length>10){
				for(String readName : readNames1){
					SAMRecord r1 = countedReadsEnd1.get(readName);
					int randomItem = generator.nextInt(readNames2.length);
					while(readNames2[randomItem].equalsIgnoreCase(readName)){
						randomItem = generator.nextInt(readNames2.length);
					}
					SAMRecord r2=countedReadsEnd2.get(readNames2[randomItem]);
					
					
					Triple<Integer,Integer, Integer> ctSummary1 = BisulfiteHicUtils.readsMethySummary(r1, minMapQ, minBaseQ, turnOffBisulfiteFilter);
					Triple<Integer,Integer, Integer> ctSummary2 = BisulfiteHicUtils.readsMethySummary(r2, minMapQ, minBaseQ, turnOffBisulfiteFilter);
					int numC1 = ctSummary1.getLeft();
					int numT1 = ctSummary1.getMiddle();
					int numCg1 = ctSummary1.getRight();
					int numC2 = ctSummary2.getLeft();
					int numT2 = ctSummary2.getMiddle();
					int numCg2 = ctSummary2.getRight();
					if((numC1 + numT1) == 0 || numCg1 < minCgNum || (numC2 + numT2) == 0 || numCg2 < minCgNum)
						continue;

					regionArray1.add((double)numC1/(double)(numC1+numT1));
					regionArray2.add((double)numC2/(double)(numC2+numT2));
					
					if(randomSampleReads>0 && readNum > randomSampleReads){
						break;
					}
					
					
				}
			}
			double[] region1 = ArrayUtils.toPrimitive(regionArray1.toArray(new Double[regionArray1.size()]));
			double[] region2 = ArrayUtils.toPrimitive(regionArray2.toArray(new Double[regionArray2.size()]));
			if(regionArray1.size()>=2 && regionArray2.size()>=2){
				return pearson.correlation(region1,region2);
			}else{
				return Double.NaN;
			}


			//System.err.println(readNum);
		}

		
		private boolean failFlagFilter(SAMRecord r){
			return r.getReadUnmappedFlag() || r.getNotPrimaryAlignmentFlag() || r.getMappingQuality() < minMapQ
					|| r.getReadFailsVendorQualityCheckFlag() || r.getDuplicateReadFlag() || (useUnpaired ? false : !r.getReadPairedFlag()) || (useBadMate ? false : !r.getProperPairFlag());
		}

		
		
		
		private void initiate(String outputPrefix) throws Exception{
			startTime = System.currentTimeMillis();

			

				methyCorWriter = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputPrefix + ".methyCor.txt.gz")), "UTF-8");

			
			generator = new MersenneTwister();
		}

		private void finish() throws Exception{
			
			
			

				methyCorWriter.close();

			
			
			long endTime   = System.currentTimeMillis();
			double totalTime = endTime - startTime;
			totalTime /= 1000;
			double totalTimeMins = totalTime/60;
			double totalTimeHours = totalTime/3600;
			
			log.info("MethyCorGlobal's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
		}
		

}
