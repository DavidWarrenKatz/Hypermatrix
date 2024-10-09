/**
 * MethyCorAcrossHiccups.java
 * Sep 12, 2016
 * 1:03:54 PM
 * yaping    lyping1986@gmail.com
 */
package main.java.edu.mit.compbio.bisulfitehic.utils;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;

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
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import edu.unc.genomics.io.BigWigFileReader;

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

public class MethyCorAcrossHiccupsSC {


	@Option(name="-bams",usage="input bams for each single cell. Default: null")
	public ArrayList<String> bams = null;

	@Option(name="-bws",usage="input big wig files for each single cell for the methylation level. Default: null")
	public ArrayList<String> bws = null;

	@Option(name="-bwCounts",usage="input big wig files for each single cell for the total C+T reads covered. Default: null")
	public ArrayList<String> bwCounts = null;
	
	@Option(name="-minCellNum",usage="minimum number of cell have methylation value there. Default: 2")
	public int minCellNum = 2;
	
	
	@Option(name="-useBadMate",usage="use not proper paired reads. Default: false")
	public boolean useBadMate = false;
	
	@Option(name="-minMapQ",usage="minimum mapping quality score required to check. Default: 20")
	public int minMapQ = 20;

	@Option(name="-minDistEnds",usage="minimum distance between two end to be considered as HiC links. Default: 20000")
	public int minDistEnds = 20000;

	@Option(name="-randomHiccupsInput",usage="when the input file is random hiccups links for the control, no need to restrict each end of fragment at the same hiccups links.., default: not enabled")
	public boolean randomHiccupsInput = false;

	@Option(name="-onlyShort",usage="only use the paired end within very short region. when enable it with positive number, it will automate disable -minDistEnds option, default: -1")
	public int onlyShort = -1;
	
	@Option(name="-skipProcessFirstRow",usage="don't process first row, default: not enabled")
	public boolean skipProcessFirstRow = false;

	@Option(name="-useGeneralBedPe",usage="use general bedpe file rather than HiCUUPS link file, default: not enabled")
	public boolean useGeneralBedPe = false;

	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "MethyCorAcrossHiccupsSC [opts] outputPrefix hiccups.bedpe.bed[.gz] cpg_list.bed[.gz]";
	
	private static Logger log = Logger.getLogger(MethyCorAcrossHiccupsSC.class);

	private OutputStreamWriter methyCorWriter = null; 
	private MersenneTwister generator = null;
	
	private static long startTime = -1;
	private static long lineNum=0;
	private static long queryNum=0;
	private static long totalReadPairNum=0;
	
	private List<BigWigFileReader> methyWigFileReaders = null;
	private List<BigWigFileReader> covWigFileReaders = null;
	
	/**
	 * @param args
	 */
	public static void main(String[] args)
		throws Exception {
		MethyCorAcrossHiccupsSC mcah = new MethyCorAcrossHiccupsSC();
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
						if(bams==null || bws==null || bwCounts== null 
								|| bams.isEmpty() || bws.isEmpty() || bwCounts.isEmpty() 
								|| bams.size() != bws.size() || bams.size() != bwCounts.size()) 
							throw new CmdLineException(parser, USAGE, new Throwable());

						//read input bed file, for each row,
						//String intervalFile = arguments.get(0);
						String outputPrefix = arguments.get(0);
						String hiccupsFile = arguments.get(1);
						String cpgListFile = arguments.get(2);
						
						initiate(outputPrefix);
						
						log.info("Processing cpg list file ... ");
						HashMap<String,IntervalTree<String>> cpgList = new HashMap<String,IntervalTree<String>>();

							GZIPInputStream gzipInputStream = null;
							BufferedReader br1;
							if(cpgListFile.endsWith(".gz")){
								gzipInputStream = new GZIPInputStream(new FileInputStream(cpgListFile));
								br1 = new BufferedReader(new InputStreamReader(gzipInputStream));
								
							}else{
								br1 = new BufferedReader(new FileReader(cpgListFile));
							}
								
								String line1;
								
								while( (line1 = br1.readLine()) != null){
									if(line1.startsWith("#"))
										continue;
									String[] splitLines = line1.split("\t");
									if(splitLines.length<3){
										continue;
									}
									String chr = splitLines[0];
									int start = Integer.parseInt(splitLines[1]);
									int end = Integer.parseInt(splitLines[2]);
									
									
									IntervalTree<String> tree;
									
									if(cpgList.containsKey(chr)){
										tree = cpgList.get(chr);
									}else{
										tree = new IntervalTree<String>();
									}
									String strand = ".";
									if(splitLines.length >= 6){
										if(splitLines[5].equalsIgnoreCase("-")){
											strand = "-";
										}else if(splitLines[5].equalsIgnoreCase("+")){
											strand = "+";
										}
										//strand = splitLines[5].equalsIgnoreCase("-") ? "-" : "+";
									}
									tree.put(start, end, strand);
									cpgList.put(chr, tree);
								}
								if(cpgListFile.endsWith(".gz")){
									gzipInputStream.close();
								}
						br1.close();
						//SamReader wgsReader = SamReaderFactory.makeDefault().setOption(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES, true).validationStringency(ValidationStringency.SILENT).open(new File(bamFile));
						//SamReader wgsReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamFile));
						//SamReaderFactory factory = SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES).validationStringency(ValidationStringency.SILENT);
						SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);

						log.info("Parsing input hiccups links BEDPE format bed file ...");
						 gzipInputStream = null;
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
								continue;
							}else{
								String[] splitin = line.split("\t");
								String chr1 = splitin[0];
								int start1 = Integer.parseInt(splitin[1]);
								int end1 = Integer.parseInt(splitin[2]);
								String chr2 = splitin[3];
								int start2 = Integer.parseInt(splitin[4]);
								int end2 = Integer.parseInt(splitin[5]);
								if((!chr1.equalsIgnoreCase(chr2)) || (onlyShort<0 && Math.abs(start1-start2)<minDistEnds) || (onlyShort<0 && Math.abs(end1-end2)<minDistEnds)){
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
								
								
									if(randomHiccupsInput){
										//log.info("randomHiccupsInput mode: ");
										methyCorrelationRandom(chr1, start1, end1, start2, end2, ob, exp, fdr, cpgList, factory);
									}else{
										if(onlyShort>0){
											//log.info("only use paired end reads within short distance: ");
											methyCorrelationShort(chr1, start1, end1, start2, end2, ob, exp, fdr, cpgList, factory);
										}else{
											methyCorrelation(chr1, start1, end1, start2, end2, ob, exp, fdr, cpgList, factory);
										}
										
									}
								
								
								
								
							}
								
							
							lineNum++;
							if(lineNum % 1000 == 0){
								log.info("Processing line: " + lineNum);
								methyCorWriter.flush();
							}
							
						}
						if(hiccupsFile.endsWith(".gz")){
							gzipInputStream.close();
						}
						br.close();
						
						finish();
		}
		
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

		
		private void methyCorrelation(String chr, int start1, int end1, int start2, int end2,  double ob, double exp, double fdr, HashMap<String,IntervalTree<String>> cpgList
				,SamReaderFactory factory) throws Exception{
			
			
			HashMap<String,ArrayList<Integer>> cpgInterEachCell = new HashMap<String,ArrayList<Integer>>();
			String regionName = chr + ":" + start1 + ":" + end1 + ":" + start2 + ":" + end2;
			IntervalTree<String> cgs = cpgList.get(chr);
			//get a template for each cell.
			ArrayList<Integer> cpgEachCellTemp = new ArrayList<Integer>(bams.size());
			for(int i = 0; i < bams.size(); i++) {
				cpgEachCellTemp.add(0);
			}
			//go over each bam, determin the cpg pairs candidate
			for(int i = 0; i < bams.size(); i++) {
				String bamFile = bams.get(i);
				
				BigWigFileReader covReader = covWigFileReaders.get(i);
				
				SamReader wgsReader = factory.open(new File(bamFile));
				
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
				wgsReader.close();
				for(String readName : countedReadsEnd1.keySet()){
					if(countedReadsEnd2.containsKey(readName)){
						
						SAMRecord r1 = countedReadsEnd1.get(readName);
						SAMRecord r2 = countedReadsEnd2.get(readName);
						int startR1 = r1.getAlignmentStart();
						int endR1 = r1.getAlignmentEnd();
						int startR2 = r2.getAlignmentStart();
						int endR2 = r2.getAlignmentEnd();
						
						Iterator<Node<String>> it1 = null;
						Iterator<Node<String>> it2 = null;
							
						if(cgs.minOverlapper(startR1, endR1)==null || cgs.minOverlapper(startR2, endR2)==null){
								continue;
						}else{
								it1 = cgs.overlappers(startR1, endR1);
								it2 = cgs.overlappers(startR2, endR2);
								ArrayList<Integer> cg1List = new ArrayList<Integer>();
								ArrayList<Integer> cg2List = new ArrayList<Integer>();
								while(it1.hasNext()) {
									int cor1 = it1.next().getEnd();
									double v = covReader.queryStats(chr, cor1, cor1).getSum();
									if(Double.isNaN(v) || Double.compare(v, 0.0) == 0) {
										continue;
									}else {
										cg1List.add(cor1);
									}
									
								}
								while(it2.hasNext()) {
									int cor2 = it2.next().getEnd();
									double v = covReader.queryStats(chr, cor2, cor2).getSum();
									if(Double.isNaN(v) || Double.compare(v, 0.0) == 0) {
										continue;
									}else {
										cg2List.add(cor2);
									}
								}
								for(int cg1 : cg1List) {
									for(int cg2 : cg2List) {
										String cpgPairName = cg1 + ":" + cg2;
										if(cpgInterEachCell.containsKey(cpgPairName)) {
											ArrayList<Integer> cpgEachCell = cpgInterEachCell.get(cpgPairName);
											cpgEachCell.set(i, 1);
											cpgInterEachCell.put(cpgPairName, cpgEachCell);
										}else {
											ArrayList<Integer> cpgEachCell = (ArrayList<Integer>) cpgEachCellTemp.clone();
											cpgEachCell.set(i, 1);
											cpgInterEachCell.put(cpgPairName, cpgEachCell);
										}
									}
								}
						}
					}
				}
			}
			
			//after determine all of the possible contacted cpg pairs in this loop, go over each of them, get methy value and output them:
			for(String cpgPairs : cpgInterEachCell.keySet()) {
				String[] cpgPair = cpgPairs.split(":");
				int cg1 = Integer.parseInt(cpgPair[0]);
				int cg2 = Integer.parseInt(cpgPair[1]);
				ArrayList<Integer> cpgEachCell = cpgInterEachCell.get(cpgPairs);
				int interactSum = 0;
				for(int cpg : cpgEachCell) {
					interactSum+=cpg;
				}
				if(interactSum < minCellNum) {
					continue;
				}
				for(int i = 0; i < cpgEachCell.size(); i++) {
					if(cpgEachCell.get(i) > 0) {
						BigWigFileReader methyReader = methyWigFileReaders.get(i);
						BigWigFileReader covReader = covWigFileReaders.get(i);
						double count1 = covReader.queryStats(chr, cg1, cg1).getSum();
						double methy1 = methyReader.queryStats(chr, cg1, cg1).getMean();
						double count2 = covReader.queryStats(chr, cg2, cg2).getSum();
						double methy2 = methyReader.queryStats(chr, cg2, cg2).getMean();
						methyCorWriter.write(chr + "\t" + (cg1-1) + "\t" + cg1 + "\t" + chr + "\t" + (cg2-1) + "\t" + cg2 + "\t" + 
								regionName + "\t" + i + "\t" + ob+ "\t" + exp+ "\t" + fdr + "\t" + methy1 + "\t" + count1
								 + "\t" + methy2 + "\t" + count2 + "\n");
					}
				}
				
			}
			
		}
		
		private void methyCorrelationRandom(String chr, int start1, int end1, int start2, int end2,  double ob, double exp, double fdr, HashMap<String,IntervalTree<String>> cpgList
				, SamReaderFactory factory) throws Exception{
			
			HashMap<String,ArrayList<Integer>> cpgInterEachCell = new HashMap<String,ArrayList<Integer>>();
			String regionName = chr + ":" + start1 + ":" + end1 + ":" + start2 + ":" + end2;
			IntervalTree<String> cgs = cpgList.get(chr);
			//get a template for each cell.
			ArrayList<Integer> cpgEachCellTemp = new ArrayList<Integer>(bams.size());
			for(int i = 0; i < bams.size(); i++) {
				cpgEachCellTemp.add(0);
			}
			//go over each bam, determin the cpg pairs candidate
			//SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
			for(int i = 0; i < bams.size(); i++) {
				String bamFile = bams.get(i);
				
				BigWigFileReader covReader = covWigFileReaders.get(i);
				
				SamReader wgsReader = factory.open(new File(bamFile));
				
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
				wgsReader.close();
			//System.err.println(countedReadsEnd1.size() + "\t" + countedReadsEnd2.size());
			int readNum = 0;
			String[] readNames1 = countedReadsEnd1.keySet().toArray(new String[countedReadsEnd1.keySet().size()]);
			String[] readNames2 = countedReadsEnd2.keySet().toArray(new String[countedReadsEnd2.keySet().size()]);
			if(readNames1.length>10 && readNames2.length>10){
				for(String readName : readNames1){
					SAMRecord r1 = countedReadsEnd1.get(readName);
					int randomItem = generator.nextInt(readNames2.length);
					while(readNames2[randomItem].equalsIgnoreCase(readName)){
						randomItem = generator.nextInt(readNames2.length);
					}
					SAMRecord r2=countedReadsEnd2.get(readNames2[randomItem]);
					
					
					int startR1 = r1.getAlignmentStart();
					int endR1 = r1.getAlignmentEnd();
					int startR2 = r2.getAlignmentStart();
					int endR2 = r2.getAlignmentEnd();
					
					Iterator<Node<String>> it1 = null;
					Iterator<Node<String>> it2 = null;
						
					if(cgs.minOverlapper(startR1, endR1)==null || cgs.minOverlapper(startR2, endR2)==null){
							continue;
					}else{
							it1 = cgs.overlappers(startR1, endR1);
							it2 = cgs.overlappers(startR2, endR2);
							ArrayList<Integer> cg1List = new ArrayList<Integer>();
							ArrayList<Integer> cg2List = new ArrayList<Integer>();
							while(it1.hasNext()) {
								int cor1 = it1.next().getEnd();
								double v = covReader.queryStats(chr, cor1, cor1).getSum();
								if(Double.isNaN(v) || Double.compare(v, 0.0) == 0) {
									continue;
								}else {
									cg1List.add(cor1);
								}
								
							}
							while(it2.hasNext()) {
								int cor2 = it2.next().getEnd();
								double v = covReader.queryStats(chr, cor2, cor2).getSum();
								if(Double.isNaN(v) || Double.compare(v, 0.0) == 0) {
									continue;
								}else {
									cg2List.add(cor2);
								}
							}
							for(int cg1 : cg1List) {
								for(int cg2 : cg2List) {
									String cpgPairName = cg1 + ":" + cg2;
									if(cpgInterEachCell.containsKey(cpgPairName)) {
										ArrayList<Integer> cpgEachCell = cpgInterEachCell.get(cpgPairName);
										cpgEachCell.set(i, 1);
										cpgInterEachCell.put(cpgPairName, cpgEachCell);
									}else {
										ArrayList<Integer> cpgEachCell = (ArrayList<Integer>) cpgEachCellTemp.clone();
										cpgEachCell.set(i, 1);
										cpgInterEachCell.put(cpgPairName, cpgEachCell);
									}
								}
							}
					}
				}
			}
		}
		
		//after determine all of the possible contacted cpg pairs in this loop, go over each of them, get methy value and output them:
		for(String cpgPairs : cpgInterEachCell.keySet()) {
			String[] cpgPair = cpgPairs.split(":");
			int cg1 = Integer.parseInt(cpgPair[0]);
			int cg2 = Integer.parseInt(cpgPair[1]);
			ArrayList<Integer> cpgEachCell = cpgInterEachCell.get(cpgPairs);
			int interactSum = 0;
			for(int cpg : cpgEachCell) {
				interactSum+=cpg;
			}
			if(interactSum < minCellNum) {
				continue;
			}
			for(int i = 0; i < cpgEachCell.size(); i++) {
				if(cpgEachCell.get(i) > 0) {
					BigWigFileReader methyReader = methyWigFileReaders.get(i);
					BigWigFileReader covReader = covWigFileReaders.get(i);
					double count1 = covReader.queryStats(chr, cg1, cg1).getSum();
					double methy1 = methyReader.queryStats(chr, cg1, cg1).getMean();
					double count2 = covReader.queryStats(chr, cg2, cg2).getSum();
					double methy2 = methyReader.queryStats(chr, cg2, cg2).getMean();
					methyCorWriter.write(chr + "\t" + (cg1-1) + "\t" + cg1 + "\t" + chr + "\t" + (cg2-1) + "\t" + cg2 + "\t" + 
							regionName + "\t" + i + "\t" + ob+ "\t" + exp+ "\t" + fdr + "\t" + methy1 + "\t" + count1
							 + "\t" + methy2 + "\t" + count2 + "\n");
				}
			}
			
		}
		
	}
	
		
		private void methyCorrelationShort(String chr, int start1, int end1, int start2, int end2, double ob, double exp, double fdr, HashMap<String,IntervalTree<String>> cpgList
				, SamReaderFactory factory) throws Exception{
			
			
			HashMap<String,ArrayList<Integer>> cpgInterEachCell = new HashMap<String,ArrayList<Integer>>();
			String regionName = chr + ":" + start1 + ":" + end1 + ":" + start2 + ":" + end2;
			IntervalTree<String> cgs = cpgList.get(chr);
			//get a template for each cell.
			ArrayList<Integer> cpgEachCellTemp = new ArrayList<Integer>(bams.size());
			for(int i = 0; i < bams.size(); i++) {
				cpgEachCellTemp.add(0);
			}
			//go over each bam, determin the cpg pairs candidate
			for(int i = 0; i < bams.size(); i++) {
				String bamFile = bams.get(i);
				
				BigWigFileReader covReader = covWigFileReaders.get(i);
				
				SamReader wgsReader = factory.open(new File(bamFile));
				
				SAMRecordIterator wgsIt = wgsReader.queryOverlapping(chr,start1,end1);
				HashMap<String, SAMRecord> countedReadsEnd1 = new HashMap<String, SAMRecord>();	
				while(wgsIt.hasNext()){
					SAMRecord r = wgsIt.next();
					if(failFlagFilter(r) || Math.abs(r.getInferredInsertSize()) > onlyShort){
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
					if(failFlagFilter(r) || Math.abs(r.getInferredInsertSize())  > onlyShort){
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
				wgsReader.close();
				for(String readName : countedReadsEnd1.keySet()){
					if(countedReadsEnd2.containsKey(readName)){
						
						SAMRecord r1 = countedReadsEnd1.get(readName);
						SAMRecord r2 = countedReadsEnd2.get(readName);
						int startR1 = r1.getAlignmentStart();
						int endR1 = r1.getAlignmentEnd();
						int startR2 = r2.getAlignmentStart();
						int endR2 = r2.getAlignmentEnd();
						
						Iterator<Node<String>> it1 = null;
						Iterator<Node<String>> it2 = null;
							
						if(cgs.minOverlapper(startR1, endR1)==null || cgs.minOverlapper(startR2, endR2)==null){
								continue;
						}else{
								it1 = cgs.overlappers(startR1, endR1);
								it2 = cgs.overlappers(startR2, endR2);
								ArrayList<Integer> cg1List = new ArrayList<Integer>();
								ArrayList<Integer> cg2List = new ArrayList<Integer>();
								while(it1.hasNext()) {
									int cor1 = it1.next().getEnd();
									double v = covReader.queryStats(chr, cor1, cor1).getSum();
									if(Double.isNaN(v) || Double.compare(v, 0.0) == 0) {
										continue;
									}else {
										cg1List.add(cor1);
									}
									
								}
								while(it2.hasNext()) {
									int cor2 = it2.next().getEnd();
									double v = covReader.queryStats(chr, cor2, cor2).getSum();
									if(Double.isNaN(v) || Double.compare(v, 0.0) == 0) {
										continue;
									}else {
										cg2List.add(cor2);
									}
								}
								for(int cg1 : cg1List) {
									for(int cg2 : cg2List) {
										String cpgPairName = cg1 + ":" + cg2;
										if(cpgInterEachCell.containsKey(cpgPairName)) {
											ArrayList<Integer> cpgEachCell = cpgInterEachCell.get(cpgPairName);
											cpgEachCell.set(i, 1);
											cpgInterEachCell.put(cpgPairName, cpgEachCell);
										}else {
											ArrayList<Integer> cpgEachCell = (ArrayList<Integer>) cpgEachCellTemp.clone();
											cpgEachCell.set(i, 1);
											cpgInterEachCell.put(cpgPairName, cpgEachCell);
										}
									}
								}
						}
					}
				}
			}
			
			//after determine all of the possible contacted cpg pairs in this loop, go over each of them, get methy value and output them:
			for(String cpgPairs : cpgInterEachCell.keySet()) {
				String[] cpgPair = cpgPairs.split(":");
				int cg1 = Integer.parseInt(cpgPair[0]);
				int cg2 = Integer.parseInt(cpgPair[1]);
				ArrayList<Integer> cpgEachCell = cpgInterEachCell.get(cpgPairs);
				int interactSum = 0;
				for(int cpg : cpgEachCell) {
					interactSum+=cpg;
				}
				if(interactSum < minCellNum) {
					continue;
				}
				for(int i = 0; i < cpgEachCell.size(); i++) {
					if(cpgEachCell.get(i) > 0) {
						BigWigFileReader methyReader = methyWigFileReaders.get(i);
						BigWigFileReader covReader = covWigFileReaders.get(i);
						double count1 = covReader.queryStats(chr, cg1, cg1).getSum();
						double methy1 = methyReader.queryStats(chr, cg1, cg1).getMean();
						double count2 = covReader.queryStats(chr, cg2, cg2).getSum();
						double methy2 = methyReader.queryStats(chr, cg2, cg2).getMean();
						methyCorWriter.write(chr + "\t" + (cg1-1) + "\t" + cg1 + "\t" + chr + "\t" + (cg2-1) + "\t" + cg2 + "\t" + 
								regionName + "\t" + i + "\t" + ob+ "\t" + exp+ "\t" + fdr + "\t" + methy1 + "\t" + count1
								 + "\t" + methy2 + "\t" + count2 + "\n");
					}
				}
				
			}
			
		}
		
		private boolean failFlagFilter(SAMRecord r){
			return r.getReadUnmappedFlag() || r.getNotPrimaryAlignmentFlag() || r.getMappingQuality() < minMapQ
					|| r.getReadFailsVendorQualityCheckFlag() || r.getDuplicateReadFlag() || !r.getReadPairedFlag() || (useBadMate ? false : !r.getProperPairFlag());
		}

		
		
		
		private void initiate(String outputPrefix) throws Exception{
			startTime = System.currentTimeMillis();
				methyCorWriter = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputPrefix + ".methyCor.txt.gz")), "UTF-8");
			
			
			generator = new MersenneTwister();
			
			methyWigFileReaders = new ArrayList<BigWigFileReader>();
			covWigFileReaders = new ArrayList<BigWigFileReader>();	
			for(int i = 0; i < bws.size(); i++){
					BigWigFileReader bbreader = new BigWigFileReader((new File(bws.get(i))).toPath()); 
					methyWigFileReaders.add(bbreader);
					BigWigFileReader bbreader2 = new BigWigFileReader((new File(bwCounts.get(i))).toPath());
					covWigFileReaders.add(bbreader2);
			}
		}
		

		private void finish() throws Exception{
			
			if(methyWigFileReaders != null){
				for(int i = 0; i < methyWigFileReaders.size(); i++){
					methyWigFileReaders.get(i).close();
				}
			}
			if(covWigFileReaders != null){
				for(int i = 0; i < covWigFileReaders.size(); i++){
					covWigFileReaders.get(i).close();
				}
			}
			
				methyCorWriter.close();
			
			
			
			long endTime   = System.currentTimeMillis();
			double totalTime = endTime - startTime;
			totalTime /= 1000;
			double totalTimeMins = totalTime/60;
			double totalTimeHours = totalTime/3600;
			
			log.info("MethyCorAcrossHiccupsSC's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
		}
		

}
