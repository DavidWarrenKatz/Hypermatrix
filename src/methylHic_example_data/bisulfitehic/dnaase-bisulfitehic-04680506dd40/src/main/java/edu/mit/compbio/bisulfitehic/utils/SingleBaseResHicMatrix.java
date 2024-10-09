/**
 * enhanceHicMatrix.java
 * Nov 8, 2016
 * 9:26:51 AM
 * yaping    lyping1986@gmail.com
 */
package main.java.edu.mit.compbio.bisulfitehic.utils;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
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
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import main.java.edu.mit.compbio.ccinference.utils.CcInferenceUtils;

import org.apache.commons.collections4.IteratorUtils;
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.util.Pair;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

/**
 * use the concordance of the methylation status at each end to increase the resolution of HiC. 
 * output as bed format (gzipped):
 * chr6    80000_120000    120000_160000   10.278
 *chr6    80000_120000    160000_200000   3.648
 *chr6    80000_120000    200000_240000   4.204
 *
 * <readname> <str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2> <mapq1> <mapq2>
 * 
 */
public class SingleBaseResHicMatrix {

	
	@Option(name="-minBaseQ",usage="minimum base quality score required to check. Default: 5")
	public int minBaseQ = 5;
	
	@Option(name="-turnOffBisulfiteFilter",usage="turn off bisulfite conversion reads filter. Default: false")
	public boolean turnOffBisulfiteFilter = false;

	@Option(name="-useBadMate",usage="useBadMate. Default: false")
	public boolean useBadMate = false;

	
	@Option(name="-minMapQ",usage="minimum mapping quality score required to check. Default: 30")
	public int minMapQ = 30;

	@Option(name="-minDistEnds",usage="minimum distance between two end to be output. Default: 1000")
	public int minDistEnds = 1000;

	@Option(name="-maxDistEnds",usage="maximum distance between two end to be output. Default: 200000000")
	public int maxDistEnds = 200000000;

	@Option(name="-bedpe",usage="bedpe file to limit the analysis only in these regions. Otherwise, it will use all of the cpg provided. Default: null")
	public String bedpe = null;
	
	@Option(name="-minReadsInWindow",usage="minimum number of reads required to calculate pearson correlation within the window at both ends Default: 3")
	public int minReadsInWindow   = 3;
	
	//@Option(name="-minCgNum",usage="For ASM mode, specify the minimum  number of CpG to call ASM. Default: 1")
	//public int minCgNum = 1;
	
	@Option(name="-bgBam",usage="provide background bam file to check the permutation p value, rather than use input.bam file. Default: null")
	public String bgBam = null;
	
	@Option(name="-permutations",usage="number of permutations to generate expectated p value within the same segments, default: 100")
	public int permutations = 100;
	
	@Option(name="-skipProcessFirstRow",usage="don't process first row, default: not enabled")
	public boolean skipProcessFirstRow = false;
	
	@Option(name="-regions",usage="only one side is overlapped with -regions is fine,bed format, default: null")
	public String regions = null;

	@Option(name="-h",usage="show option information")
	public boolean help = false;

	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "SingleBaseResHicMatrix [opts] outputPrefix cpg_list.bed[.gz] input.bam";
	
	private static Logger log = Logger.getLogger(SingleBaseResHicMatrix.class);

	private OutputStreamWriter writer = null; 
	private MersenneTwister generator = null;
	
	
	private static long startTime = -1;
	private static long queryNum=0;

	/**
	 * @param args
	 */
	public static void main(String[] args)
		throws Exception {
			SingleBaseResHicMatrix ehm = new SingleBaseResHicMatrix();
				BasicConfigurator.configure();
				ehm.doMain(args);
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
							String outputPrefix = arguments.get(0);
							String cpgFile = arguments.get(1);
							String bamFile = arguments.get(2);

							initiate(outputPrefix);
							SamReader wgsReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamFile));
							SAMFileHeader samFileHeader = wgsReader.getFileHeader();
							SAMSequenceDictionary dict = samFileHeader.getSequenceDictionary();
							SamReader wgsBgBamReader = null;
							if(bgBam != null){
								wgsBgBamReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bgBam));
							}
							
							HashMap<String,IntervalTree<String>> regionsOverlapped = new HashMap<String,IntervalTree<String>>();
							if(regions != null){
								log.info("Parsing input -regions bed file ...");
								
								GZIPInputStream gzipInputStream1 = null;
								BufferedReader br;
								if(regions.endsWith(".gz")){
									gzipInputStream1 = new GZIPInputStream(new FileInputStream(regions));
									br = new BufferedReader(new InputStreamReader(gzipInputStream1));
									
								}else{
									br = new BufferedReader(new FileReader(regions));
								}
									
									String line;
									long i = 0;
									while( (line = br.readLine()) != null){
										if(line.startsWith("#") || (skipProcessFirstRow && i == 0))
											continue;
										String[] splitin = line.split("\t");
										String chr = splitin[0];
										int start = Integer.parseInt(splitin[1]);
										int end = Integer.parseInt(splitin[2]);
										
									
										IntervalTree<String> tree;
										
										if(regionsOverlapped.containsKey(chr)){
											tree = regionsOverlapped.get(chr);
										}else{
											tree = new IntervalTree<String>();
										}
										
										tree.put(start, end, "1");
										
										
										i++;
										
									}
									if(regions.endsWith(".gz")){
										gzipInputStream1.close();
									}
									br.close();
							}
							
							HashMap<String,IntervalTree<IntervalList>> intervalCollections = new HashMap<String,IntervalTree<IntervalList>>();
							if(bedpe != null){
								log.info("Parsing input bedpe bed file ...");
								
								GZIPInputStream gzipInputStream1 = null;
								BufferedReader br;
								if(bedpe.endsWith(".gz")){
									gzipInputStream1 = new GZIPInputStream(new FileInputStream(bedpe));
									br = new BufferedReader(new InputStreamReader(gzipInputStream1));
									
								}else{
									br = new BufferedReader(new FileReader(bedpe));
								}
									
									String line;
									long i = 0;
									while( (line = br.readLine()) != null){
										if(line.startsWith("#") || (skipProcessFirstRow && i == 0))
											continue;
										String[] splitin = line.split("\t");
										String chr1 = splitin[0];
										int start1 = Integer.parseInt(splitin[1]);
										int end1 = Integer.parseInt(splitin[2]);
										String chr2 = splitin[3];
										int start2 = Integer.parseInt(splitin[4]);
										int end2 = Integer.parseInt(splitin[5]);
										if((!chr1.equalsIgnoreCase(chr2)) || Math.abs(start1-start2)<minDistEnds || Math.abs(end1-end2)<minDistEnds ||  Math.abs(start1-start2)>maxDistEnds || Math.abs(end1-end2)>maxDistEnds){
											continue;
										}
										IntervalTree<IntervalList> tree;
										
										if(intervalCollections.containsKey(chr1)){
											tree = intervalCollections.get(chr1);
										}else{
											tree = new IntervalTree<IntervalList>();
										}

										Interval interval1 = new Interval(chr1, start1, end1);
										Interval interval2 = new Interval(chr2, start2, end2);
										IntervalList list1 = null;
										if(tree.find(start1, end1) != null){
											list1 = tree.find(start1, end1).getValue();
											list1.add(interval2);
										}else{
											list1 = new IntervalList(wgsReader.getFileHeader());
											list1.add(interval2);
										}
										tree.put(start1, end1, list1);
										
										IntervalList list2 = null;
										if(tree.find(start2, end2) != null){
											list2 = tree.find(start2, end2).getValue();
											list2.add(interval1);
										}else{
											list2 = new IntervalList(wgsReader.getFileHeader());
											list2.add(interval1);
										}
										tree.put(start2, end2, list2);
										
										
										intervalCollections.put(chr1, tree);
										
										i++;
										
									}
									if(bedpe.endsWith(".gz")){
										gzipInputStream1.close();
									}
									br.close();
								
								
							}
							
							
							
							log.info("Parsing input cpgList bed file ...");
							
							
							HashMap<String,IntervalTree<String>> allCpgLocCollections = new HashMap<String,IntervalTree<String>>();
								log.info("Loading all CpG intervals ... ");
								GZIPInputStream gzipInputStream = null;
								BufferedReader br1;
								if(cpgFile.endsWith(".gz")){
									gzipInputStream = new GZIPInputStream(new FileInputStream(cpgFile));
									br1 = new BufferedReader(new InputStreamReader(gzipInputStream));
									
								}else{
									br1 = new BufferedReader(new FileReader(cpgFile));
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
										
										if(!intervalCollections.isEmpty()){
											if(intervalCollections.containsKey(chr)){
												IntervalTree<IntervalList> candidate = intervalCollections.get(chr);
												if(candidate.overlappers(start, end) == null){
													continue;
												}
											}else{
												continue;
											}
										}
										
										if(!regionsOverlapped.isEmpty()){
											if(regionsOverlapped.containsKey(chr)){
												if(regionsOverlapped.get(chr).overlappers(start-maxDistEnds > 0 ? start-maxDistEnds : 0, end+maxDistEnds) == null){
													continue;
												}
											}else{
												continue;
											}
													
							
										}
										
										
										IntervalTree<String> tree;
										
										if(allCpgLocCollections.containsKey(chr)){
											tree = allCpgLocCollections.get(chr);
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
										allCpgLocCollections.put(chr, tree);
										
										
									}
									if(cpgFile.endsWith(".gz")){
										gzipInputStream.close();
									}
									br1.close();
									
									
									HashMap<String, ArrayList<QueryInterval>> queryIntervalMap = new HashMap<String, ArrayList<QueryInterval>>();
									//ArrayList<QueryInterval> queryIntervalList = new ArrayList<QueryInterval>();
									
									//int z = 0;
									for(Entry<String,IntervalTree<String>> entry: allCpgLocCollections.entrySet()){
										String chr = entry.getKey();
										ArrayList<QueryInterval> queryIntervalList =  null;
										if(queryIntervalMap.containsKey(chr)){
											queryIntervalList = queryIntervalMap.get(chr);
										}else{
											queryIntervalList = new ArrayList<QueryInterval>();
										}
										int refIndex = dict.getSequenceIndex(entry.getKey());
										Iterator<Node<String>> it = entry.getValue().iterator();
										//System.err.println(entry.getKey() + "\t" + refIndex);
										while(it.hasNext()){
											Node<String> n = it.next();
											int start = n.getStart();
											int end = n.getEnd();
											//System.err.println(start + "\t" + end);
											queryIntervalList.add(new QueryInterval(refIndex, start+1, end));
										}
										queryIntervalMap.put(chr, queryIntervalList);
										//z++;
									}
									
									//QueryInterval[] queryIntervals = queryIntervalList.toArray(new QueryInterval[queryIntervalList.size()]);
									
									//for(QueryInterval interval : queryIntervals){
									//	System.err.println(interval);
									//}
									
					log.info("Extracting reads within intervals that contains cpg sites ... ");
					
					
					//System.err.println(cpgMethyMatrix.size() + "\t" + x + "\t" + y + "\t" + s + "\t" + a + "\t" + b + "\t" + c);
					log.info("Processing each pair of CpGs' interactin frequency ... ");
					writer.write("#chr\tstart1\tend1\tstart2\tend2\tname\tconcordCount\tdisconcordCount\tobserved_ratio\tpermutated_ratio\tcorrelation\tcorrelation_pvalue\t"
							+ "cor_random\tpvalue_random\tentropy\tentropy_random\tentropy_permuation_p_value\n");
					
					for(String chr : queryIntervalMap.keySet()){
						IntervalTree<HashMap<String,Integer>> cpgMethyMatrix = new IntervalTree<HashMap<String,Integer>>();
						ArrayList<QueryInterval> queryIntervalList =  queryIntervalMap.get(chr);
						QueryInterval[] queryIntervals = queryIntervalList.toArray(new QueryInterval[queryIntervalList.size()]);
						SAMRecordIterator wgsIt = wgsReader.queryOverlapping(queryIntervals);
						
						while(wgsIt.hasNext()){
							
							SAMRecord r = wgsIt.next();
							String readName = r.getReadName();
							
							int start = r.getAlignmentStart();
							int end = r.getAlignmentEnd();
							Iterator<IntervalTree.Node<String>> cpgs = null;
							if(allCpgLocCollections.containsKey(chr)){
								cpgs = allCpgLocCollections.get(chr).overlappers(start, end);
								if(cpgs == null){
								
									continue;
								}
								if(allCpgLocCollections.get(chr).minOverlapper(start, end)==null){
									
									continue;
								}else{
									
								}
								
							}else{
								
								continue;
							}
							
							if(failFlagFilter(r)){
								
								continue;
							}else{
								TreeMap<Integer, Integer> cpgSummary = BisulfiteHicUtils.readsMethySummaryAtEachLoc(r, minMapQ, minBaseQ, turnOffBisulfiteFilter, useBadMate);
								//System.err.print(r.getReadName() + "\t" + chr + "\t" + start + "\t" + end);
								//for(int key : cpgSummary.keySet()){
									//System.err.print("\t" + key);
								//}
								//System.err.println();
								while(cpgs.hasNext()){
									IntervalTree.Node<String> cpg = cpgs.next();
									int pos = cpg.getStart() - start;
									//System.err.println(cpg.getStart() + "\t" + pos);
									if(cpgSummary.containsKey(pos)){
										
										
											Node<HashMap<String,Integer>> tmpReadNode = cpgMethyMatrix.find(cpg.getStart(), cpg.getEnd());
											HashMap<String,Integer> tmpRead = null;
											if(tmpReadNode == null){
												tmpRead = new HashMap<String,Integer>();
												tmpRead.put(readName, cpgSummary.get(pos));
											}else{
												tmpRead = tmpReadNode.getValue();
												if(tmpRead.containsKey(readName)){
													// two ends overlapped and have not be trimmed?
												}else{
													tmpRead.put(readName, cpgSummary.get(pos));
												}
											}
											cpgMethyMatrix.put(cpg.getStart(), cpg.getEnd(), tmpRead);
										
										
										
									}
								}
							
							}

						}
						wgsIt.close();
						
						//IntervalTree<HashMap<String,Integer>> cpgs = cpgMethyMatrix.get(chr);
						List<Node<HashMap<String, Integer>>> cpgListTmp  = IteratorUtils.toList(cpgMethyMatrix.iterator()); 
						List<Node<HashMap<String, Integer>>> cpgList = new ArrayList<Node<HashMap<String, Integer>>>();
						for(int i = 0; i < cpgListTmp.size()-1; i++){
							Node<HashMap<String, Integer>> cpgi = cpgListTmp.get(i);
							if(cpgi.getValue().size() < minReadsInWindow){
								
								continue;
							}else{
								cpgList.add(cpgi);
							}
						}
						cpgListTmp.clear();
						
						IntervalTree<HashMap<String,Integer>> bgReadIntTree = new IntervalTree<HashMap<String,Integer>>();
						if(wgsBgBamReader != null){
							log.info("Extracting reads from background bam file within intervals that contains cpg sites ... ");
							for(int i = 0; i < cpgListTmp.size()-1; i++){
								Node<HashMap<String, Integer>> cpgi = cpgListTmp.get(i);
								HashMap<String,Integer> bgRead = new HashMap<String,Integer>();
								SAMRecordIterator wgsItBg = wgsReader.queryOverlapping(chr, cpgi.getStart(), cpgi.getEnd());
								while(wgsItBg.hasNext()){
									
									SAMRecord r = wgsItBg.next();
									String readName = r.getReadName();
									
									int start = r.getAlignmentStart();
									
									if(failFlagFilter(r)){
										
										continue;
									}else{
										TreeMap<Integer, Integer> cpgSummary = BisulfiteHicUtils.readsMethySummaryAtEachLoc(r, minMapQ, minBaseQ, turnOffBisulfiteFilter, useBadMate);
										int pos = cpgi.getStart() - start;
										bgRead.put(readName, cpgSummary.get(pos));
											
									}
									
								}
									wgsItBg.close();
									
									bgReadIntTree.put(cpgi.getStart(), cpgi.getEnd(), bgRead);
							}
							
						}
								
						
						
						//System.err.println(cpgList.size());
						for(int i = 0; i < cpgList.size()-1; i++){
							Node<HashMap<String, Integer>> cpgi = cpgList.get(i);
							if(cpgi.getValue().size() < minReadsInWindow){
								
								continue;
							}
							IntervalList list = null;
							Interval interval1 = null;
							if(!intervalCollections.isEmpty()){
								if(intervalCollections.containsKey(chr)){
									IntervalTree<IntervalList> candidate = intervalCollections.get(chr);
									if(candidate.overlappers(cpgi.getStart(), cpgi.getEnd()) != null){
										list = new IntervalList(samFileHeader);
										Iterator<Node<IntervalList>> listIt= candidate.overlappers(cpgi.getStart(), cpgi.getEnd());
										while(listIt.hasNext()){
											Node<IntervalList> tmp = listIt.next();
											if(interval1 == null){
												interval1 = new Interval(chr,tmp.getStart(), tmp.getEnd());
											}
											list.addall(tmp.getValue().getIntervals());
										}
										
									}else{
										continue;
									}
								}else{
									continue;
								}
							}
							Node<HashMap<String, Integer>> cpgi_bg = bgReadIntTree.find(cpgi.getStart(), cpgi.getEnd());
							
							for(int j = i+1; j < cpgList.size(); j++){
								Node<HashMap<String, Integer>> cpgj = cpgList.get(j);
								if(Math.abs(cpgi.getEnd() - cpgj.getEnd()) > maxDistEnds || Math.abs(cpgi.getEnd() - cpgj.getEnd()) < minDistEnds){
									
									continue;
								}
								if(!regionsOverlapped.isEmpty()){
									if(regionsOverlapped.containsKey(chr)){
										if(regionsOverlapped.get(chr).overlappers(cpgi.getStart(), cpgi.getEnd()) == null && regionsOverlapped.get(chr).overlappers(cpgj.getStart(), cpgj.getEnd()) == null){
											continue;
										}
									}else{
										continue;
									}
								}
								
								Interval interval2 = null;
								if(list != null){
									Interval intervalj = new Interval(chr, cpgj.getStart(), cpgj.getEnd());
									//System.err.println(intervalj);
									boolean intersected = false;
									for(Interval intervalTmp : list.getIntervals()){
										
										//if(cpgj.getStart()>intervalTmp.getStart() && cpgj.getStart()<intervalTmp.getEnd()){
											//System.err.println(intervalj + "\t" + intervalTmp + "\t" + intervalTmp.abuts(intervalj) + "\t" + intervalj.abuts(intervalTmp) + "\t" + intervalTmp.intersects(intervalj));
										//}
										if(intervalTmp.intersects(intervalj)){
											intersected = true;
											interval2 = intervalTmp.clone();
											break;
										}
									}
									//System.err.println();
									if(!intersected){
										
										continue;
									}
									
								}
								if(cpgj.getValue().size() < minReadsInWindow){
									continue;
								}
								Node<HashMap<String, Integer>> cpgj_bg = bgReadIntTree.find(cpgj.getStart(), cpgj.getEnd());
								ArrayList<Double> ratios = null;
								
								if(bgReadIntTree.size()>0){
									
									ratios = cpgContactFreqWithBgBam(cpgi.getValue(), cpgj.getValue(), cpgi_bg == null ? null : cpgi_bg.getValue(), cpgj_bg == null ? null : cpgj_bg.getValue());
								}else{
									ratios = cpgContactFreq(cpgi.getValue(), cpgj.getValue());
								}
									
								queryNum++;
								if(!ratios.get(2).isNaN()){
									//writer.write(chr + "\t" + cpgi.getStart() + "_" + cpgi.getEnd() + "\t" + cpgj.getStart() + "_" + cpgj.getEnd() + "\t" + ratios.getFirst() + "\t" + ratios.getSecond());
									String name = chr + ":" + cpgi.getStart() + ":" + cpgi.getEnd() + ":" + cpgj.getStart() + ":" + cpgj.getEnd();
									if(interval2 != null && interval1 != null){
										name = chr + ":" + interval1.getStart() + ":" + interval1.getEnd() + ":" + interval2.getStart() + ":" + interval2.getEnd(); //use bedpe's interval as row names...
									}
									writer.write(chr + "\t" + cpgi.getStart() + "\t" + cpgi.getEnd() + "\t" + cpgj.getStart() + "\t" + cpgj.getEnd() + "\t" + name);
									for(Double ratio : ratios){
										writer.write("\t" + ratio);
									}
									writer.write("\n");
								}else{
									//System.err.println(ratios.get(0));
									
								}
								if(queryNum % 1000000 == 0){
									log.info("Processing CpG pairs ... " + queryNum);
									writer.flush();
								}
							}
						}
						
					}
					wgsReader.close();
					if(bgBam != null){
						wgsBgBamReader.close();
					}
					//System.err.println(x + "\t" + y + "\t" + s + "\t" + a + "\t" + b + "\t" + c + "\t" + queryNum);
					finish();
			}
	
	//concordCount, disconcordCount, observed_ratio, permutated_ratio, correlation, p value, number of points, cor_random, pvalue_random, entropy, entropy_random, entropy_permuation_p_value
			private ArrayList<Double> cpgContactFreq(HashMap<String, Integer> x, HashMap<String, Integer> y){
				int concordCount = 0;
				int disconcordCount = 0;
				ArrayList<Double> results = new ArrayList<Double>();
				
				ArrayList<Double> xList = new ArrayList<Double>();
				ArrayList<Double> yList = new ArrayList<Double>();
				//for(String readName : x.keySet()){
					//System.err.println(readName);
				//}
				//for(String readName : y.keySet()){
				//	//System.err.println(readName);
				//}
				HashMap<String, Integer> entropyPattern = new HashMap<String, Integer>();
				for(String readName : x.keySet()){
					if(y.containsKey(readName)){
						if(x.get(readName) == y.get(readName)){
							concordCount++;
						}else{
							disconcordCount++;
						}
						xList.add((double)x.get(readName));
						yList.add((double)y.get(readName));
						String pattern = x.get(readName) + "" + y.get(readName);
						if(entropyPattern.containsKey(pattern)){
							entropyPattern.put(pattern, entropyPattern.get(pattern) + 1);
						}else{
							entropyPattern.put(pattern, 1);
						}
					}
				}
				int totalCount = concordCount + disconcordCount;
				if(totalCount == 0 || xList.size() < minReadsInWindow){
					results.add(0.);
					results.add(0.);
					
					results.add(Double.NaN);
					results.add(Double.NaN);
					
					results.add(Double.NaN);
					results.add(Double.NaN);
					//results.add(0.);
					
					results.add(Double.NaN);
					results.add(Double.NaN);
					
					results.add(Double.NaN);
					results.add(Double.NaN);
					results.add(Double.NaN);
					
					return results;
				}
				double observedRatio = Math.log10(concordCount + 0.0001) - Math.log10(disconcordCount + 0.0001 );
				//System.err.println(observedRatio);
				double[] xArray = ArrayUtils.toPrimitive(xList.toArray(new Double[xList.size()]));
				double[] yArray = ArrayUtils.toPrimitive(yList.toArray(new Double[yList.size()]));
				double[][] d = new double[xArray.length][2];
				for(int i = 0; i < xArray.length; i++){
					d[i][0] = xArray[i];
					d[i][1] = yArray[i];
				}
				PearsonsCorrelation corObserved = new PearsonsCorrelation(d);
				Double corEff = corObserved.correlation(xArray, yArray);
				Double pvalue = Double.NaN;
				if(!Double.isNaN(corEff)){
					pvalue = corObserved.getCorrelationPValues().getEntry(0, 1);
				}
				
				double entropyObserved = entropyNormalized(entropyPattern, totalCount);
				
				
				//permutation for the expectations:

				ArrayList<Double> corEffRandomList = new ArrayList<Double>();
				ArrayList<Double> pvalueRandomList = new ArrayList<Double>();
				
				ArrayList<Double> expectRatioList = new ArrayList<Double>();
				List<String> keysY  = new ArrayList<String>(y.keySet());
				List<String> keysX  = new ArrayList<String>(x.keySet());
				//ArrayList<Double> entropyRandomList = new ArrayList<Double>();
				ArrayList<Double> entropyRandomMeanList = new ArrayList<Double>();

				for(int i = 0; i < permutations; i++){
					int concordCountRandom = 0;
					int disconcordCountRandom = 0;
					
					ArrayList<Double> xListRandom = new ArrayList<Double>();
					ArrayList<Double> yListRandom = new ArrayList<Double>();
					
					HashMap<String, Integer> entropyPatternRandom = new HashMap<String, Integer>();
					
					for(int j = 0; j < totalCount; j++){
						String randomKey = keysY.get( generator.nextInt(keysY.size()));
						String randomKeyX = keysX.get( generator.nextInt(keysX.size()));
						
						if(x.get(randomKeyX) == y.get(randomKey)){
							concordCountRandom++;
						}else{
							disconcordCountRandom++;
						}
						xListRandom.add((double)x.get(randomKeyX));
						yListRandom.add((double)y.get(randomKey));
						String pattern = x.get(randomKeyX) + "" + y.get(randomKey);
						if(entropyPatternRandom.containsKey(pattern)){
							entropyPatternRandom.put(pattern, entropyPatternRandom.get(pattern) + 1);
						}else{
							entropyPatternRandom.put(pattern, 1);
						}
					}
					
					if(concordCountRandom + disconcordCountRandom == 0){
						continue;
					}
					double entropyTmp = entropyNormalized(entropyPatternRandom, concordCountRandom + disconcordCountRandom);
					entropyRandomMeanList.add(entropyTmp);
					expectRatioList.add(Math.log10(concordCountRandom + 0.0001) - Math.log10(disconcordCountRandom + 0.0001));

					double[] xArrayRandom = ArrayUtils.toPrimitive(xListRandom.toArray(new Double[xList.size()]));
					double[] yArrayRandom = ArrayUtils.toPrimitive(yListRandom.toArray(new Double[yList.size()]));
					double[][] dRandom = new double[xArrayRandom.length][2];
					for(int j = 0; j < xArrayRandom.length; j++){
						dRandom[j][0] = xArrayRandom[j];
						dRandom[j][1] = yArrayRandom[j];
					}
					PearsonsCorrelation corObservedRandom = new PearsonsCorrelation(dRandom);
					
					Double corEffRandomTmp = corObservedRandom.correlation(xArrayRandom, yArrayRandom);
					corEffRandomList.add(corEffRandomTmp);
					if(!Double.isNaN(corEffRandomTmp)){
						pvalueRandomList.add(corObservedRandom.getCorrelationPValues().getEntry(0, 1));
					}
					
				}
				
				double expectRatio = getMeanFromList(expectRatioList);
				double entropyRandomMean = getMeanFromList(entropyRandomMeanList);
				
				Collections.sort(entropyRandomMeanList);
				int pos = BisulfiteHicUtils.findPosition(entropyRandomMeanList, entropyObserved, 0, entropyRandomMeanList.size()-1);
				double entropyPermutatedPvalue = (double)(pos)/(double)entropyRandomMeanList.size();
				
				
				Double corEffRandom = getMeanFromList(corEffRandomList);
				Double pvalueRandom = getMeanFromList(pvalueRandomList);
				

				results.add((double)concordCount);
				results.add((double)disconcordCount);
				
				results.add(observedRatio);
				results.add(expectRatio);
				
				results.add(corEff);
				results.add(pvalue);
				//results.add((double)xArray.length);
				
				results.add(corEffRandom);
				results.add(pvalueRandom);
				
				results.add(entropyObserved);
				results.add(entropyRandomMean);
				results.add(entropyPermutatedPvalue);
				
				return results;
			}
	
			//concordCount, disconcordCount, observed_ratio, permutated_ratio, correlation, p value, number of points, cor_random, pvalue_random, entropy, entropy_random, entropy_permuation_p_value
			private ArrayList<Double> cpgContactFreqWithBgBam(HashMap<String, Integer> x, HashMap<String, Integer> y, HashMap<String, Integer> xb, HashMap<String, Integer> yb){
				int concordCount = 0;
				int disconcordCount = 0;
				ArrayList<Double> results = new ArrayList<Double>();
				
				ArrayList<Double> xList = new ArrayList<Double>();
				ArrayList<Double> yList = new ArrayList<Double>();
				//for(String readName : x.keySet()){
					//System.err.println(readName);
				//}
				//for(String readName : y.keySet()){
				//	//System.err.println(readName);
				//}
				HashMap<String, Integer> entropyPattern = new HashMap<String, Integer>();
				for(String readName : x.keySet()){
					if(y.containsKey(readName)){
						if(x.get(readName) == y.get(readName)){
							concordCount++;
						}else{
							disconcordCount++;
						}
						xList.add((double)x.get(readName));
						yList.add((double)y.get(readName));
						String pattern = x.get(readName) + "" + y.get(readName);
						if(entropyPattern.containsKey(pattern)){
							entropyPattern.put(pattern, entropyPattern.get(pattern) + 1);
						}else{
							entropyPattern.put(pattern, 1);
						}
					}
				}
				int totalCount = concordCount + disconcordCount;
				if(totalCount == 0 || xList.size() < minReadsInWindow){
					results.add(0.);
					results.add(0.);
					
					results.add(Double.NaN);
					results.add(Double.NaN);
					
					results.add(Double.NaN);
					results.add(Double.NaN);
					//results.add(0.);
					
					results.add(Double.NaN);
					results.add(Double.NaN);
					
					results.add(Double.NaN);
					results.add(Double.NaN);
					results.add(Double.NaN);
					
					return results;
				}
				double observedRatio = Math.log10(concordCount + 0.0001) - Math.log10(disconcordCount + 0.0001 );
				//System.err.println(observedRatio);
				double[] xArray = ArrayUtils.toPrimitive(xList.toArray(new Double[xList.size()]));
				double[] yArray = ArrayUtils.toPrimitive(yList.toArray(new Double[yList.size()]));
				double[][] d = new double[xArray.length][2];
				for(int i = 0; i < xArray.length; i++){
					d[i][0] = xArray[i];
					d[i][1] = yArray[i];
				}
				PearsonsCorrelation corObserved = new PearsonsCorrelation(d);
				Double corEff = corObserved.correlation(xArray, yArray);
				Double pvalue = Double.NaN;
				if(!Double.isNaN(corEff)){
					pvalue = corObserved.getCorrelationPValues().getEntry(0, 1);
				}
				
				double entropyObserved = entropyNormalized(entropyPattern, totalCount);
				
				
				if(xb==null || yb==null || xb.size() < minReadsInWindow || yb.size() < minReadsInWindow){
					results.add((double)concordCount);
					results.add((double)disconcordCount);
					
					results.add(observedRatio);
					results.add(Double.NaN);
					
					results.add(corEff);
					results.add(pvalue);
					//results.add(0.);
					
					results.add(Double.NaN);
					results.add(Double.NaN);
					
					results.add(entropyObserved);
					results.add(Double.NaN);
					results.add(Double.NaN);
					
					return results;
				}
				
				//permutation for the expectations:

				ArrayList<Double> corEffRandomList = new ArrayList<Double>();
				ArrayList<Double> pvalueRandomList = new ArrayList<Double>();
				
				ArrayList<Double> expectRatioList = new ArrayList<Double>();
				List<String> keysY  = new ArrayList<String>(y.keySet());
				List<String> keysX  = new ArrayList<String>(x.keySet());
				//ArrayList<Double> entropyRandomList = new ArrayList<Double>();
				ArrayList<Double> entropyRandomMeanList = new ArrayList<Double>();

				for(int i = 0; i < permutations; i++){
					int concordCountRandom = 0;
					int disconcordCountRandom = 0;
					
					ArrayList<Double> xListRandom = new ArrayList<Double>();
					ArrayList<Double> yListRandom = new ArrayList<Double>();
					
					HashMap<String, Integer> entropyPatternRandom = new HashMap<String, Integer>();
					
					for(int j = 0; j < totalCount; j++){
						String randomKey = keysY.get( generator.nextInt(yb.size()));
						String randomKeyX = keysX.get( generator.nextInt(xb.size()));
						
						if(xb.get(randomKeyX) == yb.get(randomKey)){
							concordCountRandom++;
						}else{
							disconcordCountRandom++;
						}
						xListRandom.add((double)xb.get(randomKeyX));
						yListRandom.add((double)yb.get(randomKey));
						String pattern = xb.get(randomKeyX) + "" + yb.get(randomKey);
						if(entropyPatternRandom.containsKey(pattern)){
							entropyPatternRandom.put(pattern, entropyPatternRandom.get(pattern) + 1);
						}else{
							entropyPatternRandom.put(pattern, 1);
						}
					}
					
					if(concordCountRandom + disconcordCountRandom == 0){
						continue;
					}
					double entropyTmp = entropyNormalized(entropyPatternRandom, concordCountRandom + disconcordCountRandom);
					entropyRandomMeanList.add(entropyTmp);
					expectRatioList.add(Math.log10(concordCountRandom + 0.0001) - Math.log10(disconcordCountRandom + 0.0001));

					double[] xArrayRandom = ArrayUtils.toPrimitive(xListRandom.toArray(new Double[xListRandom.size()]));
					double[] yArrayRandom = ArrayUtils.toPrimitive(yListRandom.toArray(new Double[yListRandom.size()]));
					double[][] dRandom = new double[xArrayRandom.length][2];
					for(int j = 0; j < xArrayRandom.length; j++){
						dRandom[j][0] = xArrayRandom[j];
						dRandom[j][1] = yArrayRandom[j];
					}
					PearsonsCorrelation corObservedRandom = new PearsonsCorrelation(dRandom);
					
					Double corEffRandomTmp = corObservedRandom.correlation(xArrayRandom, yArrayRandom);
					corEffRandomList.add(corEffRandomTmp);
					if(!Double.isNaN(corEffRandomTmp)){
						pvalueRandomList.add(corObservedRandom.getCorrelationPValues().getEntry(0, 1));
					}
					
				}
				
				double expectRatio = getMeanFromList(expectRatioList);
				double entropyRandomMean = getMeanFromList(entropyRandomMeanList);
				
				Collections.sort(entropyRandomMeanList);
				int pos = BisulfiteHicUtils.findPosition(entropyRandomMeanList, entropyObserved, 0, entropyRandomMeanList.size()-1);
				double entropyPermutatedPvalue = (double)(pos)/(double)entropyRandomMeanList.size();
				
				
				Double corEffRandom = getMeanFromList(corEffRandomList);
				Double pvalueRandom = getMeanFromList(pvalueRandomList);
				

				results.add((double)concordCount);
				results.add((double)disconcordCount);
				
				results.add(observedRatio);
				results.add(expectRatio);
				
				results.add(corEff);
				results.add(pvalue);
				//results.add((double)xArray.length);
				
				results.add(corEffRandom);
				results.add(pvalueRandom);
				
				results.add(entropyObserved);
				results.add(entropyRandomMean);
				results.add(entropyPermutatedPvalue);
				
				return results;
			}
			
			public double entropyNormalized(HashMap<String, Integer> counts, int total) {
				double entropy = 0;
				for (int count: counts.values()) {
					double prob = (double) count / total;
					entropy -= prob * Math.log(prob) / Math.log(2);
				}
				return entropy/(Math.log(total) / Math.log(2));
			}
			
			private boolean failFlagFilter(SAMRecord r){
				return r.getReadUnmappedFlag() || r.getNotPrimaryAlignmentFlag() || r.getMappingQuality() < minMapQ
						|| r.getReadFailsVendorQualityCheckFlag() || r.getDuplicateReadFlag() || !r.getReadPairedFlag() || !r.getProperPairFlag() ;
			}

			public double getMeanFromList(ArrayList<Double> list){
				Double sum = 0.;
				int count=0;
				for(Double a : list){
					if(!Double.isNaN(a)){
						count++;
						sum += a;
					}
				}
				if(count>0){
					return sum/count;
				}else{
					return Double.NaN;
				}
			}
			
			
			private void initiate(String outputPrefix) throws Exception{
				startTime = System.currentTimeMillis();
				writer = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputPrefix + ".singleBaseResHicMatrix.bed.gz")), "UTF-8");
				
				generator = new MersenneTwister();
			}

			private void finish() throws Exception{

				writer.close();

				long endTime   = System.currentTimeMillis();
				double totalTime = endTime - startTime;
				totalTime /= 1000;
				double totalTimeMins = totalTime/60;
				double totalTimeHours = totalTime/3600;
				
				log.info("SingleBaseResHicMatrix's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
			}
			

}
