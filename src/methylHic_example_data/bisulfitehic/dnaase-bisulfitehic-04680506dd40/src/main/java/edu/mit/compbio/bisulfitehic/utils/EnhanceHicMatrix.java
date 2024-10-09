/**
 * enhanceHicMatrix.java
 * Nov 8, 2016
 * 9:26:51 AM
 * yaping    lyping1986@gmail.com
 */
package main.java.edu.mit.compbio.bisulfitehic.utils;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.util.Pair;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

/**
 * use the concordance of the methylation status at each end to increase the resolution of HiC. 
 * use simple correlation now... later the similarity of beta-binomial distriution
 * correlation, p value, number of points, cor_random, pvalue_random, numOfPoints_random
 */
public class EnhanceHicMatrix {

	
	@Option(name="-minBaseQ",usage="minimum base quality score required to check. Default: 5")
	public int minBaseQ = 5;
	
	@Option(name="-turnOffBisulfiteFilter",usage="turn off bisulfite conversion reads filter. Default: false")
	public boolean turnOffBisulfiteFilter = false;

	@Option(name="-useBadMate",usage="use not proper paired reads. Default: false")
	public boolean useBadMate = false;

	@Option(name="-useUnpaired",usage="reads are mapped in single end mode, have not joined into pairs yet. Default: false")
	public boolean useUnpaired = false;

	@Option(name="-minMapQ",usage="minimum mapping quality score required to check. Default: 30")
	public int minMapQ = 30;

	@Option(name="-minDistEnds",usage="minimum distance between two end to be output. Default: 20000")
	public int minDistEnds = 20000;

	@Option(name="-maxDistEnds",usage="maximum distance between two end to be output. Default: 2000000")
	public int maxDistEnds = 2000000;

	@Option(name="-resolution",usage="if not use -bedpe, the window size to look at matrix. Default: 25000")
	public int resolution = 25000;

	@Option(name="-window",usage="the smaller window size to look into the resolution. Default: 5000")
	public int window = 5000;

	@Option(name="-slidingWindow",usage="use 10% of window length to slide the window rather than non-overlapped one, default: not enabled")
	public boolean slidingWindow = false;
	
	@Option(name="-minReadsInWindow",usage="minimum number of reads required to calculate pearson correlation within the window at both ends Default: 3")
	public int minReadsInWindow   = 3;

	@Option(name="-bedpe",usage="bedpe file to limit the analysis only in these regions. Otherwise, it will scan the whole genome for all of the possible resolution windows. Default: null")
	public String bedpe = null;

	@Option(name="-chrom",usage="specify the chromosome. otherwise, use all the chrs Default: null")
	public String chrom = null;
	
	@Option(name="-minCgNum",usage="specify the minimum  number of CpG/GCH. Default: 1")
	public int minCgNum = 1;
	
	@Option(name="-permutations",usage="number of permutations to generate expectated p value within the same segments, default: 10")
	public int permutations = 10;
	
	@Option(name="-skipProcessFirstRow",usage="don't process first row for bedpe file, default: not enabled")
	public boolean skipProcessFirstRow = false;

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

	final private static String USAGE = "EnhanceHicMatrix [opts] outputPrefix hg19.fa.fai cpg_list.bed[.gz] input.bam";
	
	private static Logger log = Logger.getLogger(EnhanceHicMatrix.class);

	private OutputStreamWriter writer = null; 
	private MersenneTwister generator = null;
	
	private HashMap<String,Integer> chromSize = null;
	
	private static long startTime = -1;
	//private static long lineNum=0;
	private static long queryNum=0;

	/**
	 * @param args
	 */
	public static void main(String[] args)
		throws Exception {
			EnhanceHicMatrix ehm = new EnhanceHicMatrix();
				BasicConfigurator.configure();
				ehm.doMain(args);
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
							String outputPrefix = arguments.get(0);
							String chromSizeFile = arguments.get(1);
							String cpgFile = arguments.get(2);
							String bamFile = arguments.get(3);

							initiate(outputPrefix);
							log.info("Parsing chrom size file ...");
							getChromSize(chromSizeFile);
							
							ArrayList<Pair<Interval,Interval>> regions = getPairedRange();

							log.info("The number of interval pairs:" + regions.size());

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
										if(chrom != null && !chr.equalsIgnoreCase(chrom)){
											continue;
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
							
					log.info("Extracting reads within intervals that contains cpg sites ... ");
					SamReader wgsReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamFile));
					
					for(Pair<Interval,Interval> region : regions){
						//get each intervals' cpg
						Interval interval1 = region.getFirst();
						Interval interval2 = region.getSecond();
						//String chr = interval1.getContig();
						//Iterator<IntervalTree.Node<String>> cpgsIt1 = allCpgLocCollections.get(chr).overlappers(interval1.getStart(), interval1.getEnd());
						//Iterator<IntervalTree.Node<String>> cpgsIt2 = allCpgLocCollections.get(chr).overlappers(interval2.getStart(), interval2.getEnd());
						IntervalTree<HashMap<String,Double>> cpgMethyMatrix1 = new IntervalTree<HashMap<String,Double>>();
						IntervalTree<HashMap<String,Double>> cpgMethyMatrix2 = new IntervalTree<HashMap<String,Double>>();
						//DescriptiveStatistics summary1 = new DescriptiveStatistics();
						//DescriptiveStatistics summary2 = new DescriptiveStatistics();
						//System.err.println(interval1);
						//System.err.println(interval2);
						SAMRecordIterator wgsIt = wgsReader.queryOverlapping(interval1.getContig(), interval1.getStart(), interval1.getEnd());
						
						while(wgsIt.hasNext()){
							SAMRecord r = wgsIt.next();
							String readName = r.getReadName();
							String chr = r.getContig();
							int start = r.getAlignmentStart();
							int end = r.getAlignmentEnd();
							Iterator<IntervalTree.Node<String>> cpgs = null;
							if(allCpgLocCollections.containsKey(chr)){
								cpgs = allCpgLocCollections.get(chr).overlappers(start, end);
								if(cpgs == null){
									continue;
								}
								
							}else{
								continue;
							}
							
							Triple<Integer,Integer, Integer> cpgSummary = BisulfiteHicUtils.readsMethySummaryGeneral(r, methyPatternSearch, methyPatternCytPos, bsConvPatternSearch, bsConvPatternCytPos, minMapQ, minBaseQ, turnOffBisulfiteFilter, useBadMate);
							int numC = cpgSummary.getLeft();
							int numT = cpgSummary.getMiddle();
							int numCg = cpgSummary.getRight();
							
							if(numCg<minCgNum || numC+numT==0){
								continue;
							}
							
							double methy = (double)numC/(double)(numC+numT);
							
							
							HashMap<String,Double> tmpRead = null;
							if(cpgMethyMatrix1.find(start, end) == null){
								tmpRead = new HashMap<String,Double>();
								tmpRead.put(readName, methy);
							}else{
								tmpRead = cpgMethyMatrix1.find(start, end).getValue();
								if(tmpRead.containsKey(readName)){
									// two ends overlapped and have not be trimmed?
								}else{
									tmpRead.put(readName, methy);
								}
							}
							
							cpgMethyMatrix1.put(start, end, tmpRead);
						}
						wgsIt.close();
						
						wgsIt = wgsReader.queryOverlapping(interval2.getContig(), interval2.getStart(), interval2.getEnd());
						
						while(wgsIt.hasNext()){
							SAMRecord r = wgsIt.next();
							String readName = r.getReadName();
							String chr = r.getContig();
							int start = r.getAlignmentStart();
							int end = r.getAlignmentEnd();
							Iterator<IntervalTree.Node<String>> cpgs = null;
							if(allCpgLocCollections.containsKey(chr)){
								cpgs = allCpgLocCollections.get(chr).overlappers(start, end);
								if(cpgs == null){
									continue;
								}
								
							}else{
								continue;
							}
							
							Triple<Integer,Integer, Integer> cpgSummary = BisulfiteHicUtils.readsMethySummaryGeneral(r, methyPatternSearch, methyPatternCytPos, bsConvPatternSearch, bsConvPatternCytPos, minMapQ, minBaseQ, turnOffBisulfiteFilter, useBadMate);
							int numC = cpgSummary.getLeft();
							int numT = cpgSummary.getMiddle();
							int numCg = cpgSummary.getRight();
							
							if(numCg<minCgNum || numC+numT==0){
								continue;
							}
							
							double methy = (double)numC/(double)(numC+numT);
							
							
							HashMap<String,Double> tmpRead = null;
							if(cpgMethyMatrix2.find(start, end) == null){
								tmpRead = new HashMap<String,Double>();
								tmpRead.put(readName, methy);
							}else{
								tmpRead = cpgMethyMatrix2.find(start, end).getValue();
								if(tmpRead.containsKey(readName)){
									// two ends overlapped and have not be trimmed?
								}else{
									tmpRead.put(readName, methy);
								}
							}
							
							cpgMethyMatrix2.put(start, end, tmpRead);
						}
						wgsIt.close();
						//get each cpg's reads/methyStatus in each read
						//build beta-binomila distribution
						//System.err.println(cpgMethyMatrix1.size());
						//System.err.println(cpgMethyMatrix2.size());
						//use window size to scan all within the interval. calculate distance/correlation score within each window, for each window, also permuatte to calculate permutated distance & correlation score
						if(bedpe == null){ //dont' need to go over each possible window...

								int startW1 = interval1.getStart();
								int endW1 = interval1.getEnd();
								HashMap<String, Double> reads1 = new HashMap<String, Double>();
								Iterator<Node<HashMap<String, Double>>> readsIt1 = cpgMethyMatrix1.overlappers(startW1, endW1);
								if(readsIt1 == null){
									continue;
								}else{
									while(readsIt1.hasNext()){
										reads1.putAll(readsIt1.next().getValue());
									}

								}
								if(reads1.size() < minReadsInWindow){
									continue;
								}

									int startW2 = interval2.getStart();
									int endW2 = interval2.getEnd();
									HashMap<String, Double> reads2 = new HashMap<String, Double>();
									Iterator<Node<HashMap<String, Double>>> readsIt2 = cpgMethyMatrix2.overlappers(startW2, endW2);
									if(readsIt2 == null){
										continue;
									}else{
										while(readsIt2.hasNext()){
											reads2.putAll(readsIt2.next().getValue());
										}
									}
									if(reads2.size() < minReadsInWindow){
										continue;
									}
									//System.err.println(reads1.size());
									//System.err.println(reads2.size());

									ArrayList<Double> ratios = cpgContactFreq(reads1, reads2);
									queryNum++;
									if(!ratios.get(0).isNaN()){
										String name = interval1.getContig() + ":" + startW1 + ":" + endW1 + ":" + startW2 + ":" + endW2;
										writer.write(interval1.getContig() + "\t" + startW1 + "\t" + endW1 + "\t" + interval2.getContig() + "\t" + startW2 + "\t" + endW2 + "\t" + name);
										for(Double ratio : ratios){
											writer.write("\t" + ratio);
										}
										writer.write("\n");


									}
									if(queryNum % 1000000 == 0){
										log.info("Processing CpG ... " + queryNum);
										writer.flush();
									}



						}else{
							for(int startW1 = interval1.getStart(); startW1 < interval1.getEnd()-window; startW1+=(slidingWindow ? window/10 : window)){
								int endW1 = startW1 + window;
								HashMap<String, Double> reads1 = new HashMap<String, Double>();
								Iterator<Node<HashMap<String, Double>>> readsIt1 = cpgMethyMatrix1.overlappers(startW1, endW1);
								if(readsIt1 == null){
									continue;
								}else{
									while(readsIt1.hasNext()){
										reads1.putAll(readsIt1.next().getValue());
									}

								}
								if(reads1.size() < minReadsInWindow){
									continue;
								}
								for(int startW2 = interval2.getStart(); startW2 < interval2.getEnd()-window; startW2+=(slidingWindow ? window/10 : window)){
									int endW2 = startW2 + window;
									HashMap<String, Double> reads2 = new HashMap<String, Double>();
									Iterator<Node<HashMap<String, Double>>> readsIt2 = cpgMethyMatrix2.overlappers(startW2, endW2);
									if(readsIt2 == null){
										continue;
									}else{
										while(readsIt2.hasNext()){
											reads2.putAll(readsIt2.next().getValue());
										}
									}
									if(reads2.size() < minReadsInWindow){
										continue;
									}
									//System.err.println(reads1.size());
									//System.err.println(reads2.size());

									ArrayList<Double> ratios = cpgContactFreq(reads1, reads2);
									queryNum++;
									if(!ratios.get(0).isNaN()){
										String name = interval1.getContig() + ":" + startW1 + ":" + endW1 + ":" + startW2 + ":" + endW2;
										if(bedpe != null){
											name = interval1.getContig() + ":" + interval1.getStart() + ":" + interval1.getEnd() + ":" + interval2.getStart() + ":" + interval2.getEnd(); //use bedpe's interval as row names...
										}
										writer.write(interval1.getContig() + "\t" + startW1 + "\t" + endW1 + "\t" + interval2.getContig() + "\t" + startW2 + "\t" + endW2 + "\t" + name);
										for(Double ratio : ratios){
											writer.write("\t" + ratio);
										}
										writer.write("\n");


									}
									if(queryNum % 1000000 == 0){
										log.info("Processing CpG ... " + queryNum);
										writer.flush();
									}

								}
							}
						}

						
					}
					wgsReader.close();
					
					finish();
			}
	
	
	//correlation, p value, number of points, cor_random, pvalue_random, numOfPoints_random
	private ArrayList<Double> cpgContactFreq(HashMap<String, Double> x, HashMap<String, Double> y){
		
		ArrayList<Double> results = new ArrayList<Double>();
		
		ArrayList<Double> xList = new ArrayList<Double>();
		ArrayList<Double> yList = new ArrayList<Double>();
		for(String readName : x.keySet()){
			if(y.containsKey(readName)){
				xList.add(x.get(readName));
				yList.add(y.get(readName));

			}
		}
		if(xList.size() < minReadsInWindow){
			results.add(Double.NaN);
			results.add(Double.NaN);
			results.add(0.);
			results.add(Double.NaN);
			results.add(Double.NaN);
			results.add(0.);
			return results;
		}
		
		double[] xArray = ArrayUtils.toPrimitive(xList.toArray(new Double[xList.size()]));
		double[] yArray = ArrayUtils.toPrimitive(yList.toArray(new Double[yList.size()]));
		double[][] d = new double[xArray.length][2];
		for(int i = 0; i < xArray.length; i++){
			d[i][0] = xArray[i];
			d[i][1] = yArray[i];
		}
		//System.err.println(d.length + "\t" + d[0].length);
		//System.err.println(d[0][0] + "\t" + d[0][1]);
		//System.err.println(d[1][0] + "\t" + d[1][1]);
		//System.err.println(d[2][0] + "\t" + d[2][1]);
		PearsonsCorrelation corObserved = new PearsonsCorrelation(d);
		Double corEff = corObserved.correlation(xArray, yArray);
		Double pvalue = Double.NaN;
		if(!Double.isNaN(corEff)){
			pvalue = corObserved.getCorrelationPValues().getEntry(0, 1);
		}
		
		results.add(corEff);
		results.add(pvalue);
		results.add((double)xArray.length);
		
		if(Double.isNaN(corEff) || Double.isNaN(pvalue)){
			results.add(Double.NaN);
			results.add(Double.NaN);
			results.add(0.);
			return results;
		}

		

		//permutation for the expectations:
		ArrayList<Double> xListRandom = new ArrayList<Double>();
		ArrayList<Double> yListRandom = new ArrayList<Double>();
		
		List<String> keysY  = new ArrayList<String>(y.keySet());
		for(int i = 0; i < permutations; i++){
			
			
			for(String readName : x.keySet()){
				String randomKey = keysY.get( generator.nextInt(keysY.size()));
				xListRandom.add(x.get(readName));
				yListRandom.add(y.get(randomKey));
				
			}
		}
		double[] xArrayRandom = ArrayUtils.toPrimitive(xListRandom.toArray(new Double[xList.size()]));
		double[] yArrayRandom = ArrayUtils.toPrimitive(yListRandom.toArray(new Double[yList.size()]));
		double[][] dRandom = new double[xArrayRandom.length][2];
		for(int j = 0; j < xArrayRandom.length; j++){
			dRandom[j][0] = xArrayRandom[j];
			dRandom[j][1] = yArrayRandom[j];
		}
		PearsonsCorrelation corObservedRandom = new PearsonsCorrelation(dRandom);
		Double corEffRandom = corObservedRandom.correlation(xArrayRandom, yArrayRandom);
		Double pvalueRandom = Double.NaN;
		if(!Double.isNaN(corEffRandom)){
			pvalueRandom = corObservedRandom.getCorrelationPValues().getEntry(0, 1);
		}

		results.add(corEffRandom);
		results.add(pvalueRandom);
		results.add((double)xArrayRandom.length);
		
		//expectRatio /= permutations;

		return results;
	}
	
		private ArrayList<Pair<Interval,Interval>> getPairedRange() throws FileNotFoundException, IOException{
			ArrayList<Pair<Interval,Interval>> ranges = new ArrayList<Pair<Interval,Interval>>();
			if(bedpe != null){
				log.info("Parsing input bedpe bed file ...");
				HashMap<String,IntervalTree<String>> intervalCollections = new HashMap<String,IntervalTree<String>>();
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
						if((skipProcessFirstRow && i==0) || line.startsWith("#")){
							i++;
							continue;
						}else {
							String[] splitin = line.split("\t");
							String chr1 = splitin[0];
							int start1 = Integer.parseInt(splitin[1]);
							int end1 = Integer.parseInt(splitin[2]);
							String chr2 = splitin[3];
							int start2 = Integer.parseInt(splitin[4]);
							int end2 = Integer.parseInt(splitin[5]);
							if ((!chr1.equalsIgnoreCase(chr2)) || Math.abs(start1 - end2) < minDistEnds || Math.abs(end1 - start2) < minDistEnds || Math.abs(start1 - start2) > maxDistEnds || Math.abs(end1 - end2) > maxDistEnds) {
								continue;
							}
							if(chrom != null && !chr1.equalsIgnoreCase(chrom)){
								continue;
							}

							Interval interval1 = new Interval(chr1, start1, end1);
							Interval interval2 = new Interval(chr2, start2, end2);
							ranges.add(new Pair<Interval, Interval>(interval1, interval2));
							i++;
						}
						
					}
					if(bedpe.endsWith(".gz")){
						gzipInputStream1.close();
					}
					br.close();
				
				
			}else{
				log.info("Generating windows in the whole genome ...");
				for(String chr : chromSize.keySet()){
					if(chrom != null && !chr.equalsIgnoreCase(chrom)){
						continue;
					}
					int size = chromSize.get(chr);
					for(int start1 = 0; start1 < size - 2*resolution; start1+=resolution){
						int end1 = start1 + resolution;
						for(int start2 = end1; start2 < size - resolution; start2+=resolution){
							int end2 = start2 + resolution;
							if(Math.abs(start1-end2)<minDistEnds || Math.abs(end1-start2)<minDistEnds ||  Math.abs(start1-start2)>maxDistEnds || Math.abs(end1-end2)>maxDistEnds){
								continue;
							}
							Interval interval1 = new Interval(chr, start1, end1);
							Interval interval2 = new Interval(chr, start2, end2);
							ranges.add(new Pair<Interval,Interval>(interval1, interval2));
							
						}
					}
				}
				
			}
			
			return ranges;
		}
	
		private void getChromSize(String chromSizeFile) throws FileNotFoundException, IOException{
			GZIPInputStream gzipInputStream1 = null;
			BufferedReader br;
			if(chromSizeFile.endsWith(".gz")){
				gzipInputStream1 = new GZIPInputStream(new FileInputStream(chromSizeFile));
				br = new BufferedReader(new InputStreamReader(gzipInputStream1));
				
			}else{
				br = new BufferedReader(new FileReader(chromSizeFile));
			}
				
				String line;
				
				while( (line = br.readLine()) != null){
					if(line.startsWith("#"))
						continue;
					String[] splitLines = line.split("\t");
					if(splitLines.length<2){
						continue;
					}
					String chr = splitLines[0];
					int range = Integer.parseInt(splitLines[1]);
					chromSize.put(chr,range);
				}
				if(chromSizeFile.endsWith(".gz")){
					gzipInputStream1.close();
				}
				br.close();
				
		}
			
			
			private boolean failFlagFilter(SAMRecord r){
				return r.getReadUnmappedFlag() || r.getNotPrimaryAlignmentFlag() || r.getMappingQuality() < minMapQ
						|| r.getReadFailsVendorQualityCheckFlag() || r.getDuplicateReadFlag() || (useUnpaired ? false : !r.getReadPairedFlag()) || (useBadMate ? false : !r.getProperPairFlag());
			}

			
			
			
			private void initiate(String outputPrefix) throws Exception{
				startTime = System.currentTimeMillis();
				writer = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputPrefix + ".enhancedHicMatrix.txt.gz")), "UTF-8");
				 chromSize = new HashMap<String,Integer>();
				generator = new MersenneTwister();
			}

			private void finish() throws Exception{
				
				
				
				writer.close();
				
				
				long endTime   = System.currentTimeMillis();
				double totalTime = endTime - startTime;
				totalTime /= 1000;
				double totalTimeMins = totalTime/60;
				double totalTimeHours = totalTime/3600;
				
				log.info("EnhanceHicMatrix's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
			}
			

}
