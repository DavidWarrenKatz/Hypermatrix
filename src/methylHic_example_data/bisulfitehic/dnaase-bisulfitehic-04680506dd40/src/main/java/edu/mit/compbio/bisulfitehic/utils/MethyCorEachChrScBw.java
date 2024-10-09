/**
 * MethyCorAcrossHiccups.java
 * Sep 12, 2016
 * 1:03:54 PM
 * yaping    lyping1986@gmail.com
 */
package main.java.edu.mit.compbio.bisulfitehic.utils;

import edu.unc.genomics.io.BigWigFileReader;
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

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

public class MethyCorEachChrScBw {


	@Option(name="-bws",usage="input big wig files for each single cell for the methylation level or the methy count if -bwCounts is provided Default: null")
	public ArrayList<String> bws = null;

	//@Option(name="-bwCounts",usage="input big wig files for each single cell for the total C+T reads covered. Default: null")
	//public ArrayList<String> bwCounts = null;

	@Option(name="-aveWindow",usage="the window size to average the test statistics. default value do not use window to average. Default: 1000000")
	public int aveWindow = 1000000;

	@Option(name="-random",usage="number of times to shuffle the value. Default: 0")
	public int random = 0;

	@Option(name="-seed",usage="seed used for the randomness, default not use the seed. Default: -1")
	public int seed = -1;

	@Option(name="-shuffleWithinCell",usage="shuffle the value within the cell, between allele, the order of -bws should be AlleleA.bw, AlleleB.bw. Default: false")
	public boolean shuffleWithinCell = false;

	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "MethyCorEachChrScBw [opts] outputPrefix mm9.chroms.size chr";
	
	private static Logger log = Logger.getLogger(MethyCorEachChrScBw.class);

	private OutputStreamWriter methyCorWriter = null; 
	private Random generator = null;
	
	private static long startTime = -1;
	private static long lineNum=0;
	
	private List<BigWigFileReader> methyWigFileReaders = null;
	//private List<BigWigFileReader> covWigFileReaders = null;
	
	/**
	 * @param args
	 */
	public static void main(String[] args)
		throws Exception {
		MethyCorEachChrScBw mcah = new MethyCorEachChrScBw();
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
						//if( bws==null || bws.isEmpty() || (bwCounts != null && !bwCounts.isEmpty() && bwCounts.size() != bws.size()))
						if( bws==null || bws.isEmpty() )
							throw new CmdLineException(parser, USAGE, new Throwable());

						String outputPrefix = arguments.get(0);
						String chromSizeFile = arguments.get(1);
						String chr = arguments.get(2);
						initiate(outputPrefix);


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

										double cor = (random > 0 || shuffleWithinCell)? methyCorrelationRandom(chr, start1, end1, start2, end2, pearson)
												: methyCorrelation(chr, start1, end1, start2, end2, pearson);
										methyCorWriter.write(cor + "\t");
									}
									methyCorWriter.write("\n");
									methyCorWriter.flush();
								}
							}
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

		
		private double methyCorrelation(String chr, int start1, int end1, int start2, int end2,  PearsonsCorrelation pearson) throws Exception{

			double[] region1 = new double[bws.size()];
			double[] region2 = new double[bws.size()];
			//if(covWigFileReaders.isEmpty()){

				for(int i = 0; i < bws.size(); i++){
					BigWigFileReader methyReader = methyWigFileReaders.get(i);
					double methy1 = methyReader.queryStats(chr, start1, end1).getMean();
					double methy2 = methyReader.queryStats(chr, start2, end2).getMean();
					region1[i] = methy1;
					region2[i] = methy2;
				}

		/*
		}else{
				for(int i = 0; i < bws.size(); i++){
					BigWigFileReader methCountReader = methyWigFileReaders.get(i);
					BigWigFileReader covReader = covWigFileReaders.get(i);
					double methCount1 = methCountReader.queryStats(chr, start1, end1).getSum();
					double cov1 = covReader.queryStats(chr, start1, end1).getSum();
					double methCount2 = methCountReader.queryStats(chr, start2, end2).getSum();
					double cov2 = covReader.queryStats(chr, start2, end2).getSum();
					region1[i] = methCount1/cov1;
					region2[i] = methCount2/cov2;
				}
			}
		*/
			return pearson.correlation(region1,region2);
			
		}

	private double methyCorrelationRandom(String chr, int start1, int end1, int start2, int end2,  PearsonsCorrelation pearson) throws Exception{

		ArrayList<Double> region1 = new ArrayList<Double>();
		ArrayList<Double> region2 = new ArrayList<Double>();

		for(int i = 0; i < bws.size(); i++){
			BigWigFileReader methyReader = methyWigFileReaders.get(i);
			double methy1 = methyReader.queryStats(chr, start1, end1).getMean();
			double methy2 = methyReader.queryStats(chr, start2, end2).getMean();
			region1.add(methy1);
			region2.add(methy2);
		}


		double[] r1 = ArrayUtils.toPrimitive(region1.toArray(new Double[region1.size()]));
		if(shuffleWithinCell) {//should i randomly or always  shuffle between alleles??
			for(int i = 0; i < region2.size(); i+=2){
				Double tmp = region2.get(i);
				region2.set(i,region2.get(i+1));
				region2.set(i+1,tmp);
			}
			double[] r2 = ArrayUtils.toPrimitive(region2.toArray(new Double[region2.size()]));
			return pearson.correlation(r1,r2);

		}else {

			double cor = 0;
			int j = 0;
			while(j < random){

				Collections.shuffle(region2, generator);
				double[] r2 = ArrayUtils.toPrimitive(region2.toArray(new Double[region2.size()]));
				cor += pearson.correlation(r1,r2);
				j++;
			}


			return cor/random;
		}


	}


		
		private void initiate(String outputPrefix) throws Exception{
			startTime = System.currentTimeMillis();
				methyCorWriter = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputPrefix + ".methyCor.txt.gz")), "UTF-8");
			
			if(seed != -1){
				generator = new Random(seed);
			}else{
				generator = new Random();
			}

			
			methyWigFileReaders = new ArrayList<BigWigFileReader>();
			//covWigFileReaders = new ArrayList<BigWigFileReader>();
			if(bws != null){
				for(int i = 0; i < bws.size(); i++){
					BigWigFileReader bbreader = new BigWigFileReader((new File(bws.get(i))).toPath());
					methyWigFileReaders.add(bbreader);

				}
			}

/*
			if(bwCounts != null){
				for(int i = 0; i < bwCounts.size(); i++){

					BigWigFileReader bbreader2 = new BigWigFileReader((new File(bwCounts.get(i))).toPath());
					covWigFileReaders.add(bbreader2);
				}
			}
			*/
		}
		

		private void finish() throws Exception{
			
			if(methyWigFileReaders != null){
				for(int i = 0; i < methyWigFileReaders.size(); i++){
					methyWigFileReaders.get(i).close();
				}
			}
			/*
			if(covWigFileReaders != null){
				for(int i = 0; i < covWigFileReaders.size(); i++){
					covWigFileReaders.get(i).close();
				}
			}
			*/
			
				methyCorWriter.close();
			
			
			
			long endTime   = System.currentTimeMillis();
			double totalTime = endTime - startTime;
			totalTime /= 1000;
			double totalTimeMins = totalTime/60;
			double totalTimeHours = totalTime/3600;
			
			log.info("MethyCorEachChrScBw's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
		}
		

}
