/**
 * 
 */
package main.java.edu.mit.compbio.bisulfitehic.utils;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

/**
 * @author yaping
 *
 */
public class TestIndepPoiBinom {


	@Option(name="-maxDistance",usage="maximum distance of cpg to check. default is for the whole chromosomes Default: -1")
	public int maxDistance = -1;	
	
	@Option(name="-aveWindow",usage="the window size to average the test statistics. default value do not use window to average. Default: 1")
	public int aveWindow = 1;

	@Option(name="-permutation",usage="permutation time to get p value. Default: 0")
	public int permutation = 0;

	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "TestIndepPoiBinom [opts] input.bed.gz outputPrefix";
	
	private static Logger log = Logger.getLogger(TestIndepPoiBinom.class);

	private OutputStreamWriter writer = null; 
	
	private static long startTime = -1;

	private MersenneTwister randomGenerator = null;
	
	/**
	 * @param args
	 */
	public static void main(String[] args)
		throws Exception {
		TestIndepPoiBinom tipb = new TestIndepPoiBinom();
			BasicConfigurator.configure();
			tipb.doMain(args);
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
						

						//read input bed file, for each row,
						//String intervalFile = arguments.get(0);
						String inputFile = arguments.get(0);
						String outputPrefix = arguments.get(1);
						
						initiate(outputPrefix);
						log.info("Parsing input cpg methylation level bed file ...");
						GZIPInputStream gzipInputStream = null;
						BufferedReader br;
						if(inputFile.endsWith(".gz")){
							gzipInputStream = new GZIPInputStream(new FileInputStream(inputFile));
							br = new BufferedReader(new InputStreamReader(gzipInputStream));
							
						}else{
							br = new BufferedReader(new FileReader(inputFile));
						}
						String line;
						long lineNum = 0;
						//ArrayList<double[]> dataMat = new ArrayList<double[]>();
						ArrayList<Integer> cpgPos = new ArrayList<Integer>();
						//ArrayList<Double> cpgStd = new ArrayList<Double>();
						//Variance sd = new Variance() ;
						StandardDeviation sd = new StandardDeviation() ;
						int pos = 0;
						boolean restartFlag = true;
						//ArrayList<Integer> cpgPosTmp = new ArrayList<Integer>();
						//ArrayList<ArrayList<Double>> dataMatTmp = new ArrayList<ArrayList<Double>>();
						ArrayList<DescriptiveStatistics> statsByCpgs = new ArrayList<DescriptiveStatistics>(); //record for value in the same row (cpg)
						//cells * alleles
						ArrayList<ArrayList<DescriptiveStatistics>> statsByAlleles = new ArrayList<ArrayList<DescriptiveStatistics>>(); //record for value in the same allele (cell)
						DescriptiveStatistics statInEachCpgBlock = new DescriptiveStatistics();
						ArrayList<DescriptiveStatistics> statsByAllelesInEachCpgBlock = new ArrayList<DescriptiveStatistics>();

			//DescriptiveStatistics stats = new DescriptiveStatistics();

			//System.err.println(stats.getN()+ "\t" + stats.getStandardDeviation());

						while( (line = br.readLine()) != null){
							if(line.startsWith("#")){
								continue;
							}else{
								String[] splitin = line.split("\t");
								int cPos = Integer.parseInt(splitin[2]);
								if(aveWindow>1) {
									//System.err.println(pos + "\t" + cPos);
									while((pos+aveWindow) < cPos) {
										if(statsByAllelesInEachCpgBlock.isEmpty() || statsByAllelesInEachCpgBlock.size() < (splitin.length-6)) {
											pos+=aveWindow;
											continue;
											/*
											for(int i = 6, j = 0; i < splitin.length; i++, j++) {
												statInEachCpgBlock.addValue(Double.NaN);
												DescriptiveStatistics statInEachAllele = new DescriptiveStatistics();
												statInEachAllele.addValue(Double.NaN);
												statsByAllelesInEachCpgBlock.add(statInEachAllele);
											}
											*/
										}
										statsByCpgs.add(statInEachCpgBlock.copy());
										statsByAlleles.add((ArrayList<DescriptiveStatistics>) statsByAllelesInEachCpgBlock.clone());
										statsByAllelesInEachCpgBlock.clear();
										statInEachCpgBlock.clear();
										cpgPos.add(pos);

										pos+=aveWindow;
									}
								}
								
								for(int i = 6, j = 0; i < splitin.length; i++, j++) {
									if(splitin[i].equalsIgnoreCase("NA")) {
										if(aveWindow == 1){
											statInEachCpgBlock.addValue(Double.NaN);
										}

										if(statsByAllelesInEachCpgBlock.isEmpty() || statsByAllelesInEachCpgBlock.size() < (splitin.length-6)) {
											DescriptiveStatistics statInEachAllele = new DescriptiveStatistics();
											//statInEachAllele.addValue(Double.NaN);
											statsByAllelesInEachCpgBlock.add(statInEachAllele);
										}else {
											DescriptiveStatistics statInEachAllele = statsByAllelesInEachCpgBlock.get(j);
											//statInEachAllele.addValue(Double.NaN);
											statsByAllelesInEachCpgBlock.set(j, statInEachAllele);
										}
									}else {
										double v = Double.parseDouble(splitin[i]);
										statInEachCpgBlock.addValue(v);
										if(statsByAllelesInEachCpgBlock.isEmpty() || statsByAllelesInEachCpgBlock.size() < (splitin.length-6)) {
											DescriptiveStatistics statInEachAllele = new DescriptiveStatistics();
											statInEachAllele.addValue(v);
											statsByAllelesInEachCpgBlock.add(statInEachAllele);
										}else {
											DescriptiveStatistics statInEachAllele = statsByAllelesInEachCpgBlock.get(j);
											statInEachAllele.addValue(v);
											statsByAllelesInEachCpgBlock.set(j, statInEachAllele);
										}
										
									}
									
								}
								lineNum++;
								if(aveWindow>1) {
									//System.err.println(splitin.length + "\t" + statsByAllelesInEachCpgBlock.size());
									
								}else{ //if only on single cpg level, not bin data
									statsByCpgs.add(statInEachCpgBlock.copy());
									statInEachCpgBlock.clear();
									cpgPos.add(cPos);
								}
								
								if(lineNum % 1000000 == 0){
									log.info("Processing line: " + lineNum);
								}
								
							}
							
							
						}
						if(inputFile.endsWith(".gz")){
							gzipInputStream.close();
						}
						br.close();
						log.info(statsByCpgs.size() + "\t" + statsByAlleles.size());
						log.info("Processing cpg pairs ...");
						lineNum = 0;
						long effectivPairs = 0;
						for(int i = 0; i < statsByCpgs.size(); i++) {
							DescriptiveStatistics cpg_i_d = statsByCpgs.get(i);
							double[] cpg_i = cpg_i_d.getValues();
							double cpgStd_i = cpg_i_d.getStandardDeviation();
							int cpgPos_i = cpgPos.get(i);
							for(int j = 0; j < statsByCpgs.size(); j++) {
								int cpgPos_j = cpgPos.get(j);
								if(maxDistance > 0 && Math.abs(cpgPos_j-cpgPos_i) > maxDistance) {
									writer.write("NA" + "\t");
									continue;
								}
								DescriptiveStatistics cpg_j_d = statsByCpgs.get(j);
								double[] cpg_j = cpg_j_d.getValues();
								double cpgStd_j = cpg_j_d.getStandardDeviation();
								
								ArrayList<Double> cpg_i_j_list = new  ArrayList<Double>();
								ArrayList<Double> cpg_i_j_random_sd_list = new  ArrayList<Double>();
								if(aveWindow>1) {
									ArrayList<DescriptiveStatistics> statsByAllele_nonan_i = new ArrayList<DescriptiveStatistics>();
									ArrayList<DescriptiveStatistics> statsByAllele_nonan_j = new ArrayList<DescriptiveStatistics>();
									for(int z = 0; z < statsByAlleles.get(j).size(); z++) {
										//System.err.println(statsByAlleles.get(i).size() + "\t" + statsByAlleles.get(j).size());
										double cpg_i_mean = statsByAlleles.get(i).get(z).getMean();
										double cpg_j_mean = statsByAlleles.get(j).get(z).getMean();
										if((!Double.isNaN(cpg_i_mean)) && (!Double.isNaN(cpg_j_mean))) {
												cpg_i_j_list.add((cpg_i_mean + cpg_j_mean)/2.0);	
										}
										if(permutation>0){
											if(!Double.isNaN(cpg_i_mean)){
												statsByAllele_nonan_i.add(statsByAlleles.get(i).get(z));
											}
											if(!Double.isNaN(cpg_j_mean)){
												statsByAllele_nonan_j.add(statsByAlleles.get(j).get(z));
											}
										}
									}
									int s = 0;
									int size = Math.min(statsByAllele_nonan_i.size(),statsByAllele_nonan_j.size());
									while(s<permutation){
										Collections.shuffle(statsByAllele_nonan_i);
										Collections.shuffle(statsByAllele_nonan_j);
										double[] sd_tmp = new double[size];
										for(int z = 0; z < size; z++) {
											sd_tmp[z] = (statsByAllele_nonan_i.get(z).getMean() + statsByAllele_nonan_j.get(z).getMean())/2.0;
										}

										cpg_i_j_random_sd_list.add(sd.evaluate(sd_tmp));
										s++;
									}

								}else {
									ArrayList<Double> statsByAllele_nonan_i = new ArrayList<Double>();
									ArrayList<Double> statsByAllele_nonan_j = new ArrayList<Double>();
									for(int z = 0; z < cpg_j.length; z++) {
										if((!Double.isNaN(cpg_i[z])) && (!Double.isNaN(cpg_j[z]))) {
												cpg_i_j_list.add((cpg_i[z] + cpg_j[z])/2.0);	
										}
										if(permutation>0){
											if(!Double.isNaN(cpg_i[z])){
												statsByAllele_nonan_i.add(cpg_i[z]);
											}
											if(!Double.isNaN(cpg_j[z])){
												statsByAllele_nonan_j.add(cpg_j[z]);
											}
										}
									}
									int s = 0;
									int size = Math.min(statsByAllele_nonan_i.size(),statsByAllele_nonan_j.size());
									while(s<permutation){
										Collections.shuffle(statsByAllele_nonan_i);
										Collections.shuffle(statsByAllele_nonan_j);
										double[] sd_tmp = new double[size];
										for(int z = 0; z < size; z++) {
											sd_tmp[z] = (statsByAllele_nonan_i.get(z) + statsByAllele_nonan_j.get(z))/2.0;
										}

										cpg_i_j_random_sd_list.add(sd.evaluate(sd_tmp));
										s++;
									}
								}
								
								if(cpg_i_j_list.isEmpty() || cpg_i_j_list.size() <= 2) {
									writer.write("NA" + "\t");
									continue;
								}
								effectivPairs++;
								double obsStd = sd.evaluate(ArrayUtils.toPrimitive(cpg_i_j_list.toArray(new Double[cpg_i_j_list.size()])));
								double expStd = (cpgStd_i + cpgStd_j)/2.0;
								if(permutation>0){
									Collections.sort(cpg_i_j_random_sd_list);
									int st_pos = BisulfiteHicUtils.findPosition(cpg_i_j_random_sd_list, obsStd, 0, cpg_i_j_random_sd_list.size()-1);
									//if(cpg_i_j_random_sd_list.size()-st_pos < st_pos){ //for two tails p
									//	st_pos = cpg_i_j_random_sd_list.size()-st_pos;
									//}
									double sdPermutatedPvalue = -Math.log10((double)(cpg_i_j_random_sd_list.size()-st_pos)/(double)cpg_i_j_random_sd_list.size());
									writer.write(sdPermutatedPvalue + "\t");
								}else{
									writer.write((obsStd/expStd) + "\t");
								}

							}
							writer.write("\n");
							lineNum++;
							if(lineNum % 1000 == 0){
								log.info("Processing cpgs: " + lineNum + "\t Effective pairs: " + effectivPairs);
								writer.flush();
								
							}
						}
						
						finish();
						
		}
		
		private void initiate(String outputPrefix) throws Exception{
			startTime = System.currentTimeMillis();
			writer = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputPrefix + ".matrix.txt.gz")), "UTF-8");
			randomGenerator = new MersenneTwister();
		}
		

		private void finish() throws Exception{		
			writer.close();			
		
			long endTime   = System.currentTimeMillis();
			double totalTime = endTime - startTime;
			totalTime /= 1000;
			double totalTimeMins = totalTime/60;
			double totalTimeHours = totalTime/3600;
			
			log.info("TestIndepPoiBinom's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
		}

}
