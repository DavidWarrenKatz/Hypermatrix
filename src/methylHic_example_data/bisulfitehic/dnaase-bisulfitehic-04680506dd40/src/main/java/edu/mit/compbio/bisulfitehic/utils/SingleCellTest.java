package main.java.edu.mit.compbio.bisulfitehic.utils;

import htsjdk.samtools.util.IntervalTree;
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
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
import java.util.HashMap;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;


import JavaMI.MutualInformation;

/**
 * @author yaping
 *
 */
public class SingleCellTest {


    @Option(name="-aveWindow",usage="the window size to average the test statistics. default value do not use window to average. Default: 1")
    public int aveWindow = 1;

    @Option(name="-permutation",usage="permutation time to get p value. Default: 0")
    public int permutation = 0;

    @Option(name="-regions",usage="bed file for the regions to check. limited to be the same chromosome now... Default: null")
    public String regions = null;

    @Option(name="-scoreMode",usage="the mode to give score: 0. PearsonCorrelation; 1.PoissonBinomial; 2.MutualInformation. Default: 0")
    public int scoreMode = 0;

    @Option(name="-h",usage="show option information")
    public boolean help = false;


    @Argument
    private List<String> arguments = new ArrayList<String>();

    final private static String USAGE = "SingleCellTest [opts] input.bed.gz outputPrefix";

    private static Logger log = Logger.getLogger(SingleCellTest.class);

    private OutputStreamWriter writer = null;

    private static long startTime = -1;

    private MersenneTwister randomGenerator = null;

    /**
     * @param args
     */
    public static void main(String[] args)
            throws Exception {
        SingleCellTest sct = new SingleCellTest();
        BasicConfigurator.configure();
        sct.doMain(args);
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
            while( (line = br.readLine()) != null){
                if(line.startsWith("#"))
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

            }
            if(regions.endsWith(".gz")){
                gzipInputStream1.close();
            }
            br.close();
        }


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
        ArrayList<double[]> cpgMatrix = new ArrayList<double[]>();
        ArrayList<Integer> cpgLoci = new ArrayList<Integer>();

        while( (line = br.readLine()) != null){
            if(line.startsWith("#")){
                continue;
            }else{
                String[] splitin = line.split("\t");
                String chr = splitin[0];
                int start = Integer.parseInt(splitin[1]);
                int end = Integer.parseInt(splitin[2]);
                if(!regionsOverlapped.isEmpty()){
                    if(regionsOverlapped.containsKey(chr)){
                        if(regionsOverlapped.get(chr).overlappers(start, end) == null){
                            continue;
                        }
                    }else{
                        continue;
                    }


                }
                double[] row = new double[splitin.length-6];
                for(int i = 6, j = 0; i < splitin.length; i++, j++) {
                    if(splitin[i].equalsIgnoreCase("NA")) {
                        row[j] = Double.NaN;
                    }else {
                        row[j] = Double.parseDouble(splitin[i]);
                    }

                }
                cpgMatrix.add(row);
                cpgLoci.add(end);
                lineNum++;


                if(lineNum % 1000000 == 0){
                    log.info("Processing line: " + lineNum);
                }

            }


        }
        if(inputFile.endsWith(".gz")){
            gzipInputStream.close();
        }
        br.close();
        log.info(cpgMatrix.size());
        log.info("Processing cpg pairs ...");
        lineNum = 0;
        long effectivPairs = 0;

        Variance sd = new Variance();
        PearsonsCorrelation cor = new PearsonsCorrelation();

        for(int i = 0; i < cpgMatrix.size(); i++) {
            double[] cpg_i = cpgMatrix.get(i);
            int cpg_i_loc = cpgLoci.get(i);
            writer.write(cpg_i_loc);
            for(int j = 0; j < cpgMatrix.size(); j++) {
                double[] cpg_j = cpgMatrix.get(j);
                //int cpg_j_loc = cpgLoci.get(j);
                ArrayList<Double> cpg_i_nonan_list = new ArrayList<Double>();
                ArrayList<Double> cpg_j_nonan_list = new ArrayList<Double>();
                ArrayList<Double> cpg_i_j_nonan_list = new ArrayList<Double>();
                for(int z = 0; z < cpg_j.length; z++) {
                    if((!Double.isNaN(cpg_i[z])) && (!Double.isNaN(cpg_j[z]))) {
                        cpg_i_nonan_list.add(cpg_i[z]);
                        cpg_j_nonan_list.add(cpg_j[z]);
                        cpg_i_j_nonan_list.add((cpg_i[z]+cpg_j[z])/2.0);
                    }

                }
                if(cpg_i_nonan_list.isEmpty() || cpg_i_nonan_list.size() <= 2) {
                    writer.write("NA" + "\t");
                    continue;
                }
                double[] cpg_i_nonan = ArrayUtils.toPrimitive(cpg_i_nonan_list.toArray(new Double[cpg_i_nonan_list.size()]));
                double[] cpg_j_nonan = ArrayUtils.toPrimitive(cpg_j_nonan_list.toArray(new Double[cpg_j_nonan_list.size()]));
                double[] cpg_i_j_nonan = ArrayUtils.toPrimitive(cpg_i_j_nonan_list.toArray(new Double[cpg_i_j_nonan_list.size()]));


                effectivPairs++;
                double score = Double.NaN;
                if(scoreMode == 0){ //pearson correlation
                    score = cor.correlation(cpg_i_nonan, cpg_j_nonan);
                }else if(scoreMode == 1){ //poisson binomial test
                    double obsStd = sd.evaluate(cpg_i_j_nonan);
                    double expStd = (sd.evaluate(cpg_i) + sd.evaluate(cpg_j))/2.0;
                    score = obsStd/expStd;

                }else if(scoreMode == 2){ // mutual information
                    score = MutualInformation.calculateMutualInformation(cpg_i_nonan, cpg_j_nonan);
                }else{
                    throw new IllegalArgumentException("no such a mode");
                }
                writer.write("\t" + score);

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

        log.info("SingleCellTest's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
    }

}
