/**
 *
 */
package main.java.edu.mit.compbio.bisulfitehic.utils;

import edu.unc.genomics.io.BigWigFileReader;
import htsjdk.samtools.*;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.exception.UnsortedFileException;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import main.java.edu.mit.compbio.ccinference.utils.CcInferenceUtils;
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.util.Pair;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * @author yaping
 *
 * only use the paired reads to to independence test
 * iterate all of the bam files, foreach bam files each paired reads, assign to ArrayList<ArrayList<Hashmap<String,ArrayList<Integer>>>>
 *     First array: different windows
 *     Second array: different alleles
 *     Hashmap: store by each read pairs
 *     Third array: methylation status
 */
public class TestIndepPoiBinomOnBam {

	@Option(name="-bams",usage="input bams for each single cell. Default: null")
	public ArrayList<String> bams = null;

	@Option(name="-minBaseQ",usage="minimum base quality score required to check. Default: 5")
	public int minBaseQ = 5;

	@Option(name="-minCellNum",usage="minimum number of cells. Default: 2")
	public int minCellNum = 2;

	@Option(name="-na",usage="String for NA. Default: NA")
	public String NA = "NA";

	@Option(name="-useBadMate",usage="use not proper paired reads. Default: false")
	public boolean useBadMate = false;

	@Option(name="-minMapQ",usage="minimum mapping quality score required to check. Default: 30")
	public int minMapQ = 30;

	@Option(name="-turnOffBisulfiteFilter",usage="turn off bisulfite conversion reads filter. Default: false")
	public boolean turnOffBisulfiteFilter = false;

	
	@Option(name="-aveWindow",usage="the window size to average the test statistics. default value do not use window to average. Default: 100000")
	public int aveWindow = 100000;

	@Option(name="-permutation",usage="permutation time to get p value. Default: 0")
	public int permutation = 0;

	@Option(name="-chr",usage="chromosome to detect. Default: chr19")
	public String chr = "chr19";

	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "TestIndepPoiBinomOnBam [opts] snp_vcf.gz outputPrefix -bams input1.bam -bams input2.bam ...";
	
	private static Logger log = Logger.getLogger(TestIndepPoiBinomOnBam.class);

	private OutputStreamWriter writer = null; 
	
	private static long startTime = -1;

	private MersenneTwister randomGenerator = null;
	
	/**
	 * @param args
	 */
	public static void main(String[] args)
		throws Exception {
			TestIndepPoiBinomOnBam tipbob = new TestIndepPoiBinomOnBam();
			BasicConfigurator.configure();
			tipbob.doMain(args);
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
						if(bams==null || bams.isEmpty())
								throw new CmdLineException(parser, USAGE, new Throwable());

						//read input bed file, for each row,
						//String intervalFile = arguments.get(0);
						String vcfFile = arguments.get(0);
						String outputPrefix = arguments.get(1);
						
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
								String chr_t = vc.getContig();
								//if(!chr_t.equalsIgnoreCase(chr)){
								//	continue;
								//}
								if(genotypeMapList.containsKey(chr_t)){
									tree = genotypeMapList.get(chr_t);
								}else{
									tree = new IntervalTree<Pair<Byte, Byte>>();
								}

								Byte alleleA = gt.getAllele(0).getBases()[0];
								Byte alleleB = gt.getAllele(1).getBases()[0];
								Pair<Byte, Byte> tmp = new Pair<Byte, Byte>(alleleA, alleleB);
								tree.put(vc.getStart(), vc.getEnd(), tmp);
								genotypeMapList.put(chr_t, tree);
							}


						}
						itSnp.close();
						vcfReader.close();


						log.info("Parsing input cpg methylation level from each of the single cell bam file ...");

						SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);

						/*
						*     First array: different windows
 *     Second array: different alleles: Ref allele, in bam_i+0,  Alt allele, in bam_i+1
 *     Hashmap: store by each read pairs
 *     Third array: methylation status

						 */
						ArrayList<ArrayList<HashMap<String,ArrayList<Integer>>>> methyStatus = new ArrayList<ArrayList<HashMap<String,ArrayList<Integer>>>>();
						long numCg = 0;
						long numCg_A = 0;
						long numCg_B = 0;
						for(int i = 0; i < bams.size(); i++) {
							String bamFile = bams.get(i);
							SamReader wgsReader = factory.open(new File(bamFile));
							SAMRecordIterator samIt = wgsReader.iterator();
							while (samIt.hasNext()) {
								SAMRecord r1 = samIt.next();
								SAMRecord r2 = samIt.next();
								int start1 = r1.getAlignmentStart();
								int end1 = r1.getAlignmentEnd();
								int start2 = r2.getAlignmentStart();
								int end2 = r2.getAlignmentEnd();

								String chr1 = r1.getContig();
								String chr2 = r2.getContig();
								if (!chr1.equalsIgnoreCase(chr2)) {
									continue;
								}
								//if(!chr1.equalsIgnoreCase(chr)){
								//	continue;
								//}

								String readname1 = r1.getReadName();
								String readname2 = r2.getReadName();
								if (!readname1.equalsIgnoreCase(readname2)) {
									throw new NumberFormatException("readname 1 is not equal to readname 2: " + readname1 + "\t" + readname2);
								}
								Pair<Byte, Byte> alleles = null;
								int snpStart = -1;
								//int snpEnd;
								boolean snpInRead2 = false;
								if (genotypeMapList.containsKey(chr1)) {
									IntervalTree<Pair<Byte, Byte>> cg = genotypeMapList.get(chr1);

									if (cg.minOverlapper(start1, end1) == null && cg.minOverlapper(start2, end2) == null) {
										continue;
									} else {
										if (cg.minOverlapper(start1, end1) != null) {
											IntervalTree.Node<Pair<Byte, Byte>> node = cg.minOverlapper(start1, end1);
											snpStart = node.getStart();
											//snpEnd = node.getEnd();
											alleles = node.getValue();
											snpInRead2 = false;
										} else if (cg.minOverlapper(start2, end2) != null) {
											IntervalTree.Node<Pair<Byte, Byte>> node = cg.minOverlapper(start2, end2);
											snpStart = node.getStart();
											//snpEnd = node.getEnd();
											alleles = node.getValue();
											snpInRead2 = true;
										}

									}
								} else {
									continue;
								}


								int r1pos = start1 / aveWindow;
								while (methyStatus.size() <= r1pos) {
									ArrayList<HashMap<String, ArrayList<Integer>>> alleleStatus = new ArrayList<HashMap<String, ArrayList<Integer>>>();
									methyStatus.add(alleleStatus);
								}
								ArrayList<HashMap<String, ArrayList<Integer>>> alleleStatusInRead1 = methyStatus.get(r1pos);

								int r2pos = start2 / aveWindow;
								while (methyStatus.size() <= r2pos) {
									ArrayList<HashMap<String, ArrayList<Integer>>> alleleStatus = new ArrayList<HashMap<String, ArrayList<Integer>>>();
									methyStatus.add(alleleStatus);
								}
								ArrayList<HashMap<String, ArrayList<Integer>>> alleleStatusInRead2 = methyStatus.get(r2pos);

									if(alleleStatusInRead1.size() <= bams.size()*2) {

										while (alleleStatusInRead1.size() <= bams.size()*2) {
											HashMap<String, ArrayList<Integer>> readStatus = new HashMap<String, ArrayList<Integer>>();
											alleleStatusInRead1.add(readStatus);
										}
										methyStatus.set(r1pos, alleleStatusInRead1);
									}
									if(alleleStatusInRead2.size() <= bams.size()*2) {

										while (alleleStatusInRead2.size() <= bams.size()*2) {
											HashMap<String, ArrayList<Integer>> readStatus = new HashMap<String, ArrayList<Integer>>();
											alleleStatusInRead2.add(readStatus);
										}
										methyStatus.set(r2pos, alleleStatusInRead2);
									}







								if (failFlagFilter(r1) || failFlagFilter(r2)) {
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
								if ((SequenceUtil.basesEqual(alleles.getFirst(), SequenceUtil.C) && SequenceUtil.basesEqual(alleles.getSecond(), SequenceUtil.T))
										|| (SequenceUtil.basesEqual(alleles.getFirst(), SequenceUtil.T) && SequenceUtil.basesEqual(alleles.getSecond(), SequenceUtil.C))
										) {
									if ((!secondPair && !negStrand) || (secondPair && negStrand)) {
										continue;
									}
								} else if ((SequenceUtil.basesEqual(alleles.getFirst(), SequenceUtil.G) && SequenceUtil.basesEqual(alleles.getSecond(), SequenceUtil.A))
										|| (SequenceUtil.basesEqual(alleles.getFirst(), SequenceUtil.A) && SequenceUtil.basesEqual(alleles.getSecond(), SequenceUtil.G))) {
									if ((!secondPair && negStrand) || (secondPair && !negStrand)) {
										continue;
									}
								}

								//int pos = loc.getStart() - d.getAlignmentStart();
								int pos = snpInRead2 ? (r2.getReadPositionAtReferencePosition(snpStart) - 1) : (r1.getReadPositionAtReferencePosition(snpStart) - 1);
								if (pos < 0) //no such a postion in the clipped reads
									continue;

								byte snpBase = CcInferenceUtils.toUpperCase(readBases[pos]);
								//if(negStrand){
								//snpBase = CcInferenceUtils.toUpperCase(SequenceUtil.complement(readBases[pos]));
								//}
								if (readBasesQ[pos] < minBaseQ) {
									continue;
								}

								TreeMap<Integer, Integer> ctSummary1 = BisulfiteHicUtils.readsMethySummaryAtEachLoc(r1, minMapQ, minBaseQ, turnOffBisulfiteFilter, useBadMate);
								TreeMap<Integer, Integer> ctSummary2 = BisulfiteHicUtils.readsMethySummaryAtEachLoc(r2, minMapQ, minBaseQ, turnOffBisulfiteFilter, useBadMate);

								if ((ctSummary1 == null || ctSummary1.isEmpty()) && (ctSummary2 == null || ctSummary2.isEmpty())) {
									continue;
								}


								HashMap<String, ArrayList<Integer>> read1CgInfo = new HashMap<String, ArrayList<Integer>>();
								ArrayList<Integer> read1CgList = new ArrayList<Integer>();
								read1CgList.addAll(ctSummary1.values());
								read1CgInfo.put(readname1, read1CgList);
								HashMap<String, ArrayList<Integer>> read2CgInfo = new HashMap<String, ArrayList<Integer>>();
								ArrayList<Integer> read2CgList = new ArrayList<Integer>();
								read2CgList.addAll(ctSummary2.values());
								read2CgInfo.put(readname2, read2CgList);
								numCg+=ctSummary1.size();
								numCg+=ctSummary2.size();
								if (BisulfiteHicUtils.baseMatchInBisulfiteSpace(alleles.getFirst(), snpBase, negStrand, secondPair)) {
									int allelePos = i + 0;

									HashMap<String, ArrayList<Integer>> readStatus1 = alleleStatusInRead1.get(allelePos);
									readStatus1.putAll(read1CgInfo);
									alleleStatusInRead1.set(allelePos, readStatus1);
									methyStatus.set(r1pos, alleleStatusInRead1);


									HashMap<String, ArrayList<Integer>> readStatus2 = alleleStatusInRead2.get(allelePos);
									readStatus2.putAll(read2CgInfo);
									alleleStatusInRead2.set(allelePos, readStatus2);
									methyStatus.set(r2pos, alleleStatusInRead2);
									numCg_A+=ctSummary1.size();
									numCg_A+=ctSummary2.size();

								} else if (BisulfiteHicUtils.baseMatchInBisulfiteSpace(alleles.getSecond(), snpBase, negStrand, secondPair)) {
									int allelePos = i + 1;

									HashMap<String, ArrayList<Integer>> readStatus1 = alleleStatusInRead1.get(allelePos);
									readStatus1.putAll(read1CgInfo);
									alleleStatusInRead1.set(allelePos, readStatus1);
									methyStatus.set(r1pos, alleleStatusInRead1);


									HashMap<String, ArrayList<Integer>> readStatus2 = alleleStatusInRead2.get(allelePos);
									readStatus2.putAll(read2CgInfo);
									alleleStatusInRead2.set(allelePos, readStatus2);
									methyStatus.set(r2pos, alleleStatusInRead2);
									numCg_B+=ctSummary1.size();
									numCg_B+=ctSummary2.size();
								}

							}
							samIt.close();
							wgsReader.close();
						}

			log.info("count cpg: " + numCg);
			log.info("count cpg_A: " + numCg_A);
			log.info("count cpg_B: " + numCg_B);

						log.info(methyStatus.size() + "\t" + methyStatus.get(0).size());
						log.info("Processing cpg pairs ...");
						long lineNum = 0;
						long effectivPairs = 0;
						Variance sd = new Variance() ;
						for(int i = 0; i < methyStatus.size(); i++) {//go over each window
							ArrayList<HashMap<String, ArrayList<Integer>>> alleleStatus_i = methyStatus.get(i);
							writer.write(chr + "\t" + i*aveWindow + "\t" + (i+1)*aveWindow );
							for(int j = 0; j < methyStatus.size(); j++) { //go over each window in the other end

								ArrayList<HashMap<String, ArrayList<Integer>>> alleleStatus_j = methyStatus.get(j);
								ArrayList<Pair<Integer, Integer>> cpg_i_d = new ArrayList<Pair<Integer, Integer>>();
								ArrayList<Pair<Integer, Integer>> cpg_j_d = new ArrayList<Pair<Integer, Integer>>();
								ArrayList<Double> cpg_i_j_d = new ArrayList<Double>();

								ArrayList<Double> cpg_i_j_random_sd_list = new  ArrayList<Double>();
								for (int z = 0; z < alleleStatus_i.size(); z++) { //go over each allele
									if(alleleStatus_i.isEmpty() || alleleStatus_j.isEmpty() ){
										continue;
									}
									HashMap<String, ArrayList<Integer>> readStatus_i = alleleStatus_i.get(z);
									HashMap<String, ArrayList<Integer>> readStatus_j = alleleStatus_j.get(z);

									if ((!readStatus_i.isEmpty()) && (!readStatus_j.isEmpty())) {
										//here: need to consider to put only paired reads in each of window or not?
										int sum_i = 0;
										int count_i = 0;
										int sum_j = 0;
										int count_j = 0;
										for (String readname : readStatus_i.keySet()) {
											ArrayList<Integer> read_i = readStatus_i.get(readname);
											if (readStatus_j.containsKey(readname)) {
												ArrayList<Integer> read_j = readStatus_j.get(readname);

												for (Integer t : read_i) {
													sum_i += t;
													count_i++;
												}
												for (Integer t : read_j) {
													sum_j += t;
													count_j++;
												}

											}
										}
										if (count_i + count_j > 0) {
											cpg_i_j_d.add((double)(sum_i + sum_j)/(double)(count_i + count_j));
										}
										if (count_i > 0) {
											cpg_i_d.add(new Pair<Integer,Integer>(sum_i, count_i));
										}
										if (count_j > 0) {
											cpg_j_d.add(new Pair<Integer,Integer>(sum_j, count_j));
										}

									}
									/*
									if (!readStatus_i.isEmpty()) {
										int sum = 0;
										int count = 0;
										for (ArrayList<Integer> read_i : readStatus_i.values()) {
											for (Integer t : read_i) {
												sum += t;
												count++;

											}
										}
										if (count > 0) {
											cpg_i_d.add(new Pair<Integer,Integer>(sum, count));
										}

									}
									if (!readStatus_j.isEmpty()) {
										int sum = 0;
										int count = 0;
										for (ArrayList<Integer> read_j : readStatus_j.values()) {
											for (Integer t : read_j) {
												sum += t;
												count++;
											}
										}
										if (count > 0) {
											cpg_j_d.add(new Pair<Integer,Integer>(sum, count));
										}
									}
									*/

								}

								if (cpg_i_j_d.size() <= minCellNum) {
									writer.write("\t" + NA);
									continue;
								}
								int s = 0;
								int size = Math.min(cpg_i_d.size(),cpg_j_d.size());
								while(s<permutation){
									Collections.shuffle(cpg_i_d);
									Collections.shuffle(cpg_j_d);
									double[] sd_tmp = new double[size];
									for(int z = 0; z < size; z++) {
										sd_tmp[z] = (cpg_i_d.get(z).getFirst() + cpg_j_d.get(z).getFirst())/(cpg_i_d.get(z).getSecond() + cpg_j_d.get(z).getSecond());
									}

									cpg_i_j_random_sd_list.add(sd.evaluate(sd_tmp));
									s++;
								}
								effectivPairs++;
								double obsStd = sd.evaluate(ArrayUtils.toPrimitive(cpg_i_j_d.toArray(new Double[cpg_i_j_d.size()])));

								int i_size = cpg_i_d.size();
								double[] cpg_i_sd = new double[i_size];
								for(int q = 0; q < i_size; q++){
									cpg_i_sd[q] = (double)cpg_i_d.get(q).getFirst()/(double)cpg_i_d.get(q).getSecond();
								}
								int j_size = cpg_j_d.size();
								double[] cpg_j_sd = new double[j_size];
								for(int q = 0; q < j_size; q++){
									cpg_j_sd[q] = (double)cpg_j_d.get(q).getFirst()/(double)cpg_j_d.get(q).getSecond();
								}

								double expStd = (sd.evaluate(cpg_i_sd) + sd.evaluate(cpg_j_sd))/2.0;

								if(permutation>0){
									Collections.sort(cpg_i_j_random_sd_list);
									int st_pos = BisulfiteHicUtils.findPosition(cpg_i_j_random_sd_list, obsStd, 0, cpg_i_j_random_sd_list.size()-1);
									//if(cpg_i_j_random_sd_list.size()-st_pos < st_pos){ //for two tails p
									//	st_pos = cpg_i_j_random_sd_list.size()-st_pos;
									//}
									double sdPermutatedPvalue = -Math.log10((double)(cpg_i_j_random_sd_list.size()-st_pos)/(double)cpg_i_j_random_sd_list.size());
									if(Double.isNaN(sdPermutatedPvalue)){
										writer.write("\t" + NA );
									}else{
										writer.write("\t" + sdPermutatedPvalue );
									}

								}else{
									if(Double.isNaN((obsStd/expStd))){
										writer.write("\t" + NA );
									}else{
										writer.write("\t" + (obsStd/expStd) );
									}
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

	private boolean failFlagFilter(SAMRecord r){
		return r.getReadUnmappedFlag() || r.getNotPrimaryAlignmentFlag() || r.getMappingQuality() < minMapQ
				|| r.getReadFailsVendorQualityCheckFlag() || r.getDuplicateReadFlag() || !r.getReadPairedFlag() || (useBadMate ? false : !r.getProperPairFlag());
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
			
			log.info("TestIndepPoiBinomOnBam's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
	}

}
