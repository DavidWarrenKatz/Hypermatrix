/**
 * Bhmem.java
 * Jul 21, 2016
 * 9:33:23 AM
 * yaping    lyping1986@gmail.com
 */
package main.java.edu.mit.compbio.bisulfitehic.mapping;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamPairUtil;
import htsjdk.samtools.util.CigarUtil;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.IllegalFormatException;
import java.util.List;
import java.util.zip.GZIPInputStream;

import main.java.edu.mit.compbio.bisulfitehic.utils.BisulfiteHicUtils;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.util.Pair;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import picard.sam.CreateSequenceDictionary;

import com.github.lindenb.jbwa.jni.BwaIndex;
import com.github.lindenb.jbwa.jni.BwaMem;
import com.github.lindenb.jbwa.jni.KSeq;
import com.github.lindenb.jbwa.jni.ShortRead;

/**
 *
 */
public class Bhmem {

	@Option(name="-clip5",usage="bp at 5' end to clip. default: 0")
	public int clip5 = 0;	
	@Option(name="-clip3",usage="bp at 3' end to clip. default: 0")
	public int clip3 = 0;	
	@Option(name="-buffer",usage="number of reads to buffer in memory. default: 100000")
	public int buffer = 100000;

	@Option(name="-rgId",usage="readgroup ID to add. default: null")
	public String rgId = null;	

	@Option(name="-rgSm",usage="readgroup SM to add. default: null")
	public String rgSm = null;

	@Option(name="-noSoftClipOnSuppl",usage="do not use softclip on supplemental reads. default: false")
	public boolean noSoftClipOnSuppl = false;
	
	@Option(name="-nonDirectional",usage="noDirectional protocol. default: false")
	public boolean noDirectional = false;
	
	@Option(name="-forceSameChr",usage="force to look at pairs within the same chr. default: false")
	public boolean forceSameChr = false;

	@Option(name="-outputMateDiffChr",usage="output mapping pairs mapped in different chrs. default: false")
	public boolean outputMateDiffChr = false;
	
	@Option(name="-enzymeList",usage="bed file about the location of restriction enzyme.it will be used to assist mapping also. Strongly recommmend to use it for scMethylHiC default: null")
	public String enzymeList = null;

	@Option(name="-readNameSuffix",usage="when read name in end1 and end2 are different, e.g. there is /1, /2 in the end of read name default: null")
	public String readNameSuffix = null;
	
	@Option(name="-pbat",usage="apply PBAT, which is for Joe Ecker's scWGBS protocol. default: false")
	public boolean pbat = false;
	
	
	
	@Option(name="-score",usage="minimum aligntment score to use. equal to -T in bwa mem. default: 0")
	public int score = 0;	
	@Option(name="-t",usage="number of CPU threads to use. default: 1")
	public int thread = 1;
	
	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "Bhmem [opts] dir/hg19.fa outputPrefix r1.fq[.gz] r2.fq[.gz]";
	
	private static Logger log = Logger.getLogger(Bhmem.class);

	private static long startTime = -1;
	private static long totalReads = 0L;
	private static long mappedReads = 0L;
	//private final static String ORIGINAL_SEQ_TAG = "OS";
	private final static String MISMATCHES = "NM";
	private final static String ALIGNMENT_SCORE = "AS";
	private static List<SamPairUtil.PairOrientation> ORITATION;
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		Bhmem bhmem = new Bhmem();
		BasicConfigurator.configure();
		bhmem.doMain(args);
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
					String refFile = arguments.get(0);
					String outputFile = arguments.get(1);
					String fastq1 = arguments.get(2);
					String fastq2 = arguments.get(3);

					initiate();	
					
					HashMap<String,IntervalTree<String>> regionsEnzyme = new HashMap<String,IntervalTree<String>>();
					if(enzymeList != null){
						log.info("Parsing input -enzymeList bed file ...");
						
						GZIPInputStream gzipInputStream1 = null;
						BufferedReader br;
						if(enzymeList.endsWith(".gz")){
							gzipInputStream1 = new GZIPInputStream(new FileInputStream(enzymeList));
							br = new BufferedReader(new InputStreamReader(gzipInputStream1));
							
						}else{
							br = new BufferedReader(new FileReader(enzymeList));
						}
							
							String line;
							long i = 0;
							while( (line = br.readLine()) != null){
								if(line.startsWith("#"))
									continue;
								String[] splitin = line.split("\t");
								String chr = splitin[0];
								int start = Integer.parseInt(splitin[1]);
								int end = Integer.parseInt(splitin[2]);
								
							
								IntervalTree<String> tree;
								
								if(regionsEnzyme.containsKey(chr)){
									tree = regionsEnzyme.get(chr);
								}else{
									tree = new IntervalTree<String>();
								}
								
								tree.put(start, end, "1");
								
								
								i++;
								
							}
							if(enzymeList.endsWith(".gz")){
								gzipInputStream1.close();
							}
							br.close();
					}
					
					
					SAMFileHeader samFileHeader = new SAMFileHeader();
					File sequenceDict = new File(refFile.replaceAll(".fa", ".dict"));
					if(sequenceDict.exists() && !sequenceDict.isDirectory()){
						samFileHeader.setSequenceDictionary(SAMSequenceDictionaryExtractor.extractDictionary(sequenceDict));
					}else{
						CreateSequenceDictionary csd = new CreateSequenceDictionary();
						samFileHeader.setSequenceDictionary(csd.makeSequenceDictionary(sequenceDict));
					}
					if(rgId != null && !rgId.isEmpty() && rgSm != null && !rgSm.isEmpty()){
						SAMReadGroupRecord rg = new SAMReadGroupRecord(rgId);
						rg.setLibrary(rgId);
						rg.setSample(rgSm);
						rg.setPlatform("illumina");
						rg.setPlatformUnit("barcode");
						ArrayList<SAMReadGroupRecord> rgList = new ArrayList<SAMReadGroupRecord>();
						rgList.add(rg);
						samFileHeader.setReadGroups(rgList);
					}
					
					SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
					SAMFileWriter writer = writerFactory.makeBAMWriter(samFileHeader, true, new File(outputFile));
					
					BwaIndex indexCT=new BwaIndex(new File(new File(refFile).getParent() + "/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa"));
					BwaMem memCT=new BwaMem(indexCT);
					memCT.updateMoreScoringParameters(thread, score, noSoftClipOnSuppl ? 0 : 1, clip5, clip3); //1 means use softclip supplement mapping. 0 means NO softclip supplement mapping.
					
					BwaIndex indexGA=new BwaIndex(new File(new File(refFile).getParent() + "/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"));
					BwaMem memGA=new BwaMem(indexGA);
					memGA.updateMoreScoringParameters(thread, score, noSoftClipOnSuppl ? 0 : 1, clip5, clip3); //1 means use softclip supplement mapping. 0 means NO softclip supplement mapping.

					
					KSeq kseq1=new KSeq(new File(fastq1));
					KSeq kseq2=new KSeq(new File(fastq2));
				
					List<String> L1=new ArrayList<String>();
					List<String> L2=new ArrayList<String>();
					
					List<ShortRead> L1CT=new ArrayList<ShortRead>();
					List<ShortRead> L2GA=new ArrayList<ShortRead>();
					
					for(;;){
						ShortRead read1=kseq1.next();
						ShortRead read2=kseq2.next();
						
						
						if(read1==null || read2==null || L1.size()>buffer){
							if(!L1.isEmpty()){
								String[] samsEnd1CT = pbat ? memGA.align(L1CT) : memCT.align(L1CT);
								String[] samsEnd1GA = pbat ? memCT.align(L1CT) : memGA.align(L1CT);
								String[] samsEnd2CT = pbat ? memGA.align(L2GA) : memCT.align(L2GA);
								String[] samsEnd2GA = pbat ? memCT.align(L2GA) : memGA.align(L2GA);
								HashMap<String, Pair<SAMRecord, SAMRecord>> bestPairs = joinTwoEnds(samsEnd1CT, samsEnd1GA, samsEnd2CT, samsEnd2GA, samFileHeader, L1, L2, regionsEnzyme) ;
								for(Pair<SAMRecord, SAMRecord> pair : bestPairs.values()){
									writer.addAlignment(pair.getFirst());
									writer.addAlignment(pair.getSecond());
									mappedReads+=2;
								}
							}
								
							if(read1==null || read2==null){
								break;
							}
							L1.clear();
							L2.clear();
							L1CT.clear();
							L2GA.clear();
						}
						String seq1 = new String(read1.getBases());
						String seq2 = new String(read2.getBases());
						//if(R1PBAT){
						//	seq1 = SequenceUtil.reverseComplement(seq1);
						//}
						ShortRead read1ct= pbat ? new ShortRead(read1.getName(), seq1.replace('G', 'A').getBytes() ,read1.getQualities()) : new ShortRead(read1.getName(), seq1.replace('C', 'T').getBytes() ,read1.getQualities());
						ShortRead read2ga= pbat ? new ShortRead(read2.getName(), seq2.replace('C', 'T').getBytes() ,read2.getQualities()) : new ShortRead(read2.getName(), seq2.replace('G', 'A').getBytes() ,read2.getQualities());
						
						L1.add(seq1);
						L2.add(seq2);
						L1CT.add(read1ct);
						L2GA.add(read2ga);
						totalReads+=2;
						if(totalReads % buffer == 0){
							log.info("Processing reads " + totalReads + " ...");
							log.info("Unique mapped and paired reads " + mappedReads + " ...");
							
						}
						
					}
					kseq1.dispose();
					kseq2.dispose();
					indexCT.close();
					memCT.dispose();
					indexGA.close();
					memGA.dispose();
					writer.close();
					
					finish();
	}
	
	private HashMap<String, Pair<SAMRecord, SAMRecord>> joinTwoEnds(String[] samsEnd1CT, String[] samsEnd1GA, String[] samsEnd2CT, String[] samsEnd2GA, SAMFileHeader samFileHeader, List<String> L1, List<String> L2){
		HashMap<String, Pair<SAMRecord, SAMRecord>> bestPairs = new HashMap<String, Pair<SAMRecord, SAMRecord>>();
		
		HashMap<String, SAMRecord> samsEnd1CTHash = new HashMap<String, SAMRecord>();
		HashMap<String, SAMRecord> samsEnd1GAHash = new HashMap<String, SAMRecord>();
		HashMap<String, SAMRecord> samsEnd2CTHash = new HashMap<String, SAMRecord>();
		HashMap<String, SAMRecord> samsEnd2GAHash = new HashMap<String, SAMRecord>();
		HashSet<String> readsNameCollection = new HashSet<String>();
		
		//make a hashmap about reads name
		for(int i = 0; i < samsEnd1CT.length; i++){
			for(SAMRecord r : String2SamRecord(samsEnd1CT[i], samFileHeader, L1.get(i), false)){
				if(!failFlagFilter(r)){
					String key = forceSameChr ? (r.getReadName() + "\t" + r.getReferenceName()) : r.getReadName();
					if(samsEnd1CTHash.containsKey(key)){
						samsEnd1CTHash.put(key, comparingSamRecord(samsEnd1CTHash.get(key), r));
					}else{
						readsNameCollection.add(key);
						samsEnd1CTHash.put(key, r);
					}
				}
			}
		}
		for(int i = 0; i < samsEnd1GA.length; i++){
			for(SAMRecord r : String2SamRecord(samsEnd1GA[i], samFileHeader, L1.get(i), false)){
				if(!failFlagFilter(r)){
					String key = forceSameChr ? (r.getReadName() + "\t" + r.getReferenceName()) : r.getReadName();
					if(samsEnd1GAHash.containsKey(key)){
						samsEnd1GAHash.put(key, comparingSamRecord(samsEnd1GAHash.get(key), r));
					}else{
						readsNameCollection.add(key);
						samsEnd1GAHash.put(key, r);
					}
				}
			}
		}
		for(int i = 0; i < samsEnd2CT.length; i++){
			for(SAMRecord r : String2SamRecord(samsEnd2CT[i], samFileHeader, L2.get(i), true)){
				if(!failFlagFilter(r)){
					String key = forceSameChr ? (r.getReadName() + "\t" + r.getReferenceName()) : r.getReadName();
					if(samsEnd2CTHash.containsKey(key)){
						samsEnd2CTHash.put(key, comparingSamRecord(samsEnd2CTHash.get(key), r));
					}else{
						readsNameCollection.add(key);
						samsEnd2CTHash.put(key, r);
					}
				}
			}
		}
		for(int i = 0; i < samsEnd2GA.length; i++){
			for(SAMRecord r : String2SamRecord(samsEnd2GA[i], samFileHeader, L2.get(i), true)){
				if(!failFlagFilter(r)){
					String key = forceSameChr ? (r.getReadName() + "\t" + r.getReferenceName()) : r.getReadName();
					if(samsEnd2GAHash.containsKey(key)){
						samsEnd2GAHash.put(key, comparingSamRecord(samsEnd2GAHash.get(key), r));
					}else{
						readsNameCollection.add(key);
						samsEnd2GAHash.put(key, r);
					}
				}
			}
		}
		
		
		//for each reads name, find the best pair to represent it.
		for(String name : readsNameCollection){
			String[] splitin = name.split("\\t");
			String readname = splitin[0];
			if(samsEnd1CTHash.containsKey(name) && samsEnd2CTHash.containsKey(name)){
				SAMRecord r1 = samsEnd1CTHash.get(name);
				SAMRecord r2 = samsEnd2CTHash.get(name);
				if((r1.getReadNegativeStrandFlag() && !r2.getReadNegativeStrandFlag()) || (r2.getReadNegativeStrandFlag() && !r1.getReadNegativeStrandFlag())){
					if(bestPairs.containsKey(readname)){
						bestPairs.put(readname, comparingSamRecord(bestPairs.get(readname), new Pair<SAMRecord, SAMRecord>(r1, r2)));
					}else{
						bestPairs.put(readname, new Pair<SAMRecord, SAMRecord>(r1, r2));
					}
				}
				
			}
			if(samsEnd1GAHash.containsKey(name) && samsEnd2GAHash.containsKey(name)){
				SAMRecord r1 = samsEnd1GAHash.get(name);
				SAMRecord r2 = samsEnd2GAHash.get(name);
				if((r1.getReadNegativeStrandFlag() && !r2.getReadNegativeStrandFlag()) || (r2.getReadNegativeStrandFlag() && !r1.getReadNegativeStrandFlag())){
					if(bestPairs.containsKey(readname)){
						bestPairs.put(readname, comparingSamRecord(bestPairs.get(readname), new Pair<SAMRecord, SAMRecord>(r1, r2)));
					}else{
						bestPairs.put(readname, new Pair<SAMRecord, SAMRecord>(r1, r2));
					}
				}
				
			}
			
		}
		
		//foreach best pair, compute insertion size, fix the flag to be: paired and proper paired, it is 1st/2nd end, its mate is reverse or not..
		for(String name : bestPairs.keySet()){
			SAMRecord r1 = bestPairs.get(name).getFirst();
			SAMRecord r2 = bestPairs.get(name).getSecond();
			//System.err.println(r1.getFlags() + "\t" + r1.getAlignmentStart() + "\t" + r1.getMateReferenceName() + "\t" + r1.getMateAlignmentStart() + "\t" + r1.getInferredInsertSize());
			//System.err.println(r2.getFlags() + "\t" + r2.getAlignmentStart() + "\t" + r2.getMateReferenceName() + "\t" + r2.getMateAlignmentStart() + "\t" + r2.getInferredInsertSize());
			SamPairUtil.setProperPairAndMateInfo(r1, r2, samFileHeader, ORITATION, true);
			//System.err.println(r1.getFlags() + "\t" + r1.getAlignmentStart() + "\t" + r1.getMateReferenceName() + "\t" + r1.getMateAlignmentStart() + "\t" + r1.getInferredInsertSize());
			//System.err.println(r2.getFlags() + "\t" + r2.getAlignmentStart() + "\t" + r2.getMateReferenceName() + "\t" + r2.getMateAlignmentStart() + "\t" + r2.getInferredInsertSize());
			bestPairs.put(name, new Pair<SAMRecord, SAMRecord>(r1, r2));
		}
		
		return bestPairs;
	}
	
	
	//PBAT will increase the chance to get chimeric reads, so it is even more difficulty to get the reads aligned correctly..
	private HashMap<String, Pair<SAMRecord, SAMRecord>> joinTwoEnds(String[] samsEnd1CT, String[] samsEnd1GA, String[] samsEnd2CT, String[] samsEnd2GA, SAMFileHeader samFileHeader, List<String> L1, List<String> L2, HashMap<String,IntervalTree<String>> regionsEnzyme){
		HashMap<String, Pair<SAMRecord, SAMRecord>> bestPairs = new HashMap<String, Pair<SAMRecord, SAMRecord>>();
		
		HashMap<String, ArrayList<SAMRecord>> samsEnd1CTHash = new HashMap<String, ArrayList<SAMRecord>>();
		HashMap<String, ArrayList<SAMRecord>> samsEnd1GAHash = new HashMap<String, ArrayList<SAMRecord>>();
		HashMap<String, ArrayList<SAMRecord>> samsEnd2CTHash = new HashMap<String, ArrayList<SAMRecord>>();
		HashMap<String, ArrayList<SAMRecord>> samsEnd2GAHash = new HashMap<String, ArrayList<SAMRecord>>();
		HashSet<String> readsNameCollection = new HashSet<String>();
		
		//make a hashmap about reads name
		for(int i = 0; i < samsEnd1CT.length; i++){
			for(SAMRecord r : String2SamRecord(samsEnd1CT[i], samFileHeader, L1.get(i), false)){
				if(!failFlagFilter(r)){
					String key = forceSameChr ? (r.getReadName() + "\t" + r.getReferenceName()) : r.getReadName();
					
					if(samsEnd1CTHash.containsKey(key)){
						ArrayList<SAMRecord> newCollection = samsEnd1CTHash.get(key);
						newCollection.add(r);
						samsEnd1CTHash.put(key, newCollection);
					}else{
						readsNameCollection.add(key);
						ArrayList<SAMRecord> newCollection = new ArrayList<SAMRecord>();
						newCollection.add(r);
						samsEnd1CTHash.put(key, newCollection);
					}
				}
			}
		}
		for(int i = 0; i < samsEnd1GA.length; i++){
			for(SAMRecord r : String2SamRecord(samsEnd1GA[i], samFileHeader, L1.get(i), false)){
				if(!failFlagFilter(r)){
					String key = forceSameChr ? (r.getReadName() + "\t" + r.getReferenceName()) : r.getReadName();
					if(samsEnd1GAHash.containsKey(key)){
						ArrayList<SAMRecord> newCollection = samsEnd1GAHash.get(key);
						newCollection.add(r);
						samsEnd1GAHash.put(key, newCollection);
					}else{
						readsNameCollection.add(key);
						ArrayList<SAMRecord> newCollection = new ArrayList<SAMRecord>();
						newCollection.add(r);
						samsEnd1GAHash.put(key, newCollection);
					}
				}
			}
		}
		
		for(int i = 0; i < samsEnd2CT.length; i++){
			for(SAMRecord r : String2SamRecord(samsEnd2CT[i], samFileHeader, L2.get(i), true)){
				if(!failFlagFilter(r)){
					String key = forceSameChr ? (r.getReadName() + "\t" + r.getReferenceName()) : r.getReadName();
					if(samsEnd2CTHash.containsKey(key)){
						ArrayList<SAMRecord> newCollection = samsEnd2CTHash.get(key);
						newCollection.add(r);
						samsEnd2CTHash.put(key, newCollection);
					}else{
						readsNameCollection.add(key);
						ArrayList<SAMRecord> newCollection = new ArrayList<SAMRecord>();
						newCollection.add(r);
						samsEnd2CTHash.put(key, newCollection);
					}
				}
			}
		}
		for(int i = 0; i < samsEnd2GA.length; i++){
			for(SAMRecord r : String2SamRecord(samsEnd2GA[i], samFileHeader, L2.get(i), true)){
				if(!failFlagFilter(r)){
					String key = forceSameChr ? (r.getReadName() + "\t" + r.getReferenceName()) : r.getReadName();
					if(samsEnd2GAHash.containsKey(key)){
						ArrayList<SAMRecord> newCollection = samsEnd2GAHash.get(key);
						newCollection.add(r);
						samsEnd2GAHash.put(key, newCollection);
					}else{
						readsNameCollection.add(key);
						ArrayList<SAMRecord> newCollection = new ArrayList<SAMRecord>();
						newCollection.add(r);
						samsEnd2GAHash.put(key, newCollection);
					}
				}
			}
		}
		
		
		
		
		//for each reads name, find the best pair to represent it.
		for(String name : readsNameCollection){
			String[] splitin = name.split("\\t");
			String readname = splitin[0];
			if(samsEnd1CTHash.containsKey(name) && samsEnd2CTHash.containsKey(name)){
				for(SAMRecord r1 : samsEnd1CTHash.get(name)){
					for(SAMRecord r2 : samsEnd2CTHash.get(name)){
						//if((r1.getReadNegativeStrandFlag() && !r2.getReadNegativeStrandFlag()) || (r2.getReadNegativeStrandFlag() && !r1.getReadNegativeStrandFlag())){
							if(bestPairs.containsKey(readname)){
								
								Pair<SAMRecord, SAMRecord> tmp1 = comparingSamRecordPbat(bestPairs.get(readname), new Pair<SAMRecord, SAMRecord>(r1, r2), regionsEnzyme);
								//boolean not = false;
								//if(bestPairs.get(readname).getFirst().getAlignmentStart() != tmp1.getFirst().getAlignmentStart() && bestPairs.get(readname).getSecond().getAlignmentStart() != tmp1.getSecond().getAlignmentStart()){
								//	System.err.println("original:" + bestPairs.get(readname).getFirst().getAlignmentStart() + "\t" + bestPairs.get(readname).getSecond().getAlignmentStart());
								//	System.err.println("best:" + tmp1.getFirst().getAlignmentStart() + "\t" + tmp1.getSecond().getAlignmentStart());
								//	not = true;
								//}
								
								
								bestPairs.put(readname, tmp1);
								//if(not){
								//	System.err.println("changed:" + bestPairs.get(readname).getFirst().getAlignmentStart() + "\t" + bestPairs.get(readname).getSecond().getAlignmentStart());
								//}
								
							}else{
								bestPairs.put(readname, new Pair<SAMRecord, SAMRecord>(r1, r2));
							}
						//}
					}
				}
				
				
				
			}
			if(samsEnd1GAHash.containsKey(name) && samsEnd2GAHash.containsKey(name)){
				for(SAMRecord r1 : samsEnd1GAHash.get(name)){
					for(SAMRecord r2 : samsEnd2GAHash.get(name)){
						//if((r1.getReadNegativeStrandFlag() && !r2.getReadNegativeStrandFlag()) || (r2.getReadNegativeStrandFlag() && !r1.getReadNegativeStrandFlag())){
							if(bestPairs.containsKey(readname)){
								bestPairs.put(readname, comparingSamRecordPbat(bestPairs.get(readname), new Pair<SAMRecord, SAMRecord>(r1, r2), regionsEnzyme));
							}else{
								bestPairs.put(readname, new Pair<SAMRecord, SAMRecord>(r1, r2));
							}
						//}
						
					}
				}
				
				
				
			}
			if(noDirectional){
				if(samsEnd1CTHash.containsKey(name) && samsEnd2GAHash.containsKey(name)){
					for(SAMRecord r1 : samsEnd1CTHash.get(name)){
						for(SAMRecord r2 : samsEnd2GAHash.get(name)){
							//if((r1.getReadNegativeStrandFlag() && !r2.getReadNegativeStrandFlag()) || (r2.getReadNegativeStrandFlag() && !r1.getReadNegativeStrandFlag())){
								if(bestPairs.containsKey(readname)){
									bestPairs.put(readname, comparingSamRecordPbat(bestPairs.get(readname), new Pair<SAMRecord, SAMRecord>(r1, r2), regionsEnzyme));
								}else{
									bestPairs.put(readname, new Pair<SAMRecord, SAMRecord>(r1, r2));
								}
							//}
							
						}
					}
				}
				
				if(samsEnd1GAHash.containsKey(name) && samsEnd2CTHash.containsKey(name)){
					for(SAMRecord r1 : samsEnd1GAHash.get(name)){
						for(SAMRecord r2 : samsEnd2CTHash.get(name)){
							//if((r1.getReadNegativeStrandFlag() && !r2.getReadNegativeStrandFlag()) || (r2.getReadNegativeStrandFlag() && !r1.getReadNegativeStrandFlag())){
								if(bestPairs.containsKey(readname)){
									bestPairs.put(readname, comparingSamRecordPbat(bestPairs.get(readname), new Pair<SAMRecord, SAMRecord>(r1, r2), regionsEnzyme));
								}else{
									bestPairs.put(readname, new Pair<SAMRecord, SAMRecord>(r1, r2));
								}
							//}
							
						}
					}
				}
			}
			
		}
		
		//foreach best pair, compute insertion size, fix the flag to be: paired and proper paired, it is 1st/2nd end, its mate is reverse or not..
		HashMap<String, Pair<SAMRecord, SAMRecord>> outputPairs = new HashMap<String, Pair<SAMRecord, SAMRecord>>();
		for(String name : bestPairs.keySet()){
			SAMRecord r1 = bestPairs.get(name).getFirst();
			SAMRecord r2 = bestPairs.get(name).getSecond();
			if(!r1.getContig().equalsIgnoreCase(r2.getContig())){
				if(!outputMateDiffChr){
					
					continue;
				}
			}
			
			//System.err.println(r1.getFlags() + "\t" + r1.getAlignmentStart() + "\t" + r1.getMateReferenceName() + "\t" + r1.getMateAlignmentStart() + "\t" + r1.getInferredInsertSize());
			//System.err.println(r2.getFlags() + "\t" + r2.getAlignmentStart() + "\t" + r2.getMateReferenceName() + "\t" + r2.getMateAlignmentStart() + "\t" + r2.getInferredInsertSize());
			SamPairUtil.setProperPairAndMateInfo(r1, r2, samFileHeader, ORITATION, true);
			//System.err.println(r1.getFlags() + "\t" + r1.getAlignmentStart() + "\t" + r1.getMateReferenceName() + "\t" + r1.getMateAlignmentStart() + "\t" + r1.getInferredInsertSize());
			//System.err.println(r2.getFlags() + "\t" + r2.getAlignmentStart() + "\t" + r2.getMateReferenceName() + "\t" + r2.getMateAlignmentStart() + "\t" + r2.getInferredInsertSize());
			outputPairs.put(name, new Pair<SAMRecord, SAMRecord>(r1, r2));
		}
		
		return outputPairs;
	}
	
	
	private ArrayList<SAMRecord> String2SamRecord(String lines, SAMFileHeader samFileHeader, String originalSeq, boolean SecondEnd){
		ArrayList<SAMRecord> records = new ArrayList<SAMRecord>();
		for(String line : lines.split("\\n")){
			SAMRecord record = new SAMRecord(samFileHeader);
			//System.err.println(line);
			String[] splitin = line.split("\\t");
			String readname = splitin[0];
			if(readNameSuffix != null){
				readname = readname.replaceAll(readNameSuffix, ""); // _\\d$ or _R\\d$ for _1 or _R1 in read suffix
			}
			record.setReadName(readname);
			record.setFlags(Integer.parseInt(splitin[1]));
			splitin[2]= splitin[2].replace("_CT_converted", "");
			splitin[2]= splitin[2].replace("_GA_converted", "");
			record.setReferenceName(splitin[2]);
			record.setAlignmentStart(Integer.parseInt(splitin[3]));
			record.setMappingQuality(Integer.parseInt(splitin[4]));
			record.setCigarString(splitin[5]);
			record.setMateReferenceName(splitin[6]);
			record.setMateAlignmentStart(0);
			record.setInferredInsertSize(0);
			String baseQ = splitin[10];
			if(pbat){
				record.setReadNegativeStrandFlag(!record.getReadNegativeStrandFlag());
				originalSeq = SequenceUtil.reverseComplement(originalSeq);
				baseQ = new StringBuilder(baseQ).reverse().toString();
			}
			//record.setReadString(splitin[9]);
			
			if(record.getReadNegativeStrandFlag()){
				originalSeq = SequenceUtil.reverseComplement(originalSeq);
				//originalSeq = new StringBuilder(originalSeq).reverse().toString(); //the string output by bwa need to be reversed complement when the mapping strand is negative...
				baseQ = new StringBuilder(baseQ).reverse().toString();
			}
			record.setReadString(modifySeqByCigar(originalSeq, record.getCigarString()));
			record.setBaseQualityString(baseQ);
			if(record.getReadBases().length != record.getBaseQualities().length){
				//throw new IllegalArgumentException("SAMrecord " + record.getReadName()  + " read bases length: " + record.getReadBases().length + " is different from baseQ length: " + record.getBaseQualities().length);
				//log.info("SAMrecord " + record.getReadName()  + " read bases length: " + record.getReadBases().length + " is different from baseQ length: " + record.getBaseQualities().length);

			}
			//System.err.println(record.getCigarLength() + "\t" + record.getCigarString());
			//for(char cigar : CigarUtil.cigarArrayFromString(splitin[5])){
			//	System.err.println(cigar);
			//}
			
			for(int i = 11; i < splitin.length; i++){
				if(splitin[i].isEmpty()){
					continue;
				}
				String[] splitTags = splitin[i].split(":");
				if(splitTags.length<3){
					continue;
				}
				
				if(splitTags[0].equalsIgnoreCase("SA") || splitTags[0].equalsIgnoreCase("XA")){
					continue;
				}
				//System.err.println(splitTags[0] + "\t" + splitTags[1] + "\t" + splitTags[2]);
				if(splitTags[1].equalsIgnoreCase("i")){
					record.setAttribute(splitTags[0], Integer.parseInt(splitTags[2]));
				}else{
					record.setAttribute(splitTags[0], splitTags[2]);
				}
				
			}
			//record.setAttribute(ORIGINAL_SEQ_TAG, originalSeq);
			record.setSupplementaryAlignmentFlag(false);
			record.setFirstOfPairFlag(!SecondEnd);
			record.setReadPairedFlag(true);
			record.setMateUnmappedFlag(false);
			record.setProperPairFlag(true);
			record.setSecondOfPairFlag(SecondEnd);
			record.setAttribute("RG", rgId);
			//record.setHeader(samFileHeader);
			records.add(record);
			//if(readname.equalsIgnoreCase("E00440:592:HCHF7CCX2:2:1101:10500:1327")){
			//	System.err.println(line + "\t" + record.getNotPrimaryAlignmentFlag() + "\t" + record.toString());
			//}
		}
		
		return records;
		
	}
	
	private boolean failFlagFilter(SAMRecord r){
		return r.getReadUnmappedFlag() || r.getNotPrimaryAlignmentFlag()
				|| r.getReadFailsVendorQualityCheckFlag() ;
	}
	
	private SAMRecord comparingSamRecord(SAMRecord r1, SAMRecord r2){
		if(r2.getMappingQuality() > r1.getMappingQuality()){
			return r2;
		}else if(r2.getMappingQuality() == r1.getMappingQuality()){
			if((Integer)r2.getAttribute(ALIGNMENT_SCORE) > (Integer)r1.getAttribute(ALIGNMENT_SCORE)){
				return r2;
			}else if((Integer)r2.getAttribute(ALIGNMENT_SCORE) == (Integer)r1.getAttribute(ALIGNMENT_SCORE)){
				if((Integer)r2.getAttribute(MISMATCHES) < (Integer)r1.getAttribute(MISMATCHES)){
					return r2;
				}else if((Integer)r2.getAttribute(MISMATCHES) == (Integer)r1.getAttribute(MISMATCHES)){
					int cigarMLenR1 = getCigarMLen(r1.getCigarString());
					int cigarMLenR2 = getCigarMLen(r2.getCigarString());
					
					if(cigarMLenR2 > cigarMLenR1){
						return r2;
					}
				}
			}
		}
		return r1;
	}
	
	private Pair<SAMRecord, SAMRecord> comparingSamRecord(Pair<SAMRecord, SAMRecord> r1, Pair<SAMRecord, SAMRecord> r2, HashMap<String,IntervalTree<String>> regionsEnzyme){
		if((r2.getFirst().getMappingQuality()>0 && r2.getSecond().getMappingQuality() > 0) && (r1.getFirst().getMappingQuality()==0 || r1.getSecond().getMappingQuality() == 0)){
			return r2;
		}else{
			if((r2.getFirst().getMappingQuality()+r2.getSecond().getMappingQuality()) > (r1.getFirst().getMappingQuality()+r1.getSecond().getMappingQuality())){
				return r2;
			}else if((r2.getFirst().getMappingQuality()+r2.getSecond().getMappingQuality()) == (r1.getFirst().getMappingQuality()+r1.getSecond().getMappingQuality())){
				if((!r1.getFirst().getContig().equalsIgnoreCase(r1.getSecond().getContig())) && r2.getFirst().getContig().equalsIgnoreCase(r2.getSecond().getContig())){
					return r2;
				}else{
					boolean returnR2 = false;
					if(!regionsEnzyme.isEmpty()){
						if(regionsEnzyme.containsKey(r2.getFirst().getContig()) && regionsEnzyme.containsKey(r2.getSecond().getContig())){
							if(!regionsEnzyme.containsKey(r1.getFirst().getContig()) || !regionsEnzyme.containsKey(r1.getSecond().getContig())){
								returnR2 = true;
							}else{
								if(regionsEnzyme.get(r2.getFirst().getContig()).overlappers(r2.getFirst().getAlignmentStart()-50, r2.getFirst().getAlignmentEnd()+50) != null 
										&& regionsEnzyme.get(r2.getSecond().getContig()).overlappers(r2.getSecond().getAlignmentStart()-50, r2.getSecond().getAlignmentEnd()+50) != null){
									if(regionsEnzyme.get(r1.getFirst().getContig()).overlappers(r1.getFirst().getAlignmentStart()-50, r1.getFirst().getAlignmentEnd()+50) == null 
											|| regionsEnzyme.get(r1.getSecond().getContig()).overlappers(r1.getSecond().getAlignmentStart()-50, r1.getSecond().getAlignmentEnd()+50) == null){
										returnR2 = true;
									}
								}
							}
						}
					}
					if(returnR2){
						return r2;
					}else{
						if(((Integer)r2.getFirst().getAttribute(ALIGNMENT_SCORE) + (Integer)r2.getSecond().getAttribute(ALIGNMENT_SCORE)) > ((Integer)r1.getFirst().getAttribute(ALIGNMENT_SCORE) + (Integer)r1.getSecond().getAttribute(ALIGNMENT_SCORE))){
							return r2;
						}else if(((Integer)r2.getFirst().getAttribute(ALIGNMENT_SCORE) + (Integer)r2.getSecond().getAttribute(ALIGNMENT_SCORE)) == ((Integer)r1.getFirst().getAttribute(ALIGNMENT_SCORE) + (Integer)r1.getSecond().getAttribute(ALIGNMENT_SCORE))){
							if(((Integer)r2.getFirst().getAttribute(MISMATCHES) + (Integer)r2.getSecond().getAttribute(MISMATCHES)) < ((Integer)r1.getFirst().getAttribute(MISMATCHES) + (Integer)r1.getSecond().getAttribute(MISMATCHES))){ // can not use mismatches for comparison here, it might be the reason for under-estimation of methylation level...
								return r2;
							}else if(((Integer)r2.getFirst().getAttribute(MISMATCHES) + (Integer)r2.getSecond().getAttribute(MISMATCHES)) == ((Integer)r1.getFirst().getAttribute(MISMATCHES) + (Integer)r1.getSecond().getAttribute(MISMATCHES))){
								int cigarMLenP1 = getCigarMLen(r1.getFirst().getCigarString()) + getCigarMLen(r1.getSecond().getCigarString());
								int cigarMLenP2 = getCigarMLen(r2.getFirst().getCigarString()) + getCigarMLen(r2.getSecond().getCigarString());
								if(cigarMLenP2 > cigarMLenP1){
									return r2;
								}
							}
						}
					}
				}
			}
		}
		
		return r1;
	}
	
	private Pair<SAMRecord, SAMRecord> comparingSamRecord(Pair<SAMRecord, SAMRecord> r1, Pair<SAMRecord, SAMRecord> r2){
		if(r2.getFirst().getMappingQuality() * r2.getSecond().getMappingQuality() != 0 && r1.getFirst().getMappingQuality() * r1.getSecond().getMappingQuality() == 0){ // exclude some reads that is biased, only one is 30, the other is 0...
			return r2;
		}else if(r2.getFirst().getMappingQuality() * r2.getSecond().getMappingQuality() == 0 && r1.getFirst().getMappingQuality() * r1.getSecond().getMappingQuality() == 0) {
			if ((r2.getFirst().getMappingQuality() + r2.getSecond().getMappingQuality()) > (r1.getFirst().getMappingQuality() + r1.getSecond().getMappingQuality())) {
				return r2;
			} else if ((r2.getFirst().getMappingQuality() + r2.getSecond().getMappingQuality()) == (r1.getFirst().getMappingQuality() + r1.getSecond().getMappingQuality())) {
				int cigarMLenP1 = getCigarMLen(r1.getFirst().getCigarString()) + getCigarMLen(r1.getSecond().getCigarString());
				int cigarMLenP2 = getCigarMLen(r2.getFirst().getCigarString()) + getCigarMLen(r2.getSecond().getCigarString());
				if (cigarMLenP2 > cigarMLenP1) {
					return r2;
				} else if (cigarMLenP2 == cigarMLenP1) {
					if (((Integer) r2.getFirst().getAttribute(ALIGNMENT_SCORE) + (Integer) r2.getSecond().getAttribute(ALIGNMENT_SCORE)) > ((Integer) r1.getFirst().getAttribute(ALIGNMENT_SCORE) + (Integer) r1.getSecond().getAttribute(ALIGNMENT_SCORE))) {
						return r2;
					} else if (((Integer) r2.getFirst().getAttribute(ALIGNMENT_SCORE) + (Integer) r2.getSecond().getAttribute(ALIGNMENT_SCORE)) == ((Integer) r1.getFirst().getAttribute(ALIGNMENT_SCORE) + (Integer) r1.getSecond().getAttribute(ALIGNMENT_SCORE))) {
						if (((Integer) r2.getFirst().getAttribute(MISMATCHES) + (Integer) r2.getSecond().getAttribute(MISMATCHES)) < ((Integer) r1.getFirst().getAttribute(MISMATCHES) + (Integer) r1.getSecond().getAttribute(MISMATCHES))) { // can not use mismatches for comparison here, it might be the reason for under-estimation of methylation level...
							return r2;
						}
					}
				}


			}

		}
		return r1;
	}
	
	private Pair<SAMRecord, SAMRecord> comparingSamRecordPbat(Pair<SAMRecord, SAMRecord> r1, Pair<SAMRecord, SAMRecord> r2, HashMap<String,IntervalTree<String>> regionsEnzyme){
		if((r2.getFirst().getMappingQuality()>0 && r2.getSecond().getMappingQuality() > 0) && (r1.getFirst().getMappingQuality()==0 || r1.getSecond().getMappingQuality() == 0)){
			return r2;
		}else{
			if((r2.getFirst().getMappingQuality()+r2.getSecond().getMappingQuality()) > (r1.getFirst().getMappingQuality()+r1.getSecond().getMappingQuality())){
				return r2;
			}else if((r2.getFirst().getMappingQuality()+r2.getSecond().getMappingQuality()) == (r1.getFirst().getMappingQuality()+r1.getSecond().getMappingQuality())){
				if((!r1.getFirst().getContig().equalsIgnoreCase(r1.getSecond().getContig())) && r2.getFirst().getContig().equalsIgnoreCase(r2.getSecond().getContig())){
					return r2;
				}else{
					boolean returnR2 = false;
					if(!regionsEnzyme.isEmpty()){
						if(regionsEnzyme.containsKey(r2.getFirst().getContig()) && regionsEnzyme.containsKey(r2.getSecond().getContig())){
							if(!regionsEnzyme.containsKey(r1.getFirst().getContig()) || !regionsEnzyme.containsKey(r1.getSecond().getContig())){
								returnR2 = true;
							}else{
								if(regionsEnzyme.get(r2.getFirst().getContig()).overlappers(r2.getFirst().getAlignmentStart()-50, r2.getFirst().getAlignmentEnd()+50) != null 
										&& regionsEnzyme.get(r2.getSecond().getContig()).overlappers(r2.getSecond().getAlignmentStart()-50, r2.getSecond().getAlignmentEnd()+50) != null){
									if(regionsEnzyme.get(r1.getFirst().getContig()).overlappers(r1.getFirst().getAlignmentStart()-50, r1.getFirst().getAlignmentEnd()+50) == null 
											|| regionsEnzyme.get(r1.getSecond().getContig()).overlappers(r1.getSecond().getAlignmentStart()-50, r1.getSecond().getAlignmentEnd()+50) == null){
										returnR2 = true;
									}
								}
							}
						}
					}
					if(returnR2){
						return r2;
					}else{
						if(((Integer)r2.getFirst().getAttribute(ALIGNMENT_SCORE) + (Integer)r2.getSecond().getAttribute(ALIGNMENT_SCORE)) > ((Integer)r1.getFirst().getAttribute(ALIGNMENT_SCORE) + (Integer)r1.getSecond().getAttribute(ALIGNMENT_SCORE))){
							return r2;
						}else if(((Integer)r2.getFirst().getAttribute(ALIGNMENT_SCORE) + (Integer)r2.getSecond().getAttribute(ALIGNMENT_SCORE)) == ((Integer)r1.getFirst().getAttribute(ALIGNMENT_SCORE) + (Integer)r1.getSecond().getAttribute(ALIGNMENT_SCORE))){
							if(((Integer)r2.getFirst().getAttribute(MISMATCHES) + (Integer)r2.getSecond().getAttribute(MISMATCHES)) < ((Integer)r1.getFirst().getAttribute(MISMATCHES) + (Integer)r1.getSecond().getAttribute(MISMATCHES))){ // can not use mismatches for comparison here, it might be the reason for under-estimation of methylation level...
								return r2;
							}else if(((Integer)r2.getFirst().getAttribute(MISMATCHES) + (Integer)r2.getSecond().getAttribute(MISMATCHES)) == ((Integer)r1.getFirst().getAttribute(MISMATCHES) + (Integer)r1.getSecond().getAttribute(MISMATCHES))){
								int cigarMLenP1 = getCigarMLen(r1.getFirst().getCigarString()) + getCigarMLen(r1.getSecond().getCigarString());
								int cigarMLenP2 = getCigarMLen(r2.getFirst().getCigarString()) + getCigarMLen(r2.getSecond().getCigarString());
								if(cigarMLenP2 > cigarMLenP1){
									return r2;
								}
							}
						}
					}
				}
			}
		}
		
		return r1;
	}
	
	private Pair<SAMRecord, SAMRecord> comparingSamRecordPbat(Pair<SAMRecord, SAMRecord> r1, Pair<SAMRecord, SAMRecord> r2){
		if((r2.getFirst().getMappingQuality()+r2.getSecond().getMappingQuality()) > (r1.getFirst().getMappingQuality()+r1.getSecond().getMappingQuality())){
			return r2;
		}else if((r2.getFirst().getMappingQuality()+r2.getSecond().getMappingQuality()) == (r1.getFirst().getMappingQuality()+r1.getSecond().getMappingQuality())){
			if(Math.abs(r1.getFirst().getAlignmentStart() - r1.getSecond().getAlignmentStart()) > Math.abs(r2.getFirst().getAlignmentStart() - r2.getSecond().getAlignmentStart())){
				return r2;
			}else if(Math.abs(r1.getFirst().getAlignmentStart() - r1.getSecond().getAlignmentStart()) == Math.abs(r2.getFirst().getAlignmentStart() - r2.getSecond().getAlignmentStart())){
				if(((Integer)r2.getFirst().getAttribute(ALIGNMENT_SCORE) + (Integer)r2.getSecond().getAttribute(ALIGNMENT_SCORE)) > ((Integer)r1.getFirst().getAttribute(ALIGNMENT_SCORE) + (Integer)r1.getSecond().getAttribute(ALIGNMENT_SCORE))){
					return r2;
				}else if(((Integer)r2.getFirst().getAttribute(ALIGNMENT_SCORE) + (Integer)r2.getSecond().getAttribute(ALIGNMENT_SCORE)) == ((Integer)r1.getFirst().getAttribute(ALIGNMENT_SCORE) + (Integer)r1.getSecond().getAttribute(ALIGNMENT_SCORE))){
					if(((Integer)r2.getFirst().getAttribute(MISMATCHES) + (Integer)r2.getSecond().getAttribute(MISMATCHES)) < ((Integer)r1.getFirst().getAttribute(MISMATCHES) + (Integer)r1.getSecond().getAttribute(MISMATCHES))){ // can not use mismatches for comparison here, it might be the reason for under-estimation of methylation level...
						return r2;
					}else if(((Integer)r2.getFirst().getAttribute(MISMATCHES) + (Integer)r2.getSecond().getAttribute(MISMATCHES)) == ((Integer)r1.getFirst().getAttribute(MISMATCHES) + (Integer)r1.getSecond().getAttribute(MISMATCHES))){
						int cigarMLenP1 = getCigarMLen(r1.getFirst().getCigarString()) + getCigarMLen(r1.getSecond().getCigarString());
						int cigarMLenP2 = getCigarMLen(r2.getFirst().getCigarString()) + getCigarMLen(r2.getSecond().getCigarString());
						if(cigarMLenP2 > cigarMLenP1){
							return r2;
						}
					}
				}
			}
			
			
			
			
		}
		return r1;
	}
	
	private int getCigarMLen(String cigarString){
		int cigarMLen = 0;
		for(char cigar : CigarUtil.cigarArrayFromString(cigarString)){
			if(cigar == 'M'){
				cigarMLen++;
			}
		}
		return cigarMLen;
	}
	
	private String modifySeqByCigar(String seq, String cigarString){
		char[] cigarList = CigarUtil.cigarArrayFromString(cigarString);
		byte[] seqByte = seq.getBytes();
		ArrayList<Byte> seqsNew = new ArrayList<Byte>();
		
		int offSet = 0;
		for(char cigar : cigarList){
			if(cigar == 'M' || cigar == 'S' || cigar == 'I'){
				seqsNew.add(seqByte[offSet]);
				offSet++;
			}else if(cigar == 'H'){
				offSet++;
			}
		}
		return new String(ArrayUtils.toPrimitive(seqsNew.toArray(new Byte[seqsNew.size()])));
	}
	
	private void initiate(){
		startTime = System.currentTimeMillis();
		System.loadLibrary("bwajni");
		ORITATION = new ArrayList<SamPairUtil.PairOrientation>();
		ORITATION.add(SamPairUtil.PairOrientation.FR);
		ORITATION.add(SamPairUtil.PairOrientation.RF);
		
	}

	private void finish(){
		long endTime   = System.currentTimeMillis();
		double totalTime = endTime - startTime;
		totalTime /= 1000;
		double totalTimeMins = totalTime/60;
		double totalTimeHours = totalTime/3600;
		
		log.info(totalReads + " reads in total");
		log.info("Unique Mapped and paired " + mappedReads + " reads (" + String.format("%.2f",100*(double)mappedReads/(double)totalReads) + "%)");
		log.info("Bhmem's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
	}
}
