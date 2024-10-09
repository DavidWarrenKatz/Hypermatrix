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
import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.CigarUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.random.MersenneTwister;
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
 * map long (>75bp) PE bisulfite converted reads into genome (HiC or WGBS)
 * will use Xr tag to represent the assignment of the genome to each read. 0 represent ambiguous, 1 for the first genome, 2 for the 2nd genome ...
 * Output: merge.referenceGenome.bam ref1.bam ref2.bam ... ambiguous.bam
 * Ambiguous reads in merge.bam file will be randomly assign to either genome and then LiftOver to the reference genome.
 * Due to the limitation of SAM format for multi genome mapping reads. It will output multiple bam file + 1 ambiguous reads bam file
 * For ambiguous bam file, the sam file header will just use the first reference genome's head.
 * Or use randomAmbiguous mode to randomly assign the ambiguous reads to one of the reference genome, so not output ambiguous reads bam file
 * make merg.bam file which use list of chain file between personal genome and reference genome...Then this make a new question, when new CpG only in one allele of the genome, it will cause Bis-SNP to call it 
 * as a SNP there... When a CpG is disrupted in one allele, it will also called as SNP there.. and it will be filtered out in the final output...
 * it will be better to output it when one allele is C/T SNP...
 * TODO: current AlleleSeq only looks at phased VCF and generate two .fa and chain file. Better to just use .vcf.gz file directly here, generate graph model and just do the mapping directly..
 * better to change the index method, build index for the reference genome, then extract these k-mer overlapped with SNP.vcf.gz add more indexer into it, but add tmp marker of these k-mer.
 */
public class BhmemMultRef {

	@Option(name="-refs",usage="reference genome fasta file, currently only support two references. -refs CAST.mm9.fa -refs 129.mm9.fa Default: null")
	public ArrayList<String> refs = null;
	@Option(name="-chains",usage="chain file generated by AlleleSeq or Pegasus for the conversion between two allele's genome and reference genome. Default: null")
	public ArrayList<String> chains = null;
	
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

	@Option(name="-noHicMode",usage="in normal WGBS/RRBS/NOMe-seq mode, not in BisulfiteHiC mode, not implemented yet. WGBS mode is moved to BsmemMulRef.java. default: false")
	public boolean noHicMode = false;

	@Option(name="-addAmbiguousReads",usage="output ambiguous reads into each genome's bam file. default: false")
	public boolean addAmbiguousReads = false;
	
	@Option(name="-noChimericMode",usage="not use chimeric mapped reads to increase the mapping efficiency, not implemented yet. default: false")
	public boolean noChimericMode = false;

	@Option(name="-trimReadName",usage="two ends read name are different wtih /1, /2 suffix, need to enable this option to trim it... default: false")
	public boolean trimReadName = false;
	
	@Option(name="-ratio",usage="minimum probability ratio to be assigned to the best genome, rather than the second best. Otherwise, reads will be assigned as ambiguous. default: 0")
	public double ratio = 0;
	
	@Option(name="-score",usage="minimum aligntment score to use. equal to -T in bwa mem. default: 0")
	public int score = 0;	
	
	@Option(name="-seed",usage="integer seed for the reproducibility. -1 means no seed set. default: -1")
	public int seed = -1;
	
	@Option(name="-t",usage="number of CPU threads to use. default: 1")
	public int thread = 1;
	
	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "BhmemMultRef [opts] ref.fa outputPrefix r1.fq[.gz] r2.fq[.gz]";
	
	private static Logger log = Logger.getLogger(BhmemMultRef.class);

	private static long startTime = -1;
	private static long totalReads = 0L;
	private static long mappedReads = 0L;
	private static long uniqMappedSameChrReads = 0L; //both end within the same chr
	private static long uniqMappedSameChrMapq1Reads = 0L; //both end MapQ >=1
	private static long uniqMappedSameChrMapq30Reads = 0L; //both end MapQ >=30
	private static long uniqMappedSameChrMapq1LongInsReads = 0L; //both end MapQ >=1 and insertion size > 20,000
	private static long uniqMappedSameChrMapq30LongInsReads = 0L; //both end MapQ >=30 and insertion size > 20,000
	
	//private final static String ORIGINAL_SEQ_TAG = "OS";
	private final static String MISMATCHES = "NM";
	private final static String ALIGNMENT_SCORE = "AS";
	private static List<SamPairUtil.PairOrientation> ORITATION;
	
	private MersenneTwister generator = null;
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		BhmemMultRef bhmem = new BhmemMultRef();
		BasicConfigurator.configure();
		bhmem.doMain(args);
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

					String refFile = arguments.get(0);
					String outputPrefix = arguments.get(1);
					String fastq1 = arguments.get(2);
					String fastq2 = arguments.get(3);

					if(refs == null || refs.size() < 2 || chains == null || chains.size() < 2 || refs.size() != chains.size()) throw new CmdLineException(parser, USAGE, new Throwable());
					
					initiate();	
					
					SAMFileHeader samFileHeaderMerge = new SAMFileHeader();
					File sequenceDictMerge = new File(refFile.replaceAll(".fa", ".dict"));
					if(sequenceDictMerge.exists() && !sequenceDictMerge.isDirectory()){
						samFileHeaderMerge.setSequenceDictionary(SAMSequenceDictionaryExtractor.extractDictionary(sequenceDictMerge));
					}else{
						CreateSequenceDictionary csd = new CreateSequenceDictionary();
						samFileHeaderMerge.setSequenceDictionary(csd.makeSequenceDictionary(sequenceDictMerge));
					}
					if(rgId != null && !rgId.isEmpty() && rgSm != null && !rgSm.isEmpty()){
						SAMReadGroupRecord rg = new SAMReadGroupRecord(rgId);
						rg.setLibrary(rgId);
						rg.setSample(rgSm);
						rg.setPlatform("illumina");
						rg.setPlatformUnit("barcode");
						ArrayList<SAMReadGroupRecord> rgList = new ArrayList<SAMReadGroupRecord>();
						rgList.add(rg);
						samFileHeaderMerge.setReadGroups(rgList);
					}
					
					//merge.bam file writer
					SAMFileWriter writerMerge = new SAMFileWriterFactory().makeBAMWriter(samFileHeaderMerge, true, new File(outputPrefix + ".merged.bam"));
					
					
					//ambiguous bam file
					SAMFileWriter writerAmb = new SAMFileWriterFactory().makeBAMWriter(samFileHeaderMerge, true, new File(outputPrefix + ".ambiguous.bam"));
					
					
					ArrayList<BwaMem> refCTList = new ArrayList<BwaMem>();
					ArrayList<BwaMem> refGAList = new ArrayList<BwaMem>();
					ArrayList<BwaIndex> refCTIndexList = new ArrayList<BwaIndex>();
					ArrayList<BwaIndex> refGAIndexList = new ArrayList<BwaIndex>();
					
					ArrayList<LiftOver> chainList = new ArrayList<LiftOver>();
					
					ArrayList<SAMFileHeader> samFileHeaderList =  new ArrayList<SAMFileHeader>();
					ArrayList<SAMFileWriter> samFileWriterList =  new ArrayList<SAMFileWriter>();
					
					
					
					int numRef = 1;
					for(String ref : refs){
						SAMFileHeader samFileHeader = new SAMFileHeader();
						File sequenceDict = new File(ref.replaceAll(".fa", ".dict"));
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
						
						
						samFileHeaderList.add(samFileHeader);
						SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
						SAMFileWriter writer = writerFactory.makeBAMWriter(samFileHeader, true, new File(outputPrefix + ".alignToRef_" + numRef + ".bam"));
						samFileWriterList.add(writer);
						
						LiftOver chain = new LiftOver(new File(chains.get(numRef-1)));
						chain.setLiftOverMinMatch(0.9);
						chainList.add(chain);
						
						BwaIndex indexCT=new BwaIndex(new File(new File(ref).getParent() + "/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa"));
						BwaMem memCT=new BwaMem(indexCT);
						memCT.updateMoreScoringParameters(thread, score, 1, clip5, clip3); //1 means use softclip supplement mapping. 0 means NO softclip supplement mapping.
						refCTList.add(memCT);
						refCTIndexList.add(indexCT);
						BwaIndex indexGA=new BwaIndex(new File(new File(ref).getParent() + "/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"));
						BwaMem memGA=new BwaMem(indexGA);
						memGA.updateMoreScoringParameters(thread, score, 1, clip5, clip3); //1 means use softclip supplement mapping. 0 means NO softclip supplement mapping.
						refGAList.add(memGA);
						refGAIndexList.add(indexGA);
						numRef++;
					}
					
					
					

					
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
								if(addAmbiguousReads){
									mapReadsToMultiGenomeAddAmbiguous(L1CT, L2GA, L1, L2, refCTList, refGAList, samFileHeaderList, samFileWriterList);
								}else{
									mapReadsToMultiGenome(L1CT, L2GA, L1, L2, refCTList, refGAList, samFileHeaderList, samFileWriterList, chainList, samFileHeaderMerge, writerMerge, writerAmb);
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
						String r1Name = read1.getName();
						String r2Name = read2.getName();
						if(trimReadName){
							Pattern replace = Pattern.compile("/\\d+$");
			                Matcher matcher1 = replace.matcher(r1Name);
			                r1Name=matcher1.replaceAll("");
			                Matcher matcher2 = replace.matcher(r2Name);
			                r2Name=matcher2.replaceAll("");
						}
						ShortRead read1ct=new ShortRead(r1Name, seq1.replace('C', 'T').getBytes() ,read1.getQualities());
						ShortRead read2ga=new ShortRead(r2Name, seq2.replace('G', 'A').getBytes() ,read2.getQualities());
						
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
					for(int i = 0; i < refs.size(); i++){
						refCTIndexList.get(i).close();
						refCTList.get(i).dispose();
						refGAIndexList.get(i).close();
						refGAList.get(i).dispose();
						
					}
					writerMerge.close();
					writerAmb.close();
					for(SAMFileWriter writer : samFileWriterList){
						writer.close();
					}
					
					finish();
	}
	
	private void mapReadsToMultiGenomeAddAmbiguous(List<ShortRead> L1CT, List<ShortRead> L2GA, 
			List<String> L1, List<String> L2, ArrayList<BwaMem> refCTList, ArrayList<BwaMem> refGAList, ArrayList<SAMFileHeader> samFileHeaderList, ArrayList<SAMFileWriter> samFileWriterList) throws IOException{
		
		//go over each genome for the mapping, and select the best genome for the reads
		for(int i = 0; i < refs.size(); i++){
			String[] samsEnd1CT = refCTList.get(i).align(L1CT);
			String[] samsEnd1GA = refGAList.get(i).align(L1CT);
			String[] samsEnd2CT = refCTList.get(i).align(L2GA);
			String[] samsEnd2GA = refGAList.get(i).align(L2GA);
			HashMap<String, Pair<SAMRecord, SAMRecord>> readPairs = joinTwoEnds(samsEnd1CT, samsEnd1GA, samsEnd2CT, samsEnd2GA, samFileHeaderList.get(i), L1, L2);
			for(String key : readPairs.keySet()){
				Pair<SAMRecord, SAMRecord> tmpPair = readPairs.get(key);
				SAMRecord r1 = tmpPair.getFirst();
				SAMRecord r2 = tmpPair.getSecond();
				
				samFileWriterList.get(i).addAlignment(r1);
				samFileWriterList.get(i).addAlignment(r2);
				if(i==0){
					mappedReads+=2;
					
					if(r1.getReferenceIndex() == r2.getReferenceIndex()){
						uniqMappedSameChrReads += 2;
						if(Math.min(r1.getMappingQuality(), r2.getMappingQuality()) >= 1){
							uniqMappedSameChrMapq1Reads += 2;
							if(Math.abs(r1.getInferredInsertSize())>20000){
								uniqMappedSameChrMapq1LongInsReads += 2;
							}
							if(Math.min(r1.getMappingQuality(), r2.getMappingQuality()) >= 30){
								uniqMappedSameChrMapq30Reads += 2;
								if(Math.abs(r1.getInferredInsertSize())>20000){
									uniqMappedSameChrMapq30LongInsReads += 2;
								}
							}
						}
					}
				}
				
			}
			
		}
		
		
	}
	
	private void mapReadsToMultiGenome(List<ShortRead> L1CT, List<ShortRead> L2GA, 
			List<String> L1, List<String> L2, ArrayList<BwaMem> refCTList, ArrayList<BwaMem> refGAList, ArrayList<SAMFileHeader> samFileHeaderList, ArrayList<SAMFileWriter> samFileWriterList,
			ArrayList<LiftOver> chainList, SAMFileHeader samFileHeaderMerge, SAMFileWriter writerMerge, SAMFileWriter writerAmb) throws IOException{
		
		//go over each genome for the mapping, and select the best genome for the reads
		HashMap<String, Pair<Integer, Pair<SAMRecord, SAMRecord>>> bestPair = new HashMap<String, Pair<Integer, Pair<SAMRecord, SAMRecord>>>();
		HashMap<String, Pair<Integer, Pair<SAMRecord, SAMRecord>>> secondBestPair = new HashMap<String, Pair<Integer, Pair<SAMRecord, SAMRecord>>>();
		for(int i = 0; i < refs.size(); i++){
			String[] samsEnd1CT = refCTList.get(i).align(L1CT);
			String[] samsEnd1GA = refGAList.get(i).align(L1CT);
			String[] samsEnd2CT = refCTList.get(i).align(L2GA);
			String[] samsEnd2GA = refGAList.get(i).align(L2GA);
			HashMap<String, Pair<SAMRecord, SAMRecord>> readPairs = joinTwoEnds(samsEnd1CT, samsEnd1GA, samsEnd2CT, samsEnd2GA, samFileHeaderList.get(i), L1, L2);
			for(String key : readPairs.keySet()){
				Pair<SAMRecord, SAMRecord> tmpPair = readPairs.get(key);
				if(bestPair.containsKey(key)){
					Pair<Integer, Pair<SAMRecord, SAMRecord>> p = bestPair.get(key);
					int compare = getBestSamRecordPair(tmpPair, p.getSecond());
					if(compare == -1){
						bestPair.put(key, new Pair<Integer, Pair<SAMRecord, SAMRecord>>(i, tmpPair));
						secondBestPair.put(key, p);
					}else if(compare == 0){
						if(generator.nextDouble() < 0.5){
							bestPair.put(key, new Pair<Integer, Pair<SAMRecord, SAMRecord>>(i, tmpPair));
							secondBestPair.put(key, p);
						}
					}else{
						if(secondBestPair.containsKey(key)){
							Pair<Integer, Pair<SAMRecord, SAMRecord>> second = secondBestPair.get(key);
							int compare2nd = getBestSamRecordPair(tmpPair, second.getSecond());
							if(compare2nd == -1){
								secondBestPair.put(key, new Pair<Integer, Pair<SAMRecord, SAMRecord>>(i, tmpPair));
							}else if(compare2nd == 0){
								if(generator.nextDouble() < 0.5){
									secondBestPair.put(key, new Pair<Integer, Pair<SAMRecord, SAMRecord>>(i, tmpPair));
								}
							}
						}else{
							secondBestPair.put(key, new Pair<Integer, Pair<SAMRecord, SAMRecord>>(i, tmpPair));
						}
						
					}
				}else if(secondBestPair.containsKey(key)){
					Pair<Integer, Pair<SAMRecord, SAMRecord>> second = secondBestPair.get(key);
					int compare2nd = getBestSamRecordPair(tmpPair, second.getSecond());
					if(compare2nd == -1){
						bestPair.put(key, new Pair<Integer, Pair<SAMRecord, SAMRecord>>(i, tmpPair));
					}else if(compare2nd == 0){
						if(generator.nextDouble() < 0.5){
							secondBestPair.put(key, new Pair<Integer, Pair<SAMRecord, SAMRecord>>(i, tmpPair));
						}
					}
				}else{
					bestPair.put(key, new Pair<Integer, Pair<SAMRecord, SAMRecord>>(i, tmpPair));
				}
			}
			
		}
		
		for(String key : bestPair.keySet()){
			Pair<Integer, Pair<SAMRecord, SAMRecord>> p = bestPair.get(key);
			if(secondBestPair.containsKey(key)){
				Pair<Integer, Pair<SAMRecord, SAMRecord>> second = secondBestPair.get(key);
				if(mapqRatio(p.getSecond(), second.getSecond()) > ratio){ //now it will output it to the specific genome's bam file, otherwise, just output the best to the ambiguous bam file
					samFileWriterList.get(p.getFirst()).addAlignment(p.getSecond().getFirst());
					samFileWriterList.get(p.getFirst()).addAlignment(p.getSecond().getSecond());
				}
				
			}else{ // no second best, so specifically to this genome
				samFileWriterList.get(p.getFirst()).addAlignment(p.getSecond().getFirst());
				samFileWriterList.get(p.getFirst()).addAlignment(p.getSecond().getSecond());
			}
			
			//now output it to the merge bam file and ambiguous bam file
			//liftOver the alignment start position
			SAMRecord r1 = p.getSecond().getFirst();
			SAMRecord r2 = p.getSecond().getSecond();
			
			LiftOver chain = chainList.get(p.getFirst());
			String chr1 = r1.getReferenceName();
			int end1Start = r1.getAlignmentStart();
			String chr2 = r2.getReferenceName();
			int end2Start = r2.getAlignmentStart();
			// re-infer the insertion sizes
			Interval end1LiftOver = chain.liftOver(new Interval(chr1, end1Start, end1Start));
			Interval end2LiftOver = chain.liftOver(new Interval(chr2, end2Start, end2Start));
			if(end1LiftOver == null || end2LiftOver == null){
				continue;
			}
			String chr1LiftOver = end1LiftOver.getContig();
			int end1StartLiftOver = end1LiftOver.getStart();
			String chr2LiftOver = end2LiftOver.getContig();
			int end2StartLiftOver = end2LiftOver.getStart();
			r1.setHeader(samFileHeaderMerge);
			r2.setHeader(samFileHeaderMerge);
			if(!chr1.equals(chr1LiftOver)){
				r1.setReferenceName(chr1);
			}
			if(!chr2.equals(chr2LiftOver)){
				r2.setReferenceName(chr2);
			}
			if(end1StartLiftOver != end1Start){
				r1.setAlignmentStart(end1StartLiftOver);
			}
			if(end2StartLiftOver != end2Start){
				r2.setAlignmentStart(end2StartLiftOver);
			}
			SamPairUtil.setProperPairAndMateInfo(r1, r2, samFileHeaderMerge, ORITATION, true);
			//recalculate CIGAR string?TODO: it will be super slow if reload reference genome here... need to find out how bwa access genome sequences by index...
			writerMerge.addAlignment(r1);
			writerMerge.addAlignment(r2);
			
			writerAmb.addAlignment(r1);
			writerAmb.addAlignment(r2);
			
			mappedReads+=2;
			
			if(r1.getReferenceIndex() == r2.getReferenceIndex()){
				uniqMappedSameChrReads += 2;
				if(Math.min(r1.getMappingQuality(), r2.getMappingQuality()) >= 1){
					uniqMappedSameChrMapq1Reads += 2;
					if(Math.abs(r1.getInferredInsertSize())>20000){
						uniqMappedSameChrMapq1LongInsReads += 2;
					}
					if(Math.min(r1.getMappingQuality(), r2.getMappingQuality()) >= 30){
						uniqMappedSameChrMapq30Reads += 2;
						if(Math.abs(r1.getInferredInsertSize())>20000){
							uniqMappedSameChrMapq30LongInsReads += 2;
						}
					}
				}
			}
		}
		
	}
	
	private double mapqRatio(Pair<SAMRecord, SAMRecord> r1, Pair<SAMRecord, SAMRecord> r2){
		//double q1 = (1-Math.pow(10, 0-(r1.getFirst().getMappingQuality()/10))*(1-Math.pow(10, 0-(r1.getSecond().getMappingQuality())/10)));
		//double q2 = (1-Math.pow(10, 0-(r2.getFirst().getMappingQuality()/10))*(1-Math.pow(10, 0-(r2.getSecond().getMappingQuality())/10)));
		//return q1/q2;
		return r1.getFirst().getMappingQuality() + r1.getSecond().getMappingQuality() - r2.getFirst().getMappingQuality() - r2.getSecond().getMappingQuality();
		
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
					String key = r.getReadName() + "\t" + r.getReferenceName();
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
					String key = r.getReadName() + "\t" + r.getReferenceName();
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
					String key = r.getReadName() + "\t" + r.getReferenceName();
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
					String key = r.getReadName() + "\t" + r.getReferenceName();
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
	
	
	private ArrayList<SAMRecord> String2SamRecord(String lines, SAMFileHeader samFileHeader, String originalSeq, boolean SecondEnd){
		ArrayList<SAMRecord> records = new ArrayList<SAMRecord>();
		for(String line : lines.split("\\n")){
			SAMRecord record = new SAMRecord(samFileHeader);
			//System.err.println(line);
			String[] splitin = line.split("\\t");
			record.setReadName(splitin[0]);
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
			//record.setReadString(splitin[9]);
			String baseQ = splitin[10];
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
			int cigarMLenR1 = getCigarMLen(r1.getCigarString());
			int cigarMLenR2 = getCigarMLen(r2.getCigarString());
			
			if(cigarMLenR2 > cigarMLenR1){
				return r2;
			}else if(cigarMLenR2 == cigarMLenR1){
				if((Integer)r2.getAttribute(MISMATCHES) < (Integer)r1.getAttribute(MISMATCHES)){
					return r2;
				}else if((Integer)r2.getAttribute(MISMATCHES) == (Integer)r1.getAttribute(MISMATCHES)){
					if((Integer)r2.getAttribute(ALIGNMENT_SCORE) > (Integer)r1.getAttribute(ALIGNMENT_SCORE)){
						return r2;
					}
				}
			}
		}
		return r1;
	}
	
	private Pair<SAMRecord, SAMRecord> comparingSamRecord(Pair<SAMRecord, SAMRecord> r1, Pair<SAMRecord, SAMRecord> r2){
		if(Math.min(r2.getFirst().getMappingQuality(), r2.getSecond().getMappingQuality()) > Math.min(r1.getFirst().getMappingQuality(), r1.getSecond().getMappingQuality())){
			return r2;
		}else if(Math.min(r2.getFirst().getMappingQuality(), r2.getSecond().getMappingQuality()) == Math.min(r1.getFirst().getMappingQuality(), r1.getSecond().getMappingQuality())){
			if((r2.getFirst().getMappingQuality()+r2.getSecond().getMappingQuality()) > (r1.getFirst().getMappingQuality()+r1.getSecond().getMappingQuality())){
				return r2;
			}else if((r2.getFirst().getMappingQuality()+r2.getSecond().getMappingQuality()) == (r1.getFirst().getMappingQuality()+r1.getSecond().getMappingQuality())){
				int cigarMLenP1 = getCigarMLen(r1.getFirst().getCigarString()) + getCigarMLen(r1.getSecond().getCigarString());
				int cigarMLenP2 = getCigarMLen(r2.getFirst().getCigarString()) + getCigarMLen(r2.getSecond().getCigarString());
				if(cigarMLenP2 > cigarMLenP1){
					return r2;
				}else if(cigarMLenP2 == cigarMLenP1){
					if(((Integer)r2.getFirst().getAttribute(MISMATCHES) + (Integer)r2.getSecond().getAttribute(MISMATCHES)) < ((Integer)r1.getFirst().getAttribute(MISMATCHES) + (Integer)r1.getSecond().getAttribute(MISMATCHES))){
						return r2;
					}else if(((Integer)r2.getFirst().getAttribute(MISMATCHES) + (Integer)r2.getSecond().getAttribute(MISMATCHES)) == ((Integer)r1.getFirst().getAttribute(MISMATCHES) + (Integer)r1.getSecond().getAttribute(MISMATCHES))){
						if(((Integer)r2.getFirst().getAttribute(ALIGNMENT_SCORE) + (Integer)r2.getSecond().getAttribute(ALIGNMENT_SCORE)) > ((Integer)r1.getFirst().getAttribute(ALIGNMENT_SCORE) + (Integer)r1.getSecond().getAttribute(ALIGNMENT_SCORE))){
							return r2;
						}
					}
				}
			}
		}
		
		return r1;
	}
	
	
	//1 means return r2 better, -1 means r1 better, 0 means equal
	private int getBestSamRecordPair(Pair<SAMRecord, SAMRecord> r1, Pair<SAMRecord, SAMRecord> r2){
		if(Math.min(r2.getFirst().getMappingQuality(), r2.getSecond().getMappingQuality()) > Math.min(r1.getFirst().getMappingQuality(), r1.getSecond().getMappingQuality())){
			return 1;
		}else if(Math.min(r2.getFirst().getMappingQuality(), r2.getSecond().getMappingQuality()) == Math.min(r1.getFirst().getMappingQuality(), r1.getSecond().getMappingQuality())){
			if((r2.getFirst().getMappingQuality()+r2.getSecond().getMappingQuality()) > (r1.getFirst().getMappingQuality()+r1.getSecond().getMappingQuality())){
				return 1;
			}else if((r2.getFirst().getMappingQuality()+r2.getSecond().getMappingQuality()) == (r1.getFirst().getMappingQuality()+r1.getSecond().getMappingQuality())){
				int cigarMLenP1 = getCigarMLen(r1.getFirst().getCigarString()) + getCigarMLen(r1.getSecond().getCigarString());
				int cigarMLenP2 = getCigarMLen(r2.getFirst().getCigarString()) + getCigarMLen(r2.getSecond().getCigarString());
				if(cigarMLenP2 > cigarMLenP1){
					return 1;
				}else if(cigarMLenP2 == cigarMLenP1){
					if(((Integer)r2.getFirst().getAttribute(MISMATCHES) + (Integer)r2.getSecond().getAttribute(MISMATCHES)) < ((Integer)r1.getFirst().getAttribute(MISMATCHES) + (Integer)r1.getSecond().getAttribute(MISMATCHES))){
						return 1;
					}else if(((Integer)r2.getFirst().getAttribute(MISMATCHES) + (Integer)r2.getSecond().getAttribute(MISMATCHES)) == ((Integer)r1.getFirst().getAttribute(MISMATCHES) + (Integer)r1.getSecond().getAttribute(MISMATCHES))){
						if(((Integer)r2.getFirst().getAttribute(ALIGNMENT_SCORE) + (Integer)r2.getSecond().getAttribute(ALIGNMENT_SCORE)) > ((Integer)r1.getFirst().getAttribute(ALIGNMENT_SCORE) + (Integer)r1.getSecond().getAttribute(ALIGNMENT_SCORE))){
							return 1;
						}else if(((Integer)r2.getFirst().getAttribute(ALIGNMENT_SCORE) + (Integer)r2.getSecond().getAttribute(ALIGNMENT_SCORE)) == ((Integer)r1.getFirst().getAttribute(ALIGNMENT_SCORE) + (Integer)r1.getSecond().getAttribute(ALIGNMENT_SCORE))){
							return 1;
						}else{
							return 0;
						}
					}
				}
			}
		}
		
		return -1;
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
		if(refs == null || refs.size() <=1 || refs.size() > 2){
			throw new IllegalArgumentException("Number of reference genome is not equal to two!!");
		}
		if(seed != -1){
			generator = new MersenneTwister(seed);
		}else{
			generator = new MersenneTwister();
		}
	}

	private void finish(){
		long endTime   = System.currentTimeMillis();
		double totalTime = endTime - startTime;
		totalTime /= 1000;
		double totalTimeMins = totalTime/60;
		double totalTimeHours = totalTime/3600;
		
		log.info(totalReads + " reads in total");
		log.info("Unique Mapped and paired " + mappedReads + " reads (" + String.format("%.2f",100*(double)mappedReads/(double)totalReads) + "%)");
		log.info("Unique Mapped, paired in the same chr " + uniqMappedSameChrReads + " reads (" + String.format("%.2f",100*(double)uniqMappedSameChrReads/(double)totalReads) + "%)");
		log.info("Unique Mapped, paired in the same chr  and MapQ >=1 " + uniqMappedSameChrMapq1Reads + " reads (" + String.format("%.2f",100*(double)uniqMappedSameChrMapq1Reads/(double)totalReads) + "%)");
		log.info("Unique Mapped, paired in the same chr  and MapQ >=30 " + uniqMappedSameChrMapq30Reads + " reads (" + String.format("%.2f",100*(double)uniqMappedSameChrMapq30Reads/(double)totalReads) + "%)");
		log.info("Unique Mapped, paired in the same chr, MapQ >=1 and insertionSize > 20000 " + uniqMappedSameChrMapq1LongInsReads + " reads (" + String.format("%.2f",100*(double)uniqMappedSameChrMapq1LongInsReads/(double)totalReads) + "%)");
		log.info("Unique Mapped, paired in the same chr, MapQ >=30 and insertionSize > 20000 " + uniqMappedSameChrMapq30LongInsReads + " reads (" + String.format("%.2f",100*(double)uniqMappedSameChrMapq30LongInsReads/(double)totalReads) + "%)");
		log.info("BhmemMultFa's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
	}
}

