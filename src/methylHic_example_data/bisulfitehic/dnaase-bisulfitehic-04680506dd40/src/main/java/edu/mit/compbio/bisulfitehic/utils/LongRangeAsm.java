/**
 * LongRangeAsm.java
 * Feb 27, 2017
 * 12:00:36 PM
 * yaping    lyping1986@gmail.com
 */
package main.java.edu.mit.compbio.bisulfitehic.utils;

import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.exception.CodecLineParsingException;
import htsjdk.tribble.exception.UnsortedFileException;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.TreeMap;
import java.util.zip.GZIPOutputStream;

import main.java.edu.mit.compbio.ccinference.utils.CcInferenceUtils;

import org.apache.commons.math3.util.Pair;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.GenomeLocSortedSet;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;


/**
 * output 
 */
public class LongRangeAsm {

	/**
	 *
	 * output.ASM.tab format:
	 * "#chr\tstart\tend\t"
	 * 				+ "AlleleA\tAlleleA_count\tAlleleA_C\tAlleleA_T\tAlleleA_meth_mean\t"
	 * 				+ "AlleleB\tAlleleB_count\tAlleleB_C\tAlleleB_T\tAlleleB_meth_mean\t"
	 * 				+ "End2_anchor_point\tEnd2_AlleleA_C\tEnd2_AlleleA_T\tEnd2_AlleleA_meth_mean\t"
	 * 				+ "End2_AlleleB_C\tEnd2_AlleleB_T\tEnd2_AlleleB_meth_mean\tMaxEnd2Dist\t"
	 * 				//+ "AlleleRef\tAlleleRef_count\tAlleleAlt\tAlleleAlt_count\t"
	 * 				+ "AlleleA_local_count\tAlleleA_local_C\tAlleleA_local_T\tAlleleA_local_meth_mean\t"
	 * 				+ "AlleleB_local_count\tAlleleB_local_C\tAlleleB_local_T\tAlleleB_local_meth_mean\t"
	 * 				+ "AlleleA_2nd_local_count\tAlleleA_2nd_local_C\tAlleleA_2nd_local_T\tAlleleA_2nd_local_meth_mean\t"
	 * 				+ "AlleleB_2nd_local_count\tAlleleB_2nd_local_C\tAlleleB_2nd_local_T\tAlleleB_2nd_local_meth_mean\t"
	 * 				+ "AlleleA_2nd_C\tAlleleA_2nd_T\tAlleleA_2nd_meth_mean\t"
	 * 				+ "AlleleB_2nd_C\tAlleleB_2nd_T\tAlleleB_2nd_meth_mean\t"
	 * 				+ "End2_AlleleA_2nd_C\tEnd2_AlleleA_2nd_T\tEnd2_AlleleA_2nd_meth_mean\t"
	 * 				+ "End2_AlleleB_2nd_C\tEnd2_AlleleB_2nd_T\tEnd2_AlleleB_2nd_meth_mean\n"
	 *
	 */
	@Option(name="-region",usage="specify the genomic region to look at, use interval format. e.g. chr1:1000-2000 or chr1 for the whole chromsome. default: null, the whole bam file")
	public String region = null;

	@Option(name="-turnOffBisulfiteFilter",usage="turn off bisulfite conversion reads filter. Default: false")
	public boolean turnOffBisulfiteFilter = false;

	@Option(name="-useBadMate",usage="use not proper paired reads. Default: false")
	public boolean useBadMate = false;

	@Option(name="-diffChr",usage="don't skip the region from two different chromosomes, default: not enabled")
	public boolean diffChr = false;

	@Option(name="-useUnpaired",usage="reads are mapped in single end mode, have not joined into pairs yet. Default: false")
	public boolean useUnpaired = false;

	@Option(name="-methyPatternSearch",usage="methylation pattern to search, allow IUPAC match code, such as HCG or GCH for NOMe-seq, default: CG")
	public String methyPatternSearch = "CG";

	@Option(name="-methyPatternCytPos",usage="cytosine position in -methyPatternSearch. for example, 1 for CG and 2 for HCG, default: 1")
	public int methyPatternCytPos = 1;

	@Option(name="-methyPattern2ndSearch",usage="the 2nd methylation pattern in the same reads to search, allow IUPAC match code, such as HCG or GCH for NOMe-seq, default: CG")
	public String methyPattern2ndSearch = null;

	@Option(name="-methyPattern2ndCytPos",usage="cytosine position in -methyPattern2ndSearch. for example, 1 for CG and 2 for HCG, default: 1")
	public int methyPattern2ndCytPos = 1;

	@Option(name="-bsConvPatternSearch",usage="bisulfite conversion pattern to search, allow IUPAC match code, such as WCH for NOMe-seq, default: CH")
	public String bsConvPatternSearch = "CH";

	@Option(name="-bsConvPatternCytPos",usage="cytosine position in -bsConvPatternSearch. for example, 1 for CH and 2 for WCH, default: 1")
	public int bsConvPatternCytPos = 1;

	@Option(name="-minMapQ",usage="minmum mapping quality score, default: 30")
	public int minMapQ = 30;

	@Option(name="-minBaseQ",usage="minmum mapping quality score, default: 20")
	public int minBaseQ = 20;

	@Option(name="-minDist",usage="minmum distance from two ends, default: 20000")
	public int minDist = 20000;

	@Option(name="-maxDist",usage="maximum distance from two ends, default: 1000000000")
	public int maxDist = 1000000000;
	
	@Option(name="-coverageRef",usage="specify the minimum  coverage in SNP's allele Ref, default: 5")
	public int coverageRef = 5;

	@Option(name="-coverageAlt",usage="specify the minimum coverage in SNP's allele Alt, default: 5")
	public int coverageAlt = 5;
	
	@Option(name="-qual",usage="specify the minimum genotyping quality of SNP call, default: 30")
	public int qual = 30;


	@Option(name="-partialVcf",usage="only use REF/ALT part of VCF rather than using genotype field in VCF file, default: false")
	public boolean partialVcf = false;

	@Option(name="-sample",usage="specify the sample number in the VCF file if contain multiple sample genotype. number started at 0. default: 0")
	public int sample = 0;


	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "LongRangeAsm [opts] genotype.vcf.gz input.bam output.ASM.tab.gz output.ref.tab.gz output.alt.tab.gz";
	
	private static Logger log = Logger.getLogger(LongRangeAsm.class);

	private static long startTime = -1;
	
	private OutputStreamWriter writer = null; 
	private OutputStreamWriter refWriter = null; 
	private OutputStreamWriter altWriter = null; 
	
	private final static int PURGE_INTERVAL = 1000;
	
	private long VCFCOUNT = 0L;
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		LongRangeAsm lac = new LongRangeAsm();
		BasicConfigurator.configure();
	    lac.doMain(args);

	}
	
	public void doMain(String[] args)
			throws Exception {

					CmdLineParser parser = new CmdLineParser(this);
					//parser.setUsageWidth(80);
					try
					{
						if(help || args.length < 5) throw new CmdLineException(USAGE);
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

					String vcfFile = arguments.get(0);
					String bamFile = arguments.get(1);
					String outputFile = arguments.get(2);
					String refFile = arguments.get(3);
					String altFile = arguments.get(4);
					
					VCFFileReader vcfReader = initiate(vcfFile, outputFile, refFile, altFile);
					 
					
					
					GenotypeList gl = readSnp(vcfFile, vcfReader);
					//System.err.println(gl.genomeLocsAlleleMap.size());
					processBam(bamFile, gl);
					
					finish();
	}
	

	
	private GenotypeList readSnp(String vcfFile, VCFFileReader vcfReader){
		CloseableIterator<VariantContext> itSnp = null;
		SAMSequenceDictionary seqDict = VCFFileReader.getSequenceDictionary(new File(vcfFile));
		
		 //get region
		if(region != null){
				String[] strTmp1 = region.split(":");
				if(strTmp1.length>1){
					String[] strTmp2 = strTmp1[1].split("-");
					itSnp = vcfReader.query(strTmp1[0], Integer.parseInt(strTmp2[0]), Integer.parseInt(strTmp2[1]));
					//System.err.println(strTmp1[0] + "\t" + Integer.parseInt(strTmp2[0]) + "\t" + Integer.parseInt(strTmp2[1]));
				}else{
					//process the whole chromosome;
					String chr = strTmp1[0];
					SAMSequenceRecord sm = seqDict.getSequence(chr);
					if(sm == null)
						throw new CodecLineParsingException("VCF file does not contain this sequence contig: " + chr);
					int seqEnd = sm.getSequenceLength();
					itSnp = vcfReader.query(chr, 1, seqEnd);
				}
			
		}else{
			itSnp = vcfReader.iterator();
		}
		GenomeLocParser glp = new GenomeLocParser(seqDict);
		GenomeLocSortedSet gss = new GenomeLocSortedSet(glp);
		TreeMap<GenomeLoc, Pair<Byte, Byte>> genotypeMapList = new TreeMap<GenomeLoc, Pair<Byte, Byte>>();
		while(itSnp.hasNext()){
			VariantContext vc = itSnp.next();
			//if(vc.getEnd()==221974){
			//	System.err.println(vc + "\t" + vc.isFiltered() + "\t" + vc.isIndel() + "\t" +  vc.isMixed());
			//}
				
			if(vc.isFiltered() || vc.isIndel() || vc.isMixed() || !vc.isPointEvent())
				continue;
			if(partialVcf){
				GenomeLoc genomeLoc = glp.createGenomeLoc(vc.getChr(),vc.getStart(), vc.getEnd());
				
				Byte alleleA = vc.getAlleles().get(0).getBases()[0];
				Byte alleleB = vc.getAlleles().get(1).getBases()[0];
				Pair<Byte, Byte> tmp = new Pair<Byte, Byte>(alleleA, alleleB);
				genotypeMapList.put(genomeLoc, tmp);
				gss.add(genomeLoc);
			}else{
				Genotype gt = vc.getGenotype(sample);
				//if(vc.getEnd()==221974)
					//System.err.println(vc + "\t" + vc.isFiltered() + "\t" + vc.isIndel() + "\t" +  vc.isMixed() + "\t" + gt.isFiltered() + "\t" + gt.isHet());
				if(!gt.isFiltered() && gt.isHet() && gt.isPhased()){ //only phased Het SNP is considered.
					GenomeLoc genomeLoc = glp.createGenomeLoc(vc.getChr(),vc.getStart(), vc.getEnd());
					
					Byte alleleA = gt.getAllele(0).getBases()[0]; // the first one (alleleA) is indeed the first at genotype. e,g, 1|0, the first is indeed the alternative allele. not the reference allele.
					Byte alleleB = gt.getAllele(1).getBases()[0];
					//System.err.println(vc.getChr() + "\t" + vc.getStart() + "\t" + vc.getEnd() + "\t" + String.valueOf(alleleA) + "\t" + String.valueOf(alleleB) + "\t" + vc.toStringDecodeGenotypes());
					Pair<Byte, Byte> tmp = new Pair<Byte, Byte>(alleleA, alleleB);
					if(!gss.contains(genomeLoc)){
						genotypeMapList.put(genomeLoc, tmp);
						gss.add(genomeLoc);
					}
					
				}
			}
			
		}
		itSnp.close();
		vcfReader.close();
		return new GenotypeList(gss, genotypeMapList, glp);
	}
	
	
	private void processBam(String bamFile, GenotypeList gl) throws Exception {
		File bamF = new File(bamFile);
		if(!BamFileIoUtils.isBamFile(bamF))
			throw new SAMException ("Not a Bam file: " + bamFile);
		SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(bamF);
		
		String preChr = "";
		HashMap<Integer,Pair<Integer, Integer>> methyRef = new HashMap<Integer,Pair<Integer, Integer>>();
		HashMap<Integer,Pair<Integer, Integer>> methyAlt = new HashMap<Integer,Pair<Integer, Integer>>();
		
		
		for(GenomeLoc loc : gl.genomeLocsAlleleMap.keySet()){
			if(VCFCOUNT % PURGE_INTERVAL == 0 ){
				writer.flush();
				refWriter.flush();
				altWriter.flush();
				log.info("Processing allele ... " + VCFCOUNT);
			}
				
			VCFCOUNT++;
			
			Pair<Byte, Byte> alleles = gl.genomeLocsAlleleMap.get(loc);
			
			if(preChr.equalsIgnoreCase("")){
				preChr=loc.getContig();
			}
			
			if(!preChr.equalsIgnoreCase(loc.getContig())){
				
					for(int refPos : methyRef.keySet()){
						Pair<Integer, Integer> v = methyRef.get(refPos);
						refWriter.write(preChr + "\t" + (refPos-1) + "\t" + refPos + "\t" + (double)v.getFirst()/(double)v.getSecond() + "\t" + v.getFirst() + "\t" + v.getSecond() + "\n");
					}
					for(int refPos : methyAlt.keySet()){
						Pair<Integer, Integer> v = methyAlt.get(refPos);
						altWriter.write(preChr + "\t" + (refPos-1) + "\t" + refPos + "\t" + (double)v.getFirst()/(double)v.getSecond() + "\t" + v.getFirst() + "\t" + v.getSecond() + "\n");
					}
				
				preChr=loc.getContig();
			}
			
			//System.err.println(loc.toString() + "\t" + loc.getStart() + "\t" + loc.getStop());
			SAMRecordIterator samIt = reader.queryOverlapping(loc.getContig(), loc.getStart(), loc.getStop());
			//record numC, numT of CG in each reads and their mate
			
			HashMap<String, SAMRecord> countedReadsEnd1 = new HashMap<String, SAMRecord>();	
			while(samIt.hasNext()){
				SAMRecord r = samIt.next();
				
				if(failFlagFilter(r)){
					continue;
				}else{
					if(countedReadsEnd1.containsKey(r.getReadName())){
						//countedReadsEnd1.remove(r.getReadName()); //if two ends are overlapped with the same SNPs, just randomly keep one of them there
					}else{
						countedReadsEnd1.put(r.getReadName(), r);
					}
				}
				
			}
			samIt.close();

			Pair<Integer, Integer> ctAlleleLocalRef = new Pair<Integer, Integer>(0,0);
			Pair<Integer, Integer> ctAlleleLocalAlt = new Pair<Integer, Integer>(0,0);

			int refAlleleLocalCount = 0;
			int altAlleleLocalCount = 0;

			Pair<Integer, Integer> ctAlleleRef = new Pair<Integer, Integer>(0,0);
			Pair<Integer, Integer> ctAlleleAlt = new Pair<Integer, Integer>(0,0);
			Pair<Integer, Integer> ctAlleleRefEnd2 = new Pair<Integer, Integer>(0,0);
			Pair<Integer, Integer> ctAlleleAltEnd2 = new Pair<Integer, Integer>(0,0);
			
			int refAlleleCount = 0;
			int altAlleleCount = 0;
			HashSet<Integer> end2Anchors = new HashSet<Integer>();
			int maxEnd2Dist = -1;

			//for 2nd methylation patterns
			Pair<Integer, Integer> ctAllele2ndLocalRef = new Pair<Integer, Integer>(0,0);
			Pair<Integer, Integer> ctAllele2ndLocalAlt = new Pair<Integer, Integer>(0,0);
			Pair<Integer, Integer> ctAllele2ndRef = new Pair<Integer, Integer>(0,0);
			Pair<Integer, Integer> ctAllele2ndAlt = new Pair<Integer, Integer>(0,0);
			Pair<Integer, Integer> ctAllele2ndRefEnd2 = new Pair<Integer, Integer>(0,0);
			Pair<Integer, Integer> ctAllele2ndAltEnd2 = new Pair<Integer, Integer>(0,0);
			
			for(String readName : countedReadsEnd1.keySet()){
				SAMRecord d = countedReadsEnd1.get(readName);
				SAMRecord r2 = reader.queryMate(d);
				boolean badr2 = failFlagFilter(r2);

					boolean negStrand = d.getReadNegativeStrandFlag();
					boolean secondPair = d.getReadPairedFlag() && d.getSecondOfPairFlag();
					
					byte[] readBasesQ = BisulfiteHicUtils.getClippedReadsBaseQuality(d);
					byte[] readBases = CcInferenceUtils.toUpperCase(BisulfiteHicUtils.getClippedReadsBase(d));
					
					//define which allele the reads belonged to
					
					//in C/T or A/G SNP only consider bisulfite uncoverted strand
					if((SequenceUtil.basesEqual(alleles.getFirst(), SequenceUtil.C) && SequenceUtil.basesEqual(alleles.getSecond(), SequenceUtil.T))
							|| (SequenceUtil.basesEqual(alleles.getFirst(), SequenceUtil.T) && SequenceUtil.basesEqual(alleles.getSecond(), SequenceUtil.C))
							){
						if((!secondPair && !negStrand) || (secondPair && negStrand)){
							continue;
						}
					}else if((SequenceUtil.basesEqual(alleles.getFirst(), SequenceUtil.G) && SequenceUtil.basesEqual(alleles.getSecond(), SequenceUtil.A))
							|| (SequenceUtil.basesEqual(alleles.getFirst(), SequenceUtil.A) && SequenceUtil.basesEqual(alleles.getSecond(), SequenceUtil.G))){
						if((!secondPair && negStrand) || (secondPair && !negStrand)){
							continue;
						}
					}
					
					//int pos = loc.getStart() - d.getAlignmentStart();
					int pos = d.getReadPositionAtReferencePosition(loc.getStart())-1;
					if(pos < 0 || pos >= readBases.length) //no such a postion in the clipped reads
						continue;
						
					byte snpBase = readBases[pos];
					if(readBasesQ[pos] < minBaseQ){
						continue;
					}
					//if(loc.getContig().equalsIgnoreCase("15") && loc.getStop() == 26068633){
						//System.err.println(d.getReadName() + "\t" + d.getAlignmentStart() + "\t" + d.getReadNegativeStrandFlag() + "\t" + d.getSecondOfPairFlag() + "\t" +  new String(readBases) + "\t" + pos + "\t" + (char)snpBase + "\t" + d.getBaseQualities()[pos]);
					//	System.err.println(d.getReadName() + "\t" +  new String(readBases) + "\t" +  new String(SequenceUtil.makeReferenceFromAlignment(d, false)));
					//	System.err.println(new String(BisSnp2Utils.getClippedReadBases(d)) + "\t" + d.getAlignmentStart() + "\t" + d.getAlignmentEnd() + "\t" + d.getUnclippedStart() + "\t" + d.getUnclippedEnd());
					//}
					int r1Start = d.getAlignmentStart();
					int r2Start = -1;
					if(!badr2){
						r2Start = r2.getAlignmentStart();
					}


					TreeMap<Integer, Integer> ctSummary = BisulfiteHicUtils.readsMethySummaryAtEachLocGeneral(d, methyPatternSearch, methyPatternCytPos, bsConvPatternSearch, bsConvPatternCytPos, minMapQ, minBaseQ, turnOffBisulfiteFilter, useBadMate);
					TreeMap<Integer, Integer> ctSummaryEnd2 = new TreeMap<Integer, Integer>();
					if(!badr2){
						ctSummaryEnd2 = BisulfiteHicUtils.readsMethySummaryAtEachLocGeneral(r2, methyPatternSearch, methyPatternCytPos, bsConvPatternSearch, bsConvPatternCytPos, minMapQ, minBaseQ, turnOffBisulfiteFilter, useBadMate);
					}
					if(ctSummary.isEmpty())
						continue;

					TreeMap<Integer, Integer> ctSummary2nd = new TreeMap<Integer, Integer>();
					TreeMap<Integer, Integer> ctSummary2ndEnd2 = new TreeMap<Integer, Integer>();
					if(methyPattern2ndSearch != null){
						ctSummary2nd = BisulfiteHicUtils.readsMethySummaryAtEachLocGeneral(d, methyPattern2ndSearch, methyPattern2ndCytPos, bsConvPatternSearch, bsConvPatternCytPos, minMapQ, minBaseQ, turnOffBisulfiteFilter, useBadMate);
						if(!badr2){
							ctSummary2ndEnd2 = BisulfiteHicUtils.readsMethySummaryAtEachLocGeneral(r2, methyPattern2ndSearch, methyPattern2ndCytPos, bsConvPatternSearch, bsConvPatternCytPos, minMapQ, minBaseQ, turnOffBisulfiteFilter, useBadMate);
						}
					}

					boolean refStrand = BisulfiteHicUtils.baseMatchInBisulfiteSpace(alleles.getFirst(), snpBase, negStrand, secondPair);
					boolean altStrand = BisulfiteHicUtils.baseMatchInBisulfiteSpace(alleles.getSecond(), snpBase, negStrand, secondPair);
					int numC = 0;
					int numCT = 0;
					int numCLocal = 0;
					int numCTLocal = 0;
					for(int rpos : ctSummary.keySet()){
						int methyStat = ctSummary.get(rpos);
						int cpos = r1Start + rpos;
						numCLocal += methyStat;
						numCTLocal++;
						if(!badr2){
							numC += methyStat;
							numCT++;
						}

						if(refStrand){
							if(methyRef.containsKey(cpos)){
								Pair<Integer, Integer> stat = methyRef.get(cpos);
								methyRef.put(cpos, new Pair<Integer, Integer>(stat.getFirst() + methyStat,stat.getSecond() + 1));
								
							}else{
								methyRef.put(cpos, new Pair<Integer, Integer>(methyStat,1));
							}
						}else if(altStrand){
							if(methyAlt.containsKey(cpos)){
								Pair<Integer, Integer> stat = methyAlt.get(cpos);
								methyAlt.put(cpos, new Pair<Integer, Integer>(stat.getFirst() + methyStat,stat.getSecond() + 1));
								
							}else{
								methyAlt.put(cpos, new Pair<Integer, Integer>(methyStat,1));
							}
						}
					}
					
					int numCEnd2 = 0;
					int numCTEnd2 = 0;
					if(!badr2){
						for(int rpos : ctSummaryEnd2.keySet()){
							int methyStat = ctSummaryEnd2.get(rpos);
							int cpos = r2Start + rpos;
							numCEnd2 += methyStat;
							numCTEnd2++;
							if(refStrand){
								if(methyRef.containsKey(cpos)){
									Pair<Integer, Integer> stat = methyRef.get(cpos);
									methyRef.put(cpos, new Pair<Integer, Integer>(stat.getFirst() + methyStat,stat.getSecond() + 1));

								}else{
									methyRef.put(cpos, new Pair<Integer, Integer>(methyStat,1));
								}
							}else if(altStrand){
								if(methyAlt.containsKey(cpos)){
									Pair<Integer, Integer> stat = methyAlt.get(cpos);
									methyAlt.put(cpos, new Pair<Integer, Integer>(stat.getFirst() + methyStat,stat.getSecond() + 1));

								}else{
									methyAlt.put(cpos, new Pair<Integer, Integer>(methyStat,1));
								}
							}
						}
					}

				int numC2nd = 0;
				int numCT2nd = 0;
				int numC2ndLocal = 0;
				int numCT2ndLocal = 0;
				int numC2ndEnd2 = 0;
				int numCT2ndEnd2 = 0;
					if(methyPattern2ndSearch != null){

						for(int rpos : ctSummary2nd.keySet()){
							int methyStat = ctSummary2nd.get(rpos);
							numC2ndLocal += methyStat;
							numCT2ndLocal++;
							if(!badr2){
								numC2nd += methyStat;
								numCT2nd++;
							}
						}


						if(!badr2){
							for(int rpos : ctSummary2ndEnd2.keySet()){
								int methyStat = ctSummary2ndEnd2.get(rpos);
								numC2ndEnd2 += methyStat;
								numCT2ndEnd2++;
							}
						}
					}

					
					if(refStrand){
						ctAlleleRef = new Pair<Integer, Integer>(ctAlleleRef.getFirst() + numC,ctAlleleRef.getSecond() + numCT);
						ctAlleleRefEnd2 = new Pair<Integer, Integer>(ctAlleleRefEnd2.getFirst() + numCEnd2,ctAlleleRefEnd2.getSecond() + numCTEnd2);
						ctAlleleLocalRef = new Pair<Integer, Integer>(ctAlleleLocalRef.getFirst() + numCLocal,ctAlleleLocalRef.getSecond() + numCTLocal);
						refAlleleLocalCount++;
						if(!badr2){
							refAlleleCount++;
						}
						if(methyPattern2ndSearch != null){
							ctAllele2ndRef = new Pair<Integer, Integer>(ctAllele2ndRef.getFirst() + numC2nd,ctAllele2ndRef.getSecond() + numCT2nd);
							ctAllele2ndRefEnd2 = new Pair<Integer, Integer>(ctAllele2ndRefEnd2.getFirst() + numC2ndEnd2,ctAllele2ndRefEnd2.getSecond() + numCT2ndEnd2);
							ctAllele2ndLocalRef = new Pair<Integer, Integer>(ctAllele2ndLocalRef.getFirst() + numC2ndLocal,ctAllele2ndLocalRef.getSecond() + numCT2ndLocal);
						}
					}else if(altStrand){
						ctAlleleAlt = new Pair<Integer, Integer>(ctAlleleAlt.getFirst() + numC,ctAlleleAlt.getSecond() + numCT);
						ctAlleleAltEnd2 = new Pair<Integer, Integer>(ctAlleleAltEnd2.getFirst() + numCEnd2,ctAlleleAltEnd2.getSecond() + numCTEnd2);
						ctAlleleLocalAlt = new Pair<Integer, Integer>(ctAlleleLocalAlt.getFirst() + numCLocal,ctAlleleLocalAlt.getSecond() + numCTLocal);
						altAlleleLocalCount++;
						if(!badr2){
							altAlleleCount++;
						}
						if(methyPattern2ndSearch != null){
							ctAllele2ndAlt = new Pair<Integer, Integer>(ctAllele2ndAlt.getFirst() + numC2nd,ctAllele2ndAlt.getSecond() + numCT2nd);
							ctAllele2ndAltEnd2 = new Pair<Integer, Integer>(ctAllele2ndAltEnd2.getFirst() + numC2ndEnd2,ctAllele2ndAltEnd2.getSecond() + numCT2ndEnd2);
							ctAllele2ndLocalAlt = new Pair<Integer, Integer>(ctAllele2ndLocalAlt.getFirst() + numC2ndLocal,ctAllele2ndLocalAlt.getSecond() + numCT2ndLocal);
						}
					}
					if(!badr2){
						if(end2Anchors.isEmpty()){
							end2Anchors.add(r2Start);
							maxEnd2Dist = Math.abs(r2Start-loc.getStart());
						}else{
							boolean addEntry = true;
							for(int p : end2Anchors){
								if(Math.abs(p-r2Start)<500){ //too far away from original, add it
									addEntry = false;
									break;
								}

							}
							if(addEntry){
								end2Anchors.add(r2Start);
								if(Math.abs(r2Start-loc.getStart()) > maxEnd2Dist){
									maxEnd2Dist = Math.abs(r2Start-loc.getStart());
								}
							}
						}
					}


			}

			//if there is no significant amount of long-range interacted reads, it will skip this SNPs..
			if(refAlleleLocalCount < coverageRef || altAlleleLocalCount< coverageAlt || (ctAlleleLocalRef.getSecond() ==0 &&ctAlleleLocalAlt.getSecond() ==0) ) //for the end 2, let the result figure out.. when no CG in the end2, it could still be useful to used for the ASE linking
				continue;
			String end2AchorString = "";
			for(int p : end2Anchors){
				if(end2AchorString.isEmpty()){
					end2AchorString = end2AchorString + p;
				}else{
					end2AchorString = end2AchorString + ":" + p ;
				}
				
			}

			if(methyPattern2ndSearch == null){
				writer.write(loc.getContig() + "\t" + loc.getStart() + "\t" + loc.getStop()  + "\t" + String.format("%c", alleles.getFirst()) + "\t" + refAlleleCount
						+ "\t" + ctAlleleRef.getFirst() + "\t" + ctAlleleRef.getSecond() + "\t" + (ctAlleleRef.getSecond() == 0 ? Double.NaN : (double)ctAlleleRef.getFirst()/(double)ctAlleleRef.getSecond())
						+ "\t" + String.format("%c", alleles.getSecond()) + "\t" + altAlleleCount
						+ "\t" + ctAlleleAlt.getFirst() + "\t" + ctAlleleAlt.getSecond() + "\t" + (ctAlleleAlt.getSecond() == 0 ? Double.NaN : (double)ctAlleleAlt.getFirst()/(double)ctAlleleAlt.getSecond())
						+ "\t" + end2AchorString + "\t" + ctAlleleRefEnd2.getFirst() + "\t" + ctAlleleRefEnd2.getSecond() + "\t" + (ctAlleleRefEnd2.getSecond()==0 ? Double.NaN : (double)ctAlleleRefEnd2.getFirst()/(double)ctAlleleRefEnd2.getSecond())
						+ "\t" + ctAlleleAltEnd2.getFirst() + "\t" + ctAlleleAltEnd2.getSecond() + "\t" + (ctAlleleAltEnd2.getSecond()==0 ? Double.NaN : (double)ctAlleleAltEnd2.getFirst()/(double)ctAlleleAltEnd2.getSecond()) + "\t" + maxEnd2Dist
						+ "\t" + refAlleleLocalCount + "\t" + ctAlleleLocalRef.getFirst() + "\t" + ctAlleleLocalRef.getSecond() + "\t" + (ctAlleleLocalRef.getSecond() == 0 ? Double.NaN : (double)ctAlleleLocalRef.getFirst()/(double)ctAlleleLocalRef.getSecond())
						+ "\t" + altAlleleLocalCount + "\t" + ctAlleleLocalAlt.getFirst() + "\t" + ctAlleleLocalAlt.getSecond() + "\t" + (ctAlleleLocalAlt.getSecond() == 0 ? Double.NaN : (double)ctAlleleLocalAlt.getFirst()/(double)ctAlleleLocalAlt.getSecond()) + "\n");
			}else{
				writer.write(loc.getContig() + "\t" + loc.getStart() + "\t" + loc.getStop()  + "\t" + String.format("%c", alleles.getFirst()) + "\t" + refAlleleCount
						+ "\t" + ctAlleleRef.getFirst() + "\t" + ctAlleleRef.getSecond() + "\t" + (ctAlleleRef.getSecond() == 0 ? Double.NaN : (double)ctAlleleRef.getFirst()/(double)ctAlleleRef.getSecond())
						+ "\t" + String.format("%c", alleles.getSecond()) + "\t" + altAlleleCount
						+ "\t" + ctAlleleAlt.getFirst() + "\t" + ctAlleleAlt.getSecond() + "\t" + (ctAlleleAlt.getSecond() == 0 ? Double.NaN : (double)ctAlleleAlt.getFirst()/(double)ctAlleleAlt.getSecond())
						+ "\t" + end2AchorString + "\t" + ctAlleleRefEnd2.getFirst() + "\t" + ctAlleleRefEnd2.getSecond() + "\t" + (ctAlleleRefEnd2.getSecond()==0 ? Double.NaN : (double)ctAlleleRefEnd2.getFirst()/(double)ctAlleleRefEnd2.getSecond())
						+ "\t" + ctAlleleAltEnd2.getFirst() + "\t" + ctAlleleAltEnd2.getSecond() + "\t" + (ctAlleleAltEnd2.getSecond()==0 ? Double.NaN : (double)ctAlleleAltEnd2.getFirst()/(double)ctAlleleAltEnd2.getSecond()) + "\t" + maxEnd2Dist
						+ "\t" + refAlleleLocalCount + "\t" + ctAlleleLocalRef.getFirst() + "\t" + ctAlleleLocalRef.getSecond() + "\t" + (ctAlleleLocalRef.getSecond() == 0 ? Double.NaN : (double)ctAlleleLocalRef.getFirst()/(double)ctAlleleLocalRef.getSecond())
						+ "\t" + altAlleleLocalCount + "\t" + ctAlleleLocalAlt.getFirst() + "\t" + ctAlleleLocalAlt.getSecond() + "\t" + (ctAlleleLocalAlt.getSecond() == 0 ? Double.NaN : (double)ctAlleleLocalAlt.getFirst()/(double)ctAlleleLocalAlt.getSecond())
						+ "\t" + ctAllele2ndLocalRef.getFirst() + "\t" + ctAllele2ndLocalRef.getSecond() + "\t" + (ctAllele2ndLocalRef.getSecond() == 0 ? Double.NaN : (double)ctAllele2ndLocalRef.getFirst()/(double)ctAllele2ndLocalRef.getSecond())
						+ "\t" + ctAllele2ndLocalAlt.getFirst() + "\t" + ctAllele2ndLocalAlt.getSecond() + "\t" + (ctAllele2ndLocalAlt.getSecond() == 0 ? Double.NaN : (double)ctAllele2ndLocalAlt.getFirst()/(double)ctAllele2ndLocalAlt.getSecond())
						+ "\t" + ctAllele2ndRef.getFirst() + "\t" + ctAllele2ndRef.getSecond() + "\t" + (ctAllele2ndRef.getSecond() == 0 ? Double.NaN : (double)ctAllele2ndRef.getFirst()/(double)ctAllele2ndRef.getSecond())
						+ "\t" + ctAllele2ndAlt.getFirst() + "\t" + ctAllele2ndAlt.getSecond() + "\t" + (ctAllele2ndAlt.getSecond() == 0 ? Double.NaN : (double)ctAllele2ndAlt.getFirst()/(double)ctAllele2ndAlt.getSecond())
						+ "\t" + ctAllele2ndRefEnd2.getFirst() + "\t" + ctAllele2ndRefEnd2.getSecond() + "\t" + (ctAllele2ndRefEnd2.getSecond() == 0 ? Double.NaN : (double)ctAllele2ndRefEnd2.getFirst()/(double)ctAllele2ndRefEnd2.getSecond())
						+ "\t" + ctAllele2ndAltEnd2.getFirst() + "\t" + ctAllele2ndAltEnd2.getSecond() + "\t" + (ctAllele2ndAltEnd2.getSecond() == 0 ? Double.NaN : (double)ctAllele2ndAltEnd2.getFirst()/(double)ctAllele2ndAltEnd2.getSecond()) + "\n");
			}

			
			
		}
		for(int refPos : methyRef.keySet()){
			Pair<Integer, Integer> v = methyRef.get(refPos);
			refWriter.write(preChr + "\t" + (refPos-1) + "\t" + refPos + "\t" + (double)v.getFirst()/(double)v.getSecond() + "\t" + v.getFirst() + "\t" + v.getSecond() + "\n");
		}
		for(int refPos : methyAlt.keySet()){
			Pair<Integer, Integer> v = methyAlt.get(refPos);
			altWriter.write(preChr + "\t" + (refPos-1) + "\t" + refPos + "\t" + (double)v.getFirst()/(double)v.getSecond() + "\t" + v.getFirst() + "\t" + v.getSecond() + "\n");
		}
		reader.close();
	}


	private boolean failFlagFilter(SAMRecord r){
		return r.getReadUnmappedFlag() || r.getNotPrimaryAlignmentFlag() || r.getMappingQuality() < minMapQ
				|| r.getReadFailsVendorQualityCheckFlag() || r.getDuplicateReadFlag() || (useUnpaired ? false : !r.getReadPairedFlag()) || (useBadMate ? false : !r.getProperPairFlag()) || Math.abs(r.getInferredInsertSize()) > maxDist || (diffChr ? (Math.abs(r.getInferredInsertSize()) < minDist && Math.abs(r.getInferredInsertSize())!=0 ): Math.abs(r.getInferredInsertSize()) < minDist);
	}
	
	private VCFFileReader initiate(String vcfFile, String outputFile, String refFile, String altFile) throws IOException{
		startTime = System.currentTimeMillis();
		
		File indexFile = Tribble.tabixIndexFile(new File(vcfFile));
		if(!indexFile.exists() || !indexFile.canRead()){
			throw new UnsortedFileException(vcfFile + " file's index " + indexFile.getName() + " does not exist, please use tabix to index it ...");
		}
		
		writer =  new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFile)), "UTF-8");
		//write the file header
		if(methyPattern2ndSearch == null){
			writer.write("#chr\tstart\tend\tAlleleA\tAlleleA_count\tAlleleA_C\tAlleleA_T\tAlleleA_meth_mean\t"
					+ "AlleleB\tAlleleB_count\tAlleleB_C\tAlleleB_T\tAlleleB_meth_mean\tEnd2_anchor_point\tEnd2_AlleleA_C\tEnd2_AlleleA_T\tEnd2_AlleleA_meth_mean\t"
					+ "End2_AlleleB_C\tEnd2_AlleleB_T\tEnd2_AlleleB_meth_mean\tMaxEnd2Dist\t"
					+ "AlleleA_local_count\tAlleleA_local_C\tAlleleA_local_T\tAlleleA_local_meth_mean\t"
					+ "AlleleB_local_count\tAlleleB_local_C\tAlleleB_local_T\tAlleleB_local_meth_mean\n");
		}else{
			writer.write("#chr\tstart\tend\tAlleleA\tAlleleA_count\tAlleleA_C\tAlleleA_T\tAlleleA_meth_mean\t"
					+ "AlleleB\tAlleleB_count\tAlleleB_C\tAlleleB_T\tAlleleB_meth_mean\tEnd2_anchor_point\tEnd2_AlleleA_C\tEnd2_AlleleA_T\tEnd2_AlleleA_meth_mean\t"
					+ "End2_AlleleB_C\tEnd2_AlleleB_T\tEnd2_AlleleB_meth_mean\tMaxEnd2Dist\t"
					+ "AlleleA_local_count\tAlleleA_local_C\tAlleleA_local_T\tAlleleA_local_meth_mean\t"
					+ "AlleleB_local_count\tAlleleB_local_C\tAlleleB_local_T\tAlleleB_local_meth_mean\t"
					+ "AlleleA_2nd_local_C\tAlleleA_2nd_local_T\tAlleleA_2nd_local_meth_mean\t"
					+ "AlleleB_2nd_local_C\tAlleleB_2nd_local_T\tAlleleB_2nd_local_meth_mean\t"
					+ "AlleleA_2nd_C\tAlleleA_2nd_T\tAlleleA_2nd_meth_mean\t"
					+ "AlleleB_2nd_C\tAlleleB_2nd_T\tAlleleB_2nd_meth_mean\t"
					+ "End2_AlleleA_2nd_C\tEnd2_AlleleA_2nd_T\tEnd2_AlleleA_2nd_meth_mean\t"
					+ "End2_AlleleB_2nd_C\tEnd2_AlleleB_2nd_T\tEnd2_AlleleB_2nd_meth_mean\n");
		}

		refWriter =  new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(refFile)), "UTF-8");
		altWriter =  new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(altFile)), "UTF-8");
		
		 return new VCFFileReader(new File(vcfFile), true);
	}

	private void finish() throws IOException{
		writer.close();
		refWriter.close();
		altWriter.close();
		long endTime   = System.currentTimeMillis();
		double totalTime = endTime - startTime;
		totalTime /= 1000;
		double totalTimeMins = totalTime/60;
		double totalTimeHours = totalTime/3600;
		
		log.info("Processing " + VCFCOUNT  + " SNPs in total");
		log.info("LongRangeAsm's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
	}

}
