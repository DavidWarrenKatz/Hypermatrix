/**
 * LongRangeAsm.java
 * Feb 27, 2017
 * 12:00:36 PM
 * yaping    lyping1986@gmail.com
 */
package main.java.edu.mit.compbio.bisulfitehic.utils;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.exception.CodecLineParsingException;
import htsjdk.tribble.exception.UnsortedFileException;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import main.java.edu.mit.compbio.ccinference.utils.CcInferenceUtils;
import org.apache.commons.math3.stat.inference.AlternativeHypothesis;
import org.apache.commons.math3.util.Pair;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.stat.inference.BinomialTest;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.GenomeLocSortedSet;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;


/**
 * output if not merged snp within the same genes
 * chr(gene) start(gene) end(gene) gene_name . gene_strand SNP_pos alleleA alleleA_count alleleB alleleB_count
 * if merged
 * chr(gene) start(gene) end(gene) gene_name . gene_strand SNP_pos1:SNP_pos2:...:SNP_posN alleleA alleleA_count alleleB alleleB_count
 */
public class AlleImbalanceExpr {


	@Option(name="-chrom",usage="only treat the specific chromosome. default: null, the whole bam file")
	public String chrom = null;

	@Option(name="-mergedInSameGene",usage="merge the allelic imbalance summary statistics within the same gene.  default: false")
	public boolean mergedInSameGene = false;

	@Option(name="-minMapQ",usage="minmum mapping quality score, default: 30")
	public int minMapQ = 30;

	@Option(name="-minBaseQ",usage="minmum mapping quality score, default: 5")
	public int minBaseQ = 5;
	
	@Option(name="-minCov",usage="specify the minimum coverage in at least one of the SNP's allele , default: 5")
	public int minCov = 5;

	@Option(name="-qual",usage="specify the minimum genotyping quality of SNP call, default: 30")
	public int qual = 30;

	@Option(name="-partialVcf",usage="only use REF/ALT part of VCF rather than using genotype field in VCF file, default: false")
	public boolean partialVcf = false;

	@Option(name="-consistentPrefix",usage="check the consistent of chr prefix, basically, add chr at the beginning of GRch37 genome in VCF or bed gene file's location, default: false")
	public boolean consistentPrefix = false;

	@Option(name="-waspFilter",usage="use WASP filter when enabling --waspOutputMode SAMtag in STAR aligner. Basically, when reads have vW flag, only keep reads with vW:i:1, default: false")
	public boolean waspFilter = false;

	@Option(name="-sample",usage="specify the sample number in the VCF file if contain multiple sample genotype. number started at 0. default: 0")
	public int sample = 0;


	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "AlleImbalanceExpr [opts] gene.bed[.gz] genotype.vcf.gz input.bam output.ASE.tab.gz";
	
	private static Logger log = Logger.getLogger(AlleImbalanceExpr.class);

	private static long startTime = -1;
	
	private OutputStreamWriter writer = null;
	
	private final static int PURGE_INTERVAL = 1000;
	
	private long VCFCOUNT = 0L;

	private final static String WASP_FLAG = "vW";
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		AlleImbalanceExpr lac = new AlleImbalanceExpr();
		BasicConfigurator.configure();
	    lac.doMain(args);

	}
	
	public void doMain(String[] args)
			throws Exception {

					CmdLineParser parser = new CmdLineParser(this);
					//parser.setUsageWidth(80);
					try
					{
						if(help || args.length < 4) throw new CmdLineException(USAGE);
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

					String geneFile = arguments.get(0);
					String vcfFile = arguments.get(1);
					String bamFile = arguments.get(2);
					String outputFile = arguments.get(3);

					VCFFileReader vcfReader = initiate(vcfFile, outputFile);
					SAMSequenceDictionary seqDict = VCFFileReader.getSequenceDictionary(new File(vcfFile));
					GenomeLocParser glp = new GenomeLocParser(seqDict);

					log.info("Parsing input gene bed file ...");
					HashMap<GenomeLoc, String> geneList = new HashMap<GenomeLoc, String>();
					GZIPInputStream gzipInputStream = null;
					BufferedReader br;
					if(geneFile.endsWith(".gz")){
							gzipInputStream = new GZIPInputStream(new FileInputStream(geneFile));
							br = new BufferedReader(new InputStreamReader(gzipInputStream));

					}else{
							br = new BufferedReader(new FileReader(geneFile));
					}
					String line;

					while( (line = br.readLine()) != null) {
							if (line.startsWith("#")) {
								continue;
							} else {
								String[] splitin = line.split("\t");
								String chr = splitin[0];
								int start = Integer.parseInt(splitin[1]);
								int end = Integer.parseInt(splitin[2]);

								//only deal with the regions...
								if(chrom != null && !chrom.equalsIgnoreCase(chr)){
									continue;
								}
								GenomeLoc genomeLoc = glp.createGenomeLoc(chr, start, end);
								geneList.put(genomeLoc, line);
							}
					}
					if(geneFile.endsWith(".gz")){
						gzipInputStream.close();
					}
					br.close();

					log.info("Parsing input vcf file ...");
					GenotypeList gl = readSnp(vcfReader, seqDict, glp);

					log.info("Parsing bam file ...");
					//System.err.println(gl.genomeLocsAlleleMap.size());
					processBam(bamFile, gl, geneList);



					finish();
	}
	

	
	private GenotypeList readSnp(VCFFileReader vcfReader, SAMSequenceDictionary seqDict, GenomeLocParser glp){
		CloseableIterator<VariantContext> itSnp = null;

		
		 //get region
		if(chrom != null){
					SAMSequenceRecord sm = seqDict.getSequence(chrom);
					if(sm == null)
						throw new CodecLineParsingException("VCF file does not contain this sequence contig: " + chrom);
					int seqEnd = sm.getSequenceLength();
					itSnp = vcfReader.query(chrom, 1, seqEnd);
		}else{
			itSnp = vcfReader.iterator();
		}

		GenomeLocSortedSet gss = new GenomeLocSortedSet(glp);
		TreeMap<GenomeLoc, Triple<Byte, Byte, Byte>> genotypeMapList = new TreeMap<GenomeLoc, Triple<Byte, Byte, Byte>>();
		while(itSnp.hasNext()){
			VariantContext vc = itSnp.next();
			//if(vc.getEnd()==221974){
			//	System.err.println(vc + "\t" + vc.isFiltered() + "\t" + vc.isIndel() + "\t" +  vc.isMixed());
			//}
				
			if(vc.isFiltered() || vc.isIndel() || vc.isMixed() || !vc.isPointEvent())
				continue;
			if(partialVcf){

				GenomeLoc genomeLoc = glp.createGenomeLoc(vc.getChr(), vc.getStart(), vc.getEnd());
				Byte alleleA = vc.getAlleles().get(0).getBases()[0];
				Byte alleleB = vc.getAlleles().get(1).getBases()[0];
				Byte alleleRef = vc.getReference().getBases()[0];
				Triple<Byte, Byte, Byte> tmp = Triple.of(alleleA, alleleB, alleleRef);
				genotypeMapList.put(genomeLoc, tmp);
				gss.add(genomeLoc);
			}else{
				Genotype gt = vc.getGenotype(sample);
				//if(vc.getEnd()==221974)
					//System.err.println(vc + "\t" + vc.isFiltered() + "\t" + vc.isIndel() + "\t" +  vc.isMixed() + "\t" + gt.isFiltered() + "\t" + gt.isHet());
				if(!gt.isFiltered() && gt.isHet() && gt.isPhased()){ //only phased Het SNP is considered.
					GenomeLoc genomeLoc = glp.createGenomeLoc(vc.getChr(), vc.getStart(), vc.getEnd());
					
					Byte alleleA = gt.getAllele(0).getBases()[0]; // the first one (alleleA) is indeed the first at genotype. e,g, 1|0, the first is indeed the alternative allele. not the reference allele.
					Byte alleleB = gt.getAllele(1).getBases()[0];
					Byte alleleRef = vc.getReference().getBases()[0];
					//System.err.println(vc.getChr() + "\t" + vc.getStart() + "\t" + vc.getEnd() + "\t" + String.valueOf(alleleA) + "\t" + String.valueOf(alleleB) + "\t" + vc.toStringDecodeGenotypes());
					Triple<Byte, Byte, Byte> tmp = Triple.of(alleleA, alleleB, alleleRef);
					if(!gss.contains(genomeLoc)){
						genotypeMapList.put(genomeLoc, tmp);
						gss.add(genomeLoc);
					}
					
				}
			}
			
		}
		itSnp.close();
		vcfReader.close();

		return new GenotypeList(gss, glp, genotypeMapList);
	}
	
	
	private void processBam(String bamFile, GenotypeList gl, HashMap<GenomeLoc, String> geneList) throws Exception {
		File bamF = new File(bamFile);
		if(!BamFileIoUtils.isBamFile(bamF))
			throw new SAMException ("Not a Bam file: " + bamFile);
		SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(bamF);
		BinomialTest bTest = new BinomialTest();
		for(GenomeLoc locGene : geneList.keySet()){
			String geneString = geneList.get(locGene);
			List<GenomeLoc> locsOveralpped = gl.genomeLocsCol.getOverlapping(locGene);
			int alleleRefCount = 0;
			int alleleACount = 0;
			int alleleBCount = 0;
			int alleleNCount = 0;
			String allelePos = "";
			String alleleAStr = "";
			String alleleBStr = "";
			String alleleRefStr = "";
			for(GenomeLoc loc : locsOveralpped){
				if(VCFCOUNT % PURGE_INTERVAL == 0 ){
					writer.flush();
					log.info("Processing allele ... " + VCFCOUNT);
				}
				if(!mergedInSameGene){
					alleleRefCount = 0;
					alleleACount = 0;
					alleleBCount = 0;
					alleleNCount = 0;
					allelePos = "";
					alleleAStr = "";
					alleleBStr = "";
					alleleRefStr = "";
				}
				VCFCOUNT++;

				Triple<Byte, Byte, Byte> alleles = gl.genomeLocsAlleleMapWithRef.get(loc);

				SAMRecordIterator samIt = reader.queryOverlapping(consistentPrefix ? (loc.getContig().startsWith("chr") ? loc.getContig() : ("chr" + loc.getContig()) ): loc.getContig(), loc.getStart(), loc.getStop());

				HashMap<String, SAMRecord> countedReads = new HashMap<String, SAMRecord>();
				while(samIt.hasNext()){
					SAMRecord r = samIt.next();

					if(failFlagFilter(r)){
						continue;
					}else{
						countedReads.put(r.getReadName(), r); //two ends overlapped with the same SNPs are counted once only.

					}

				}
				samIt.close();

				for(String readName : countedReads.keySet()){
					SAMRecord d = countedReads.get(readName);

					byte[] readBasesQ = BisulfiteHicUtils.getClippedReadsBaseQuality(d);
					byte[] readBases = CcInferenceUtils.toUpperCase(BisulfiteHicUtils.getClippedReadsBase(d));
					//int pos = loc.getStart() - d.getAlignmentStart();
					int pos = d.getReadPositionAtReferencePosition(loc.getStart())-1;
					if(pos < 0 || pos >= readBases.length) //no such a postion in the clipped reads
						continue;

					byte snpBase = readBases[pos];
					if(readBasesQ[pos] < minBaseQ){
						continue;
					}
					boolean AlleleA = SequenceUtil.basesEqual(alleles.getLeft(), snpBase);
					boolean AlleleB = SequenceUtil.basesEqual(alleles.getMiddle(), snpBase);
					boolean AlleleRef = SequenceUtil.basesEqual(alleles.getRight(), snpBase);

					if(AlleleA){
						alleleACount++;
					}else if(AlleleB){
						alleleBCount++;
					}else{
						alleleNCount++;
					}

					if(AlleleRef){
						alleleRefCount++;
					}
				}
				if(!mergedInSameGene){ // write each SNP position.
					if(alleleACount < minCov && alleleBCount< minCov)
						continue;
					allelePos = String.format("%d", loc.getStart());
					alleleAStr = String.format("%c", alleles.getLeft());
					alleleBStr = String.format("%c", alleles.getMiddle());
					alleleRefStr = String.format("%c", alleles.getRight());
					double pvalue = bTest.binomialTest(alleleACount+alleleBCount, alleleACount,0.5, AlternativeHypothesis.TWO_SIDED);
					writer.write(geneString + "\t" + allelePos  + "\t" + alleleRefStr + "\t" + alleleRefCount + "\t" + alleleAStr + "\t" + alleleACount
							+ "\t" + alleleBStr + "\t" + alleleBCount + "\t" + String.format("%.3f", pvalue) + "\t" + alleleNCount + "\n");
				}else{
					if(allelePos.isEmpty()){
						allelePos = String.format("%d", loc.getStart());
						alleleAStr = String.format("%c", alleles.getLeft());
						alleleBStr = String.format("%c", alleles.getMiddle());
						alleleRefStr = String.format("%c", alleles.getRight());
					}else{
						allelePos = allelePos + ":" + String.format("%d", loc.getStart());
						alleleAStr = alleleAStr + ":" + String.format("%c", alleles.getLeft());
						alleleBStr = alleleBStr + ":" + String.format("%c", alleles.getMiddle());
						alleleRefStr = alleleRefStr + ":" + String.format("%c", alleles.getRight());
					}
				}
			}
			if(mergedInSameGene) { // write each gene position.
				if(alleleACount < minCov && alleleBCount< minCov)
					continue;
				double pvalue = bTest.binomialTest(alleleACount+alleleBCount, alleleACount,0.5, AlternativeHypothesis.TWO_SIDED);
				writer.write(geneString + "\t" + allelePos  + "\t" + alleleACount
						+ "\t" + alleleBCount + "\t" + String.format("%.3f", pvalue) + "\t" + alleleNCount + "\n");
			}
		}

		reader.close();
	}


	private boolean failFlagFilter(SAMRecord r){
		return r.getReadUnmappedFlag() || r.getNotPrimaryAlignmentFlag() || r.getMappingQuality() < minMapQ
				|| r.getReadFailsVendorQualityCheckFlag() || (waspFilter && r.getAttribute(WASP_FLAG)!=null && ((Integer)r.getAttribute(WASP_FLAG))!=1);
	}
	
	private VCFFileReader initiate(String vcfFile, String outputFile) throws IOException{
		startTime = System.currentTimeMillis();
		
		File indexFile = Tribble.tabixIndexFile(new File(vcfFile));
		if(!indexFile.exists() || !indexFile.canRead()){
			throw new UnsortedFileException(vcfFile + " file's index " + indexFile.getName() + " does not exist, please use tabix to index it ...");
		}
		
		writer =  new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFile)), "UTF-8");
		//write the file header
		if(mergedInSameGene){
			writer.write("#chr\tgene_start\tgene_end\tgene_name\t.\tstrand\tsnp_pos\tAlleleRef_count\tAlleleA_count\t"
					+ "AlleleB_count\tBinomial_p\tAlleleN_count\n");
		}else{
			writer.write("#chr\tgene_start\tgene_end\tgene_name\t.\tstrand\tsnp_pos\tAlleleRef\tAlleleRef_count\tAlleleA\tAlleleA_count\t"
					+ "AlleleB\tAlleleB_count\tBinomial_p\tAlleleN_count\n");
		}



		 return new VCFFileReader(new File(vcfFile), true);
	}

	private void finish() throws IOException{
		writer.close();

		long endTime   = System.currentTimeMillis();
		double totalTime = endTime - startTime;
		totalTime /= 1000;
		double totalTimeMins = totalTime/60;
		double totalTimeHours = totalTime/3600;
		
		log.info("Processing " + VCFCOUNT  + " SNPs in total");
		log.info("AlleImbalanceExpr's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
	}

}
