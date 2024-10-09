/**
 * BisulfiteHicUtils.java
 * Sep 13, 2016
 * 10:36:16 AM
 * yaping    lyping1986@gmail.com
 */
package main.java.edu.mit.compbio.bisulfitehic.utils;

import java.util.ArrayList;
import java.util.TreeMap;

import main.java.edu.mit.compbio.ccinference.utils.BaseUtilsMore;
import main.java.edu.mit.compbio.ccinference.utils.CcInferenceUtils;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.util.Pair;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CigarUtil;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.gatk.utils.BaseUtils;

/**
 *
 */
public class BisulfiteHicUtils {
	
	//it looks like Bhmem has some problem for the reads mapped to the negative strand... need to reverse it here...when in the future, this is correct in Bhmem, need to get rid of this reverse here
	static public byte[] getClippedReadsBase(SAMRecord r){
		int truncateStart = r.getReadPositionAtReferencePosition(r.getAlignmentStart())-1;
		int truncateEnd = r.getReadPositionAtReferencePosition(r.getAlignmentEnd());
		byte[] truncatedBases = new byte[truncateEnd-truncateStart];
		byte[] bases = r.getReadBases();
		//if(r.getReadNegativeStrandFlag()){ //TODO: remove the whole if, after using the correct bam file
		//	ArrayUtils.reverse(bases);
		//}
		
		for(int i = truncateStart, index=0; i < truncateEnd; i++, index++){
			truncatedBases[index] = bases[i];
		}
		
		return truncatedBases;
	}
	
	static public byte[] getClippedReadsBaseQuality(SAMRecord r){
		int truncateStart = r.getReadPositionAtReferencePosition(r.getAlignmentStart())-1;
		int truncateEnd = r.getReadPositionAtReferencePosition(r.getAlignmentEnd());
		byte[] truncatedBaseQs = new byte[truncateEnd-truncateStart];
		byte[] baseQs = r.getBaseQualities();
		//if(r.getReadNegativeStrandFlag()){//TODO: remove the whole if, after using the correct bam file
		//	ArrayUtils.reverse(baseQs);
		//}
		for(int i = truncateStart, index=0; i < truncateEnd; i++, index++){
			truncatedBaseQs[index] = baseQs[i];
		}
		
		return truncatedBaseQs;
	}
	
	static public byte[] complementArray(byte[] a){
		byte[] b = new byte[a.length];
		for(int i = 0; i < a.length; i++){
			b[i]=SequenceUtil.complement(a[i]);
		}
		return b;
	}
	
	static public String complementArray(String as){
		byte[] a = as.getBytes();
		byte[] b = new byte[a.length];
		for(int i = 0; i < a.length; i++){
			b[i]=SequenceUtil.complement(a[i]);
		}
		return new String(b);
	}
	
	//TODO: check the CIGAR string problem in original bwa... CIGAR should not start with I...
	static public byte[] modifyRefSeqByCigar(byte[] seqByte, String cigarString){ //mainly for "I" and "D" CIGAR string
		char[] cigarList = CigarUtil.cigarArrayFromString(cigarString);

		ArrayList<Byte> seqsNew = new ArrayList<Byte>();
		//System.err.println(new String(seqByte) + "\t" + seqByte.length + "\t" + cigarString);
		int offSet = 0;
		boolean startWithInsertion = true;
		for(char cigar : cigarList){
			
			if(cigar == 'M' || (cigar == 'I' && !startWithInsertion)){
				seqsNew.add(seqByte[offSet]);
				offSet++;
				startWithInsertion = false;
			}else if(cigar == 'D'){
				//seqsNew.add(SequenceUtil.N);
				
			}else if(cigar == 'S' || (cigar == 'I' && startWithInsertion)){
				offSet++;
			}
		}
		return ArrayUtils.toPrimitive(seqsNew.toArray(new Byte[seqsNew.size()]));
	}
	
	static public Triple<Integer,Integer, Integer> readsMethySummary(SAMRecord read, int minMapQ, int minBaseQ, boolean turnOffBisulfiteFilter) throws Exception{
		boolean negStrand = read.getReadNegativeStrandFlag();
		boolean secondPair = read.getReadPairedFlag() && read.getSecondOfPairFlag();
		//filter if it is bad reads
		if(read.getReadUnmappedFlag() || read.getNotPrimaryAlignmentFlag() || read.getDuplicateReadFlag() 
				|| read.getReadFailsVendorQualityCheckFlag() || (read.getReadPairedFlag() && !read.getProperPairFlag()) 
				|| read.getMappingQuality() < minMapQ){
			return Triple.of(0,0,0);
		}
		
		byte[] refBases = CcInferenceUtils.toUpperCase(modifyRefSeqByCigar(SequenceUtil.makeReferenceFromAlignment(read, false), read.getCigarString()));
		
		byte[] readBasesQ = getClippedReadsBaseQuality(read);
		byte[] readBases = CcInferenceUtils.toUpperCase(getClippedReadsBase(read));
		//if(negStrand){
			//System.err.println(new String(refBases) + "\n" + new String(readBases) + "\t" +read.getCigarString() + "\t" + secondPair);
			//readBases = complementArray(readBases);
		//}
		
		Triple<Integer,Integer, Integer> ctSummary = Triple.of(0,0,0);
		int pos = 0;
		if(!turnOffBisulfiteFilter){
			//System.err.println(new String(refBases) + "\t" + refBases.length + "\n" + new String(readBases) + "\t" + readBases.length + "\n" + new String(SequenceUtil.makeReferenceFromAlignment(read, false)) + "\t" + read.getCigarString());
			pos = CcInferenceUtils.bisulfiteIncompleteReads(readBases, refBases, negStrand, secondPair, false, true, true, 0.1, 0.4, 0, "CH", 1);
			if(pos<0){
				return ctSummary;
			}
		}
		
		ctSummary = readsMethySummaryWithFivePrimeFilter(refBases, readBases, readBasesQ, negStrand, secondPair, minBaseQ, pos);

		return ctSummary; //first is C, second is T
	}


	static public Triple<Integer,Integer, Integer> readsMethySummaryGeneral(SAMRecord read, String patSearchMeth, int cPos, String patCheckBsConv, int cPosBs, int minMapQ, int minBaseQ, boolean turnOffBisulfiteFilter, boolean useBadMate) throws Exception{
		boolean negStrand = read.getReadNegativeStrandFlag();
		boolean secondPair = read.getReadPairedFlag() && read.getSecondOfPairFlag();
		//filter if it is bad reads
		if(read.getReadUnmappedFlag() || read.getNotPrimaryAlignmentFlag() || read.getDuplicateReadFlag()
				|| read.getReadFailsVendorQualityCheckFlag() || !read.getReadPairedFlag() || (useBadMate ? false : !read.getProperPairFlag())
				|| read.getMappingQuality() < minMapQ){
			return Triple.of(0,0,0);
		}

		byte[] refBases = CcInferenceUtils.toUpperCase(modifyRefSeqByCigar(SequenceUtil.makeReferenceFromAlignment(read, false), read.getCigarString()));

		byte[] readBasesQ = getClippedReadsBaseQuality(read);
		byte[] readBases = CcInferenceUtils.toUpperCase(getClippedReadsBase(read));
		//if(negStrand){
		//System.err.println(new String(refBases) + "\n" + new String(readBases) + "\t" +read.getCigarString() + "\t" + secondPair);
		//readBases = complementArray(readBases);
		//}

		Triple<Integer,Integer, Integer> ctSummary = Triple.of(0,0,0);
		int pos = 0;
		if(!turnOffBisulfiteFilter){
			//System.err.println(new String(refBases) + "\t" + refBases.length + "\n" + new String(readBases) + "\t" + readBases.length + "\n" + new String(SequenceUtil.makeReferenceFromAlignment(read, false)) + "\t" + read.getCigarString());
			pos = CcInferenceUtils.bisulfiteIncompleteReads(readBases, refBases, negStrand, secondPair, false, true, true, 0.1, 0.4, 0, patCheckBsConv, cPosBs);
			if(pos<0){
				return ctSummary;
			}
		}
		byte[] baseSearchMeth = patSearchMeth.getBytes();
		ctSummary = readsMethySummaryWithFivePrimeFilterGeneral(refBases, readBases, readBasesQ, baseSearchMeth, cPos, negStrand, secondPair, minBaseQ, pos);

		return ctSummary; //first is C, second is T, third is number of CG
	}


	static public Triple<Integer,Integer, Integer> readsMethySummary(SAMRecord read, int minMapQ, int minBaseQ, boolean turnOffBisulfiteFilter, boolean useBadMate) throws Exception{
		boolean negStrand = read.getReadNegativeStrandFlag();
		boolean secondPair = read.getReadPairedFlag() && read.getSecondOfPairFlag();
		//filter if it is bad reads
		if(read.getReadUnmappedFlag() || read.getNotPrimaryAlignmentFlag() || read.getDuplicateReadFlag() 
				|| read.getReadFailsVendorQualityCheckFlag() || !read.getReadPairedFlag() || (useBadMate ? false : !read.getProperPairFlag()) 
				|| read.getMappingQuality() < minMapQ){
			return Triple.of(0,0,0);
		}
		
		byte[] refBases = CcInferenceUtils.toUpperCase(modifyRefSeqByCigar(SequenceUtil.makeReferenceFromAlignment(read, false), read.getCigarString()));
		
		byte[] readBasesQ = getClippedReadsBaseQuality(read);
		byte[] readBases = CcInferenceUtils.toUpperCase(getClippedReadsBase(read));
		//if(negStrand){
			//System.err.println(new String(refBases) + "\n" + new String(readBases) + "\t" +read.getCigarString() + "\t" + secondPair);
			//readBases = complementArray(readBases);
		//}
		
		Triple<Integer,Integer, Integer> ctSummary = Triple.of(0,0,0);
		int pos = 0;
		if(!turnOffBisulfiteFilter){
			//System.err.println(new String(refBases) + "\t" + refBases.length + "\n" + new String(readBases) + "\t" + readBases.length + "\n" + new String(SequenceUtil.makeReferenceFromAlignment(read, false)) + "\t" + read.getCigarString());
			pos = CcInferenceUtils.bisulfiteIncompleteReads(readBases, refBases, negStrand, secondPair, false, true, true, 0.1, 0.4, 0, "CH", 1);
			if(pos<0){
				return ctSummary;
			}
		}
		
		ctSummary = readsMethySummaryWithFivePrimeFilter(refBases, readBases, readBasesQ, negStrand, secondPair, minBaseQ, pos);

		return ctSummary; //first is C, second is T
	}

	static public TreeMap<Integer, Integer> readsMethySummaryAtEachLocGeneral(SAMRecord read, String patSearchMeth, int cPos, String patCheckBsConv, int cPosBs, int minMapQ, int minBaseQ, boolean turnOffBisulfiteFilter, boolean useBadMate) throws Exception{
		boolean negStrand = read.getReadNegativeStrandFlag();
		boolean secondPair = read.getReadPairedFlag() && read.getSecondOfPairFlag();
		//filter if it is bad reads
		if(read.getReadUnmappedFlag() || read.getNotPrimaryAlignmentFlag() || read.getDuplicateReadFlag()
				|| read.getReadFailsVendorQualityCheckFlag() || !read.getReadPairedFlag() || (useBadMate ? false : !read.getProperPairFlag())
				|| read.getMappingQuality() < minMapQ){
			return null;
		}

		byte[] refBases = CcInferenceUtils.toUpperCase(modifyRefSeqByCigar(SequenceUtil.makeReferenceFromAlignment(read, false), read.getCigarString()));

		byte[] readBasesQ = getClippedReadsBaseQuality(read);
		byte[] readBases = CcInferenceUtils.toUpperCase(getClippedReadsBase(read));

		if(refBases.length != readBases.length || refBases.length != readBasesQ.length){
			throw new IllegalArgumentException("read base is different from reference base length:" + read.getReadName() + "\t" + read.getCigarString() + "\t" + read.getContig() + "\t" + read.getAlignmentStart()
					+ "\n" + new String(refBases) + "\n" + new String(readBases) + "\n" + new String(readBasesQ));
		}

		TreeMap<Integer, Integer> ctSummary = new TreeMap<Integer, Integer>();
		int pos = 0;
		if(!turnOffBisulfiteFilter){

			pos = CcInferenceUtils.bisulfiteIncompleteReads(readBases, refBases, negStrand, secondPair, false, true, true, 0.1, 0.4, 0, patCheckBsConv, cPosBs);
			if(pos<0){
				return ctSummary;
			}
		}
		byte[] baseSearchMeth = patSearchMeth.getBytes();
		ctSummary = readsMethyAndLocWithFivePrimeFilterGeneral(refBases, readBases, readBasesQ, baseSearchMeth, cPos, negStrand, secondPair, minBaseQ, pos);

		return ctSummary; //first is location, second is methylation status: 0 unmethy, 1 methy.
	}
	
	static public TreeMap<Integer, Integer> readsMethySummaryAtEachLoc(SAMRecord read, int minMapQ, int minBaseQ, boolean turnOffBisulfiteFilter, boolean useBadMate) throws Exception{
		boolean negStrand = read.getReadNegativeStrandFlag();
		boolean secondPair = read.getReadPairedFlag() && read.getSecondOfPairFlag();
		//filter if it is bad reads
		if(read.getReadUnmappedFlag() || read.getNotPrimaryAlignmentFlag() || read.getDuplicateReadFlag() 
				|| read.getReadFailsVendorQualityCheckFlag() || !read.getReadPairedFlag() || (useBadMate ? false : !read.getProperPairFlag()) 
				|| read.getMappingQuality() < minMapQ){
			return null;
		}
		
		byte[] refBases = CcInferenceUtils.toUpperCase(modifyRefSeqByCigar(SequenceUtil.makeReferenceFromAlignment(read, false), read.getCigarString()));
		
		byte[] readBasesQ = getClippedReadsBaseQuality(read);
		byte[] readBases = CcInferenceUtils.toUpperCase(getClippedReadsBase(read));
		
		if(refBases.length != readBases.length || refBases.length != readBasesQ.length){
			throw new IllegalArgumentException("read base is different from reference base length:" + read.getReadName() + "\t" + read.getCigarString() + "\t" + read.getContig() + "\t" + read.getAlignmentStart()
					+ "\n" + new String(refBases) + "\n" + new String(readBases) + "\n" + new String(readBasesQ));
		}
		
		TreeMap<Integer, Integer> ctSummary = new TreeMap<Integer, Integer>();
		int pos = 0;
		if(!turnOffBisulfiteFilter){
			
			pos = CcInferenceUtils.bisulfiteIncompleteReads(readBases, refBases, negStrand, secondPair, false, true, true, 0.1, 0.4, 0, "CH", 1);
			if(pos<0){
				return ctSummary;
			}
		}
		
		ctSummary = readsMethyAndLocWithFivePrimeFilter(refBases, readBases, readBasesQ, negStrand, secondPair, minBaseQ, pos);
		
		return ctSummary; //first is location, second is methylation status: 0 unmethy, 1 methy.
	}
	
	static public TreeMap<Integer, Integer> readsMethySummaryAtEachLoc(SAMRecord read, int minMapQ, int minBaseQ, boolean turnOffBisulfiteFilter) throws Exception{
		boolean negStrand = read.getReadNegativeStrandFlag();
		boolean secondPair = read.getReadPairedFlag() && read.getSecondOfPairFlag();
		//filter if it is bad reads
		if(read.getReadUnmappedFlag() || read.getNotPrimaryAlignmentFlag() || read.getDuplicateReadFlag() 
				|| read.getReadFailsVendorQualityCheckFlag() || !read.getReadPairedFlag()
				|| read.getMappingQuality() < minMapQ){
			return null;
		}
		
		byte[] refBases = CcInferenceUtils.toUpperCase(modifyRefSeqByCigar(SequenceUtil.makeReferenceFromAlignment(read, false), read.getCigarString()));
		
		byte[] readBasesQ = getClippedReadsBaseQuality(read);
		byte[] readBases = CcInferenceUtils.toUpperCase(getClippedReadsBase(read));
		
		if(refBases.length != readBases.length || refBases.length != readBasesQ.length){
			throw new IllegalArgumentException("read base is different from reference base length:" + read.getReadName() + "\t" + read.getCigarString() + "\t" + read.getContig() + "\t" + read.getAlignmentStart()
					+ "\n" + new String(refBases) + "\n" + new String(readBases) + "\n" + new String(readBasesQ));
		}
		
		TreeMap<Integer, Integer> ctSummary = new TreeMap<Integer, Integer>();
		int pos = 0;
		if(!turnOffBisulfiteFilter){
			
			pos = CcInferenceUtils.bisulfiteIncompleteReads(readBases, refBases, negStrand, secondPair, false, true, true, 0.1, 0.4, 0, "CH", 1);
			if(pos<0){
				return ctSummary;
			}
		}
		
		ctSummary = readsMethyAndLocWithFivePrimeFilter(refBases, readBases, readBasesQ, negStrand, secondPair, minBaseQ, pos);
		
		return ctSummary; //first is location, second is methylation status: 0 unmethy, 1 methy.
	}

	
	static public Triple<Integer,Integer, Integer> readsMethySummary(SAMRecord read, int minMapQ, int minBaseQ){
		boolean negStrand = read.getReadNegativeStrandFlag();
		boolean secondPair = read.getReadPairedFlag() && read.getSecondOfPairFlag();

		byte[] refBases = CcInferenceUtils.toUpperCase(modifyRefSeqByCigar(SequenceUtil.makeReferenceFromAlignment(read, false), read.getCigarString()));
		
		byte[] readBasesQ = getClippedReadsBaseQuality(read);
		byte[] readBases = CcInferenceUtils.toUpperCase(getClippedReadsBase(read));
		if(negStrand){
			readBases = complementArray(readBases);
		}
		Triple<Integer,Integer, Integer> ctSummary = readsMethySummary(refBases, readBases, readBasesQ, negStrand, secondPair, minBaseQ);
		
		
		
		return ctSummary; //first is methylation level, second is Cpg number
	}
	
	static public Triple<Integer,Integer, Integer> readsMethySummary(byte[] refBases, byte[] bases, byte[] basesQ, boolean negStrand, boolean secondPair, int minBaseQ){
		
		int numC = 0;
		int numT = 0;
		int numCg = 0;
		if((!secondPair && !negStrand) || (secondPair && negStrand)){
			for(int i = 0; i < refBases.length -1; i ++){
				if(SequenceUtil.basesEqual(refBases[i], SequenceUtil.C) && SequenceUtil.basesEqual(refBases[i+1], SequenceUtil.G)){
					if(SequenceUtil.basesEqual(bases[i], SequenceUtil.C) && basesQ[i] >= minBaseQ){
						numC++;
						numCg++;
					}else if(SequenceUtil.basesEqual(bases[i], SequenceUtil.T) && basesQ[i] >= minBaseQ){
						numT++;
						numCg++;
					}
				}
			}
		}else{
			for(int i = 1; i < refBases.length; i ++){
				if(SequenceUtil.basesEqual(refBases[i], SequenceUtil.G) && SequenceUtil.basesEqual(refBases[i-1], SequenceUtil.C)){
					if(SequenceUtil.basesEqual(bases[i], SequenceUtil.G) && basesQ[i] >= minBaseQ){
						numC++;
						numCg++;
					}else if(SequenceUtil.basesEqual(bases[i], SequenceUtil.A) && basesQ[i] >= minBaseQ){
						numT++;
						numCg++;
					}
				}
			}
		}
		//Double methy = Double.NaN;
		//if(numC + numT > 0){
		//	methy = (double)numC/(double)(numC + numT);
		//}
		 
		return Triple.of(numC, numT,numCg);  //first is numC, second is numT, third is Cpg number
	}

	//cPos: for GCH, it is 2. for CG, it is 1.
	//searchBases allow IUAUC match: H: A,C,T
	static public Triple<Integer,Integer, Integer> readsMethySummaryWithFivePrimeFilterGeneral(byte[] refBases, byte[] bases, byte[] basesQ, byte[] searchBases, int cPos, boolean negStrand, boolean secondPair, int minBaseQ, int bisulfiteConvertStart){
		if(bisulfiteConvertStart<0){
			bisulfiteConvertStart=0;
		}
		if((!secondPair && !negStrand) || (secondPair && negStrand)) {

		}else {
			//reverse completment refBases, bases, reverse basesQ
			bases = BaseUtils.simpleReverseComplement(bases);
			refBases = BaseUtils.simpleReverseComplement(refBases);
			basesQ = BaseUtilsMore.simpleReverse(basesQ);
		}
		int numC = 0;
		int numT = 0;
		int numCg = 0;
		for(int i = bisulfiteConvertStart; i < refBases.length - (searchBases.length-1); i++){
			int j = i;
			boolean found = true;
			for(int index = 0; index < searchBases.length; j++,index++) {
				if (!BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(searchBases[index], refBases[j])) {
					found = false;
					break;
				}
			}
			if(found){
				int currentCytPos = i+cPos-1;
				//System.err.println(new String(refBases) + "\t" + bisulfiteConvertStart + "\t" + i + "\t" + new String(bases) + "\t" + j + "\t" + new String(basesQ) + "\t" + negStrand + "\t" + secondPair);
				if(SequenceUtil.basesEqual(bases[currentCytPos], SequenceUtil.C) && basesQ[currentCytPos] >= minBaseQ){
					numC++;
					numCg++;
				}else if(SequenceUtil.basesEqual(bases[currentCytPos], SequenceUtil.T) && basesQ[currentCytPos] >= minBaseQ){
					numT++;
					numCg++;
				}
			}
		}
		return Triple.of(numC, numT,numCg); //first is numC, second is numT, third is Cpg number
	}


	static public Triple<Integer,Integer, Integer> readsMethySummaryWithFivePrimeFilter(byte[] refBases, byte[] bases, byte[] basesQ, boolean negStrand, boolean secondPair, int minBaseQ, int bisulfiteConvertStart){
		if(bisulfiteConvertStart<0){
			bisulfiteConvertStart=0;
		}
		int numC = 0;
		int numT = 0;
		int numCg = 0;
		if((!secondPair && !negStrand) || (secondPair && negStrand)){
			for(int i = bisulfiteConvertStart; i < refBases.length -1; i++){
				if(SequenceUtil.basesEqual(refBases[i], SequenceUtil.C) && SequenceUtil.basesEqual(refBases[i+1], SequenceUtil.G)){
					if(SequenceUtil.basesEqual(bases[i], SequenceUtil.C) && basesQ[i] >= minBaseQ){
						numC++;
						numCg++;
					}else if(SequenceUtil.basesEqual(bases[i], SequenceUtil.T) && basesQ[i] >= minBaseQ){
						numT++;
						numCg++;
					}
				}
			}
		}else{
			for(int i = 1; i < refBases.length-bisulfiteConvertStart; i++){
				//System.err.println(refBases.length + "\t" + bisulfiteConvertStart + "\t" + i + "\t" + bases.length);
				if(SequenceUtil.basesEqual(refBases[i], SequenceUtil.G) && SequenceUtil.basesEqual(refBases[i-1], SequenceUtil.C)){
					if(SequenceUtil.basesEqual(bases[i], SequenceUtil.G) && basesQ[i] >= minBaseQ){
						numC++;
						numCg++;
					}else if(SequenceUtil.basesEqual(bases[i], SequenceUtil.A) && basesQ[i] >= minBaseQ){
						numT++;
						numCg++;
					}
				}
			}
		}
		//Double methy = Double.NaN;
		//if(numC + numT > 0){
		//	methy = (double)numC/(double)(numC + numT);
		//}
		 
		return Triple.of(numC, numT,numCg); //first is numC, second is numT, third is Cpg number
	}

	static public TreeMap<Integer, Integer> readsMethyAndLocWithFivePrimeFilterGeneral(byte[] refBases, byte[] bases, byte[] basesQ, byte[] searchBases, int cPos, boolean negStrand, boolean secondPair, int minBaseQ, int bisulfiteConvertStart){
		if(bisulfiteConvertStart<0){
			bisulfiteConvertStart=0;
		}
		if((!secondPair && !negStrand) || (secondPair && negStrand)) {

		}else {
			//reverse completment refBases, bases, reverse basesQ
			bases = BaseUtils.simpleReverseComplement(bases);
			refBases = BaseUtils.simpleReverseComplement(refBases);
			basesQ = BaseUtilsMore.simpleReverse(basesQ);
		}
		TreeMap<Integer, Integer> cpgMethy = new TreeMap<Integer, Integer>(); //key: relative location from beginning of fragment, value: 0:unmethy, 1: methy

		for(int i = bisulfiteConvertStart; i < refBases.length - (searchBases.length-1); i++){
			int j = i;
			boolean found = true;
			for(int index = 0; index < searchBases.length; j++,index++) {
				if (!BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(searchBases[index], refBases[j])) {
					found = false;
					break;
				}
			}
			if(found){
				int currentCytPos = i+cPos-1;
				//System.err.println(new String(refBases) + "\t" + bisulfiteConvertStart + "\t" + i + "\t" + new String(bases) + "\t" + j + "\t" + new String(basesQ) + "\t" + negStrand + "\t" + secondPair);
				if(SequenceUtil.basesEqual(bases[currentCytPos], SequenceUtil.C) && basesQ[currentCytPos] >= minBaseQ){
					cpgMethy.put(currentCytPos, 1);
				}else if(SequenceUtil.basesEqual(bases[currentCytPos], SequenceUtil.T) && basesQ[currentCytPos] >= minBaseQ){
					cpgMethy.put(currentCytPos, 0);
				}
			}
		}

		return cpgMethy; //key: relative location from beginning of fragment, value: 0:unmethy, 1: methy
	}
	
	static public TreeMap<Integer, Integer> readsMethyAndLocWithFivePrimeFilter(byte[] refBases, byte[] bases, byte[] basesQ, 
			boolean negStrand, boolean secondPair, int minBaseQ, int bisulfiteConvertStart){
		//ArrayList<Interval> cpgLoc = new ArrayList<Interval>();
		//ArrayList<Integer> cpgMethy = new ArrayList<Integer>();
		//new Interval(chr, start, end);
		TreeMap<Integer, Integer> cpgMethy = new TreeMap<Integer, Integer>(); //key: relative location from beginning of fragment, value: 0:unmethy, 1: methy

		
		if((!secondPair && !negStrand) || (secondPair && negStrand)){
			for(int i = bisulfiteConvertStart; i < refBases.length -1; i++){
				if(SequenceUtil.basesEqual(refBases[i], SequenceUtil.C) && SequenceUtil.basesEqual(refBases[i+1], SequenceUtil.G)){
					if(SequenceUtil.basesEqual(bases[i], SequenceUtil.C) && basesQ[i] >= minBaseQ){
						cpgMethy.put(i, 1);
					}else if(SequenceUtil.basesEqual(bases[i], SequenceUtil.T) && basesQ[i] >= minBaseQ){
						cpgMethy.put(i, 0);
					}
				}
			}
		}else{
			for(int i = 1; i < refBases.length-bisulfiteConvertStart; i++){
				if(SequenceUtil.basesEqual(refBases[i], SequenceUtil.G) && SequenceUtil.basesEqual(refBases[i-1], SequenceUtil.C)){
					if(SequenceUtil.basesEqual(bases[i], SequenceUtil.G) && basesQ[i] >= minBaseQ){
						cpgMethy.put(i-1, 1);
					}else if(SequenceUtil.basesEqual(bases[i], SequenceUtil.A) && basesQ[i] >= minBaseQ){
						cpgMethy.put(i-1, 0);
					}
				}
				
			}
		}
		
		
		return cpgMethy; //key: relative location from beginning of fragment, value: 0:unmethy, 1: methy
	}
	
	
	static public boolean baseMatchInBisulfiteSpace(byte refBase, byte base, boolean negStrand, boolean secondPair) {
		if (!negStrand) {
			if (secondPair) {
				if (SequenceUtil.basesEqual(refBase, SequenceUtil.G) && SequenceUtil.basesEqual(base, SequenceUtil.A)) {
					return true;
				} else {
					return SequenceUtil.basesEqual(refBase, base);
				}
			} else {
				if (SequenceUtil.basesEqual(refBase, SequenceUtil.C) && SequenceUtil.basesEqual(base, SequenceUtil.T)) {
					return true;
				} else {
					return SequenceUtil.basesEqual(refBase, base);
				}
			}

		} else {
			if (secondPair) {
				if (SequenceUtil.basesEqual(refBase, SequenceUtil.C) && SequenceUtil.basesEqual(base, SequenceUtil.T)) {
					return true;
				} else {
					return SequenceUtil.basesEqual(refBase, base);
				}
			} else {
				if (SequenceUtil.basesEqual(refBase, SequenceUtil.G) && SequenceUtil.basesEqual(base, SequenceUtil.A)) {
					return true;
				} else {
					return SequenceUtil.basesEqual(refBase, base);
				}
			}

		}
	}
	
	static public int findPosition(ArrayList<Double> permutatedProbList, double observedProb, int lower, int upper){
		if(lower == 0 && observedProb <= permutatedProbList.get(0)){
			return lower;
		}else if(upper == permutatedProbList.size()-1 && observedProb > permutatedProbList.get(permutatedProbList.size()-1)){
			return upper;
		}else if((upper-lower)==1 && (observedProb>=permutatedProbList.get(lower) && observedProb<=permutatedProbList.get(upper))){
			return upper;
		}else if(observedProb>=permutatedProbList.get(lower) && observedProb<=permutatedProbList.get(upper)){
			return findPosition(permutatedProbList, observedProb, lower+(upper-lower)/2, upper);
		}else if(observedProb < permutatedProbList.get(lower)){
			return findPosition(permutatedProbList, observedProb, Math.max(lower-(upper-lower),0), lower);
		}else{
			return findPosition(permutatedProbList, observedProb, upper, Math.max(upper+(upper-lower),permutatedProbList.size()-1));
		}
		
	}
	
}
