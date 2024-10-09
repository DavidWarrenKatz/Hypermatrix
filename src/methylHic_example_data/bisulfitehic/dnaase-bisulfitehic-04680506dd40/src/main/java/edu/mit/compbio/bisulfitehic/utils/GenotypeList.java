/**
 * GenotypeList.java
 * Nov 17, 2015
 * 3:10:44 PM
 * yaping    lyping1986@gmail.com
 */
package main.java.edu.mit.compbio.bisulfitehic.utils;

import java.io.Serializable;
import java.util.HashMap;
import java.util.TreeMap;

import org.apache.commons.math3.util.Pair;
import org.apache.commons.lang3.tuple.Triple;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.GenomeLocSortedSet;

/**
 *
 */
public class GenotypeList implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = -1073731400420098203L;
	/**
	 * 
	 */
	public GenomeLocSortedSet genomeLocsCol = null;
	public TreeMap<GenomeLoc, Pair<Byte, Byte>> genomeLocsAlleleMap = null;
	public TreeMap<GenomeLoc, Triple<Byte, Byte, Byte>> genomeLocsAlleleMapWithRef = null;
	public GenomeLocParser glp = null;
	public GenotypeList(GenomeLocSortedSet genomeLocsCol, TreeMap<GenomeLoc, Pair<Byte, Byte>> genomeLocsAlleleMap, GenomeLocParser glp){
		this.genomeLocsCol = genomeLocsCol;
		this.genomeLocsAlleleMap = genomeLocsAlleleMap;
		this.glp = glp;
	}

	public GenotypeList(GenomeLocSortedSet genomeLocsCol, GenomeLocParser glp, TreeMap<GenomeLoc, Triple<Byte, Byte, Byte>> genomeLocsAlleleMapWithRef){
		this.genomeLocsCol = genomeLocsCol;
		this.genomeLocsAlleleMapWithRef = genomeLocsAlleleMapWithRef;
		this.glp = glp;
	}

}
