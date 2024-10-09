/**
 * 
 */
package main.java.edu.mit.compbio.bisulfitehic.utils;

import java.math.BigDecimal;
import java.util.Random;

import org.apache.log4j.BasicConfigurator;

/**
 * @author yaping
 * interpreted from scala code at https://deaktator.github.io/2016/02/12/poisson-binomial-revisited/
 */
public class PoissonBinomial2 {
	//private BigDecimal zero;
	//private BigDecimal one;
	double[] pmfList;
	double[] cdfList;
	public PoissonBinomial2(double[] cpg1) {
		//zero = new BigDecimal(0);
		//one = new BigDecimal(1);
		pmfList = pmf(cpg1);
		cdf();
	}
	
	/*
			public static void main(String[] args)
					throws Exception {
		double[] test = new double[10];
		Random r = new Random();
		for(int i = 0; i < test.length; i++) {
			test[i] = r.nextdouble();
			System.err.print(test[i] + "\t");
		}
		System.err.println();
		PoissonBinomial2 pb = new PoissonBinomial2(test);
		double[] pmf = pb.pmf(test);
		System.err.println(pb.get_pmf(10) + "\t" + pb.get_cdf(10));
		}
	*/
	
	public double get_pmf(int numSucess) {
		if(numSucess>=pmfList.length) {
			return Double.NaN;
		}
		return pmfList[numSucess];
	}
	
	public double get_cdf(int numSucess) {
		if(numSucess>=cdfList.length) {
			return Double.NaN;
		}
		return cdfList[numSucess];
	}
				
	public void cdf() {
		cdfList = new double[pmfList.length];
		double cum = 0;
		for(int i = 0; i < cdfList.length; i++) {
			cum += pmfList[i];
			cdfList[i] = cum;
		}
		
	}
	
	public double[] pmf(double[] pr) {
		return pmf(pr, Integer.MAX_VALUE, 2);
	}
	
	private double[] pmf(double[] pr, int maxN, int maxCumPr) {
		double[] w = new double[pr.length];
		for(int i = 0; i < pr.length; i++) {
			w[i] = pr[i]/(1-pr[i]);
		}
		int n = w.length;
		int mN = Math.min(maxN, n);
		double z = 1;
		for( double s : w) {
			z = z / (s+1);
		}
		double[] r = new double[n+1];
		for(int i = 0; i< r.length; i++) {
			r[i] = 1.;
		}
		r[n] = z;
		int i = 1;
		int j = 0;
		int k = 0;
		int m = 0;
		double s = 0;
		double cumPr = r[n];
		while(cumPr < maxCumPr && i <= mN) {
		      s = 0;  j = 0;  m = n - i;  k = i - 1;
		      
		      while (j <= m) {
		        s += r[j] * w[k + j];
		        r[j] = s;
		        j += 1;
		      }
		      
		      r[j - 1] *= z;
		      cumPr += r[j - 1];
		      
		      i += 1;
		  }
		
		return finalizeR(r, i, n);
				
	}
	
	private double[] finalizeR(double[] r, int i, int n){
		    if (i <= n) {
		    		double[] smallerR = new double[i];
		      System.arraycopy(r, n - i + 1, smallerR, 0, i);
		       reverse(smallerR);
		      return smallerR;
		    }
		    else{
		    		reverse(r);
		    		return r;
		    }
	}
	
	private void reverse(double[] a){
		    int n = a.length / 2;
		    int i = 0;
		    double t = a[0];
		    int j = 0;
		    while (i <= n) {
		      j = a.length - i - 1;
		      t = a[i];
		      a[i] = a[j];
		      a[j] = t;
		      i += 1;
		    }
		  }
}
