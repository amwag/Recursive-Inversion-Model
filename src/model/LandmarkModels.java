package model;

import java.io.IOException;
import java.util.Arrays;
import java.util.Random;

public class LandmarkModels {
	int nitems;
	int nlandmarks;
	int ntotal;
	double[] thetas;
	int[] perm;
	int[][] GIDs; //N by nitems+nlandmarks
	int[][] Q; //nitems+nlandmarks square
	int[][][] Ms; //N by nlandmarks(-1) by nitems (Which items after nth landmark?)
	
	double epsilon=0.0001;
	
	
	public static void main(String[] args) throws IOException{
		final long startTime = System.currentTimeMillis();
		
		int[] perm= {0,1,8,2,3,9,4,5,10,6,7};
		int m=3;
		double[] thetas= {0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,0.99};
		int N=10;
		
		Random rand=new Random(101);
		
		for(int[] q1:sampleLGMMPerms(N,m,thetas,perm,rand)) {
			System.out.println(Arrays.toString(q1));
		}
	}
	
	static int[][] sampleLGMMPerms(int N, int m, double[] thetas, int[] perm){
		//RANDOMIZATION IS FIXED - REMOVE LATER
		
		Random rand=new Random();
		return(sampleLGMMPerms(N,m,thetas,perm,rand));
	}
	
	static int[][] sampleLGMMPerms(int N, int m, double[] thetas, int[] perm, Random rand){
		/*
		 * @param N number of perms to be sampled
		 * @param m number of landmarks in the perm
		 * @param thetas Parameters of the model
		 * @param perm Perm to sample from, with [0...n-1] as items and [n...n+m-1] as landmarks
		 */
		
		int[][] svals=new int[N][perm.length];
		int[] sind=new int[N];
		for(int i=0;i<perm.length;i++) {
			sind=sampleDiscreteExp(thetas[i],perm.length-i,N,rand);
			for(int k=0;k<N;k++)
				svals[k][i]=sind[k];
		}
		
		int[][] perms=new int[N][perm.length];
		
		int[] tmpperm=new int[perm.length];
		int lmstart=perm.length-m;
		int lmitem;
		for(int k=0;k<N;k++) {
			System.arraycopy(perm, 0, tmpperm, 0, perm.length);
			lmitem=lmstart;
			for(int i=0;i<perm.length-1;i++) {
				if(tmpperm[svals[k][i]]<lmstart) {
					//If normal item, add to perm, remove from tempperm
					perms[k][i]=tmpperm[svals[k][i]];
					for(int j=svals[k][i];j<perm.length-i-1;j++)
						tmpperm[j]=tmpperm[j+1];
				}
				else {
					//If landmark, add first lm to perm, remove first lm from tempperm
					perms[k][i]=lmitem++;
					int j=0;
					while(tmpperm[j]<lmstart)
						j++;
					for(;j<perm.length-i-1;j++)
						tmpperm[j]=tmpperm[j+1];
				}
			}
			perms[k][perm.length-1]=tmpperm[0];
			System.out.println("s: "+Arrays.toString(svals[k]));
			System.out.println("p: "+Arrays.toString(perm));
			System.out.println("f: "+Arrays.toString(perms[k]));
		}
		
		/*
		//We reindex everything to the correct index of perm before filling the perm matrix
		for(int k=0;k<N;k++) {
			System.out.println("a: "+Arrays.toString(svals[k]));
			for(int i=0;i<perm.length-1;i++) {
				int mod=0;
				for(int j=i+1;j<perm.length;j++) {
					if(svals[k][i]<=svals[k][j]+mod)
						svals[k][j]++;
					else mod++;
				}
			}
			System.out.println("b: "+Arrays.toString(svals[k]));
			System.out.println("p: "+Arrays.toString(perm));
			
			//These are now indexed to the correct value, but landmarks are not in order - need to 
			
			int lm=perm.length-m;
			for(int i=0;i<perm.length;i++) {
				if(perm[svals[k][i]]<perm.length-m) {
					perms[k][i]=perm[svals[k][i]];
				}
				else {
					perms[k][i]=lm++;
				}
			}
			System.out.println("q: "+Arrays.toString(perms[k]));
		}
		*/ //OLD ATTEMPT DID NOT WORK AS INTENDED
		
		return perms;
	}
	

	
	static int[] sampleDiscreteExp(double theta, int n, int N, Random rand) {
		//Creates the CDF for a discrete exponential, then samples by calling sampleDiscreteExpCDF
		//Samples N exp(theta,n) values
		double[] CDF=new double[n];
		CDF[0]=1;
		for(int i=1;i<n;i++) {
			CDF[i]=CDF[i-1]+Math.pow(theta, i);
		}
		for(int i=0;i<n;i++)
			CDF[i]=CDF[i]/CDF[n-1];
		
		return(sampleDiscreteExpCDF(CDF,N,rand));
	}
	
	static int[] sampleDiscreteExpCDF(double[] CDF, int N, Random rand) {
		//Samples N exp(theta,n) values
		for(int i=0;i<CDF.length;i++)
			System.out.println(CDF[i]);
		
		int[] sample=new int[N];
		int k=0;
		double randval;
		for(int i=0;i<N;i++) {
			randval=rand.nextDouble();
			while(CDF[k]<randval)
				k++;
			sample[i]=k;
			k=0;
		}
		
		return(sample);
	}
}
