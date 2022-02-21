package model;

import java.io.IOException;
import java.util.Arrays;
import java.util.Random;

public class LandmarkGMM {
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
		
		//String filename = "C:\\Users\\Anna\\Dropbox\\AdaptiveRim\\data\\ISSPMilestone\\ISSP_LM_f_12345.txt";
		
		//int[][][] GID = new int[1][6][6];
		
		
		
		//GID[0][0]= {1,2,3,4,5,6,7};
		//GID[0] = AbstractModel.loadGroupIDs(filename);
		
		//CRIM generator = new CRIM(GID, 0, 0);
		
		// There are 4 landmarks, so 1 for each landmark and 'end'
		// Marks if item i is left of landmark k for perm n at [n][k-1][i]
		/*
		int[][][] Mleft=new int[5][generator.N_pi0[0]][20];
		for(int n=0;n<generator.N_pi0[0];n++) {
			for(int k=0;k<4;k++) {
				for(int i=0;i<20;i++) {
					Mleft[k][n][i]=GID[0][n][i]<GID[0][n][20+k] ? 0:1;
				}
			}
			int k=4;
			for(int i=0;i<20;i++) {
				Mleft[k][n][i]=GID[0][n][i]>GID[0][n][19+k] ? 0:1;
			}
		}
		*/
		//int k=20;
		//int[] inds= {1,2};
		//int[] M=getMs(k,inds,20,Mleft);
		//double vl=getVl(k,inds,generator.Qobs_pi0[0]);
		
		//int l=5;
		//int[] inds2={10,12,14};
		//double[][] LR=getLRCount(l,inds2,GID[0]);
		//double pVbar=getVl(l,inds2,generator.Qobs_pi0[0]);
		
		//System.out.println("Hello");
		//System.out.println(getThetaStandard(LR, pVbar, generator));
		
		//System.out.println("Hello again");
		//System.out.println(getThetaMilestone(M, vl));
		
		
		
		// 7, 9, and 2 values are worrisome. 
		
		
		
		int[] items={0,1,2,3,4,5,6};
		int[] landmarks={7,8};
		
		int[] myi={0};
		double[] myt={0};
		myLL head=new myLL(myi,myt,100);
		
		int maxsize=2;
		
		double[][] Z=new double[2][2];
		
		int[] gmmitems = {0,1,2,7,3,4,8,5,6};
		double[] thetavec = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};
		int nitems=7;
		int nlandmarks=2;
		int samplesize=10000;
		//int[][] tmpIntArray=sampleLandmarkVs(gmmitems, thetavec, nitems, samplesize);
		

		
		/*
		int[][][] GIDs=new int[1][samplesize][nlandmarks+nitems];
		int[] perm={0,1,2,7,3,4,8,5,6};
		
		for(int p=0;p<samplesize;p++) {
			int[] test=tmpIntArray[p];
			int[] inds=new int[9];
			for(int i=perm.length-2;i>=0;i--) {
				for(int k=i+1;k<perm.length;k++)
					if(inds[k]<test[i])
						inds[i]++;
					else
						inds[k]++;
			}
			for(int i=0;i<perm.length;i++)
				GIDs[0][p][perm[i]]=inds[i];
		}
		
		int[] test= {-1,-1,-1,-1,-1,-1,-1,5,6};
		
		System.out.println("Hello friends");
		System.out.println("Perm  - "+Arrays.toString(perm));
		System.out.println(" Vs   - "+Arrays.toString(tmpIntArray[0]));
		
		int[] Vs=tmpIntArray[0];
		if(Vs[7]==1) {
			test[7]=6;
			test[8]=5;
		}
		for(int i=6;i>=0;i--) {
			int tmpi=i;
			int v=Vs[i];
			for(v=Vs[i];v>0;v--) {
				test[tmpi]=test[tmpi+1];
				tmpi++;
			}
			test[tmpi]=perm[i];
		}
		System.out.println("Operm - "+Arrays.toString(test));
		
		System.out.println(" GIDs - "+Arrays.toString(GIDs[0][0]));
		System.out.println("Save me from myself");
		
		*/

		int[] modelperm={0,1,7,2,3,4,8,5,6};
		double[] modeltheta={1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8};
		LandmarkGMM test=new LandmarkGMM(7,2,modelperm,modeltheta);
		
		double[][] vcdf=test.getModelVCDF();
		for(int n=0;n<vcdf.length;n++)
			System.out.println(Arrays.toString(vcdf[n]));
		
		int[][] perms=test.samplePerms(10000);
		int[][][] GIDs=new int[1][perms.length][perms[0].length];
		
		for(int n=0;n<perms.length;n++) {
			for(int k=0;k<perms[0].length;k++) {
				GIDs[0][n][perms[n][k]]=k;
			}
		}
		
		int[] ritems1= {0,2,4,6,8};
		System.out.println(Arrays.toString(ritems1));
		
		System.out.println("Perm: "+Arrays.toString(perms[0]));
		test.setGIDs(perms);
		System.out.println("GIDs: "+Arrays.toString(test.GIDs[0]));	
		test.setMs();
		System.out.println("Ms: "+Arrays.toString(test.Ms[0][0]));

		int[] M1=test.getMVec(7, ritems1);
		System.out.println("M: "+M1[0]);
		
		System.out.println("Perm: "+Arrays.toString(perms[50]));
		test.setGIDs(perms);
		System.out.println("GIDs: "+Arrays.toString(test.GIDs[50]));	
		test.setMs();
		System.out.println("Ms: "+Arrays.toString(test.Ms[50][0]));
		M1=test.getMVec(7, ritems1);
		System.out.println("M: "+M1[50]);
		
		for(int[] q:test.Q) {
			System.out.println("Q: "+Arrays.toString(q));
		}
		
		System.out.println("Node for 7");
		
		int[] ritems2= {2,3,4,8,5,6};
		System.out.println(test.getNegLogLike(7,ritems2,0)+" "+test.getNegdLogLike(7,ritems2,0.28));
		System.out.println(test.getNegLogLike(7,ritems2,0.3)+" "+test.getNegdLogLike(7,ritems2,0.3));
		System.out.println(test.getNegLogLike(7,ritems2,0.35)+" "+test.getNegdLogLike(7,ritems2,0.32));

		System.out.println("Node for 3");
		
		int[] ritems3= {4,8,5,6};
		System.out.println(test.getNegLogLike(3,ritems3,0)+" "+test.getNegdLogLike(3,ritems3,0.48));
		System.out.println(test.getNegLogLike(3,ritems3,0.5)+" "+test.getNegdLogLike(3,ritems3,0.5));
		System.out.println(test.getNegLogLike(3,ritems3,0.55)+" "+test.getNegdLogLike(3,ritems3,0.52));
		
		int[] perm={0,1,7,2,3,4,8,5,6};
		test.fitPerm(perm);
		
		int[] perm1={0,1,2,3,4,7,8,5,6};
		test.fitPerm(perm1);
		
		int[] perm2={1,0,7,2,3,4,8,5,6};
		test.fitPerm(perm2);
		
		int[] perm3={7,8,0,1,2,3,4,5,6};
		test.fitPerm(perm3);
		
		int[] perm4={0,1,2,3,4,5,6,7,8};
		test.fitPerm(perm4);
		
		/*
		CRIM generator = new CRIM(GIDs, 0, 0);
		double[][] help=new double[2][2];
		generator.set_dLogZ(help, 0.8);
		System.out.println(Arrays.toString(help[1]));
		*/
		/*
		int[][][] Mleft=new int[nlandmarks+1][generator.N_pi0[0]][nitems];
		for(int n=0;n<generator.N_pi0[0];n++) {
			for(int k=0;k<nlandmarks;k++) {
				for(int i=0;i<nitems;i++) {
					Mleft[k][n][i]=GIDs[0][n][i]<GIDs[0][n][nitems+k] ? 1:0;
				}
			}
			int k=nlandmarks;
			for(int i=0;i<nitems;i++) {
				Mleft[k][n][i]=GIDs[0][n][i]>GIDs[0][n][nitems-1+k] ? 1:0;
			}
		}
		System.out.println(Arrays.toString(GIDs[0][0]));
		System.out.println(Arrays.toString(Mleft[0][0]));
		System.out.println(Arrays.toString(Mleft[1][0]));
		System.out.println(Arrays.toString(Mleft[2][0]));
		
		
		for(double[] q:generator.Qobs_pi0[0])
			System.out.println(Arrays.toString(q));
		
		for(int i:items) {
			for(int j:items) {
				if(i!=j) {
					int[] ritems={j};
					int[] fullitems= {i,j};
					double[][] LRCount=getLRCount(i,ritems,GIDs[0]);
					double pVbar=getVl(i,ritems,generator.Qobs_pi0[0]);
					double[] thetas= {getThetaStandard(LRCount, pVbar, generator)};
					if(thetas[0]>0) {
						double logLike=generator.logLikelihood(Z, LRCount, pVbar, thetas[0], samplesize);
						myLL tmp=new myLL(fullitems,thetas,logLike);
						tmp.appendTo(head);
					}
				}
			}
			
			int[] ritems={landmarks[landmarks.length-1]};
			int[] fullitems= {i,landmarks[landmarks.length-1]};
			double[][] LRCount=getLRCount(i,ritems,GIDs[0]);
			double pVbar=getVl(i,ritems,generator.Qobs_pi0[0]);
			double[] thetas= {getThetaStandard(LRCount, pVbar, generator)};
			if(thetas[0]>0) {
				double logLike=generator.logLikelihood(Z, LRCount, pVbar, thetas[0], samplesize);
				myLL tmp=new myLL(fullitems,thetas,logLike);
				System.out.println(tmp+" "+logLike);
				tmp.appendTo(head);
			}
			
			ritems[0]=i;
			int[] M=getMs(landmarks[landmarks.length-1],ritems,nitems,Mleft);
			double vl=getVl(landmarks[landmarks.length-1],ritems,generator.Qobs_pi0[0]);
			double[] thetas2= {getThetaMilestone(M, vl)};
			int[] fullitems2= {landmarks[landmarks.length-1],i};
			if(thetas2[0]>0) {
				double logLike=logPMilestone(M,vl,thetas2[0]);
				myLL tmp=new myLL(fullitems2,thetas2,logLike);
				System.out.println(tmp+" "+logLike);
				tmp.appendTo(head);
			}
		}
		
		System.out.println("Initialized");

		myLL newll=head;
		while(newll.hasNext()) {
			newll=newll.next;
			System.out.println(newll);
		}	
		
		int stop=9;//nitems+nlandmarks;
		
		while(head.next.items.length<stop) {
			myLL node=head.next;
			head.next=node.next;
			head.next.prev=head;
			
			double[][] Z2=new double[2][node.items.length+1];
			
			for(int i:items) {
				if(!node.contains(i)) {
						int[] fullitems=new int[node.items.length+1];
						fullitems[0]=i;
						System.arraycopy(node.items, 0, fullitems, 1, node.items.length);
						double[][] LRCount=getLRCount(i,node.items,GIDs[0]);
						double pVbar=getVl(i,node.items,generator.Qobs_pi0[0]);
						double[] thetas=new double[node.items.length];
						thetas[0]=getThetaStandard(LRCount, pVbar, generator);
						System.arraycopy(node.thetas, 0, thetas, 1, node.items.length-1);
						if(thetas[0]>=0) {
							double logLike=generator.logLikelihood(Z2, LRCount, pVbar, thetas[0], samplesize);
							logLike+=node.loglikelihood;
							myLL tmp=new myLL(fullitems,thetas,logLike);
							tmp.appendTo(head);
						}
				}
			}
				
			int k=nlandmarks;
			while(k>0 && node.contains(landmarks[k-1]))
					k--;
					
			if(k>0) {
				int[] M=getMs(landmarks[k-1],node.items,nitems,Mleft);
				double vl=getVl(landmarks[k-1],node.items,generator.Qobs_pi0[0]);
				double[] thetas2=new double[node.items.length];		
				thetas2[0]=getThetaMilestone(M, vl);
				System.arraycopy(node.thetas, 0, thetas2, 1, node.items.length-1);
				int[] fullitems2=new int[node.items.length+1];
				fullitems2[0]=landmarks[k-1];
				System.arraycopy(node.items, 0, fullitems2, 1, node.items.length);
				
				if(thetas2[0]>=0) {
					double logLike=logPMilestone(M,vl,thetas2[0]);
					logLike+=node.loglikelihood;
					myLL tmp=new myLL(fullitems2,thetas2,logLike);
					tmp.appendTo(head);
				}
			}
			if(maxsize==node.items.length) {
				maxsize++;
				System.out.println("Branched node with "+(maxsize-1)+" items after "+((System.currentTimeMillis()-startTime)/1000)+" seconds of runtime");
			}
		}
			
		System.out.println("Finished?");
		
		int i=0;
		myLL newll2=head;
		while(newll2.hasNext() && i++<10) {
			newll2=newll2.next;
			System.out.println(newll2);
		}
		*/
	}
	
	LandmarkGMM(int iitems, int ilandmarks, int[] iperm, double[] ithetas){
		if(iitems+ilandmarks!=iperm.length) {
			System.out.println("Incorrect perm length for number of items");
			return;
		}
		if(iperm.length!=ithetas.length+1) {
			System.out.println("Error: incorrect number of thetas for perm length");
			return;
		}
				
		nitems=iitems;
		nlandmarks=ilandmarks;
		ntotal=nitems+nlandmarks;
		perm=new int[iperm.length];
		System.arraycopy(iperm,0,perm,0,iperm.length);
		perm=iperm;
		thetas=new double[ithetas.length];
		System.arraycopy(ithetas,0,thetas,0,ithetas.length);
		thetas=ithetas;
		System.out.println("Initialized Landmark GMM:");
		String outstr="(";
		for(int i=0;i<iperm.length-1;i++) {
			if(perm[i]<nitems)
				outstr+=perm[i]+" ("+thetas[i]+"); ";
			else
				outstr+="|"+perm[i]+"| ("+thetas[i]+"); ";
		}
		if(perm[perm.length-1]<nitems)
			outstr+=perm[perm.length-1]+")";
		else
			outstr+="|"+perm[perm.length-1]+"|)";
		System.out.println(outstr);
	}
	
	LandmarkGMM(int[][] iGIDs, int iitems, int ilandmarks){
		nitems=iitems;
		nlandmarks=ilandmarks;
		ntotal=nitems+nlandmarks;
		if(ntotal!=iGIDs[0].length) {
			System.out.println("Incorrect perm length for number of items");
			return;
		}
		GIDs=new int[iGIDs.length][ntotal];
		for(int n=0;n<iGIDs.length;)
		GIDs=iGIDs;
	}
	
	void fitPerm(int[] perm){
		if(GIDs==null) {
			System.out.println("No GIDs: rerun after using setGIDs");
		}
		int n=perm.length;
		double[] theta=new double[n-1];
		double negLogLike=0;
		for(int i=1;i<n;i++) {
			int[] ritems=new int[i];
			System.arraycopy(perm, n-i, ritems, 0, i);
			int litem=perm[n-i-1];
			System.out.println(litem+" "+Arrays.toString(ritems));
			if(getNegdLogLike(litem,ritems,epsilon)<0) {
				theta[theta.length-i]=epsilon;
			}
			else {
				double bot=epsilon;
				double top=1;
				while(getNegdLogLike(litem,ritems,top)>0) {
					bot=top;
					top=bot+1;
				}
				while(Math.abs(top-bot)>epsilon) {
					if(getNegdLogLike(litem,ritems,(top+bot)/2)>0)
						bot=(top+bot)/2;
					else
						top=(top+bot)/2;
				}
				theta[theta.length-i]=(top+bot)/2;
				negLogLike+=getNegLogLike(litem,ritems,theta[theta.length-i]);
			}
		}
		System.out.println(Arrays.toString(theta)+" "+negLogLike);
	}
	
	double getNegLogLike(int litem, int[] ritems, double theta) {
		double v=0;
		for(int r:ritems) {
			v+=Q[r][litem];
		}
		if(litem<nitems || litem==nitems+nlandmarks-1) {
			double denom=(1-Math.exp(-theta*(ritems.length+1)))/(1-Math.exp(-theta));
			double negLogLike=-theta*v-GIDs.length*Math.log(denom);
			return(negLogLike);
		}
		else {
			double[] denoms=new double[ritems.length+1];
			denoms[0]=1;
			int i=1;
			double[] denoms2=new double[ritems.length+1];
			denoms2[0]=1;
			int i2=1;
			for(int r:ritems) {
				denoms[i]=denoms[i-1]+Math.exp(-theta*i);
				denoms2[i]=-i*theta+Math.log(Math.exp(theta*(i+1))-1)-Math.log(Math.exp(theta)-1);
				i++;
			}
			
			int M;
			int[] Mcounts=new int[ritems.length];
			for(int n=0;n<Ms.length;n++) {
				M=0;
				for(int r:ritems) {
					if(r<nitems) {
						M+=Ms[n][litem-nitems][r];
					}
				}
				Mcounts[M]++;
			}
			double denom=0;
			double denom2=0;
			for(int m=0;m<Mcounts.length;m++) {
				denom+=Mcounts[m]*Math.log(denoms[m]);
				denom2+=Mcounts[m]*denoms2[m];
			}
			double negLogLike=-theta*v-denom;//ONLY NEED ONE OF THESE
			double negLogLike2=-theta*v-denom2;//CHOOSE!!!!
			return(negLogLike2);
		}
	}
	
	double getNegdLogLike(int litem, int[] ritems, double theta) {
		double v=0;
		for(int r:ritems) {
			v+=Q[r][litem];
		}
		if(litem<nitems || litem==nitems+nlandmarks-1) {
			double denom=-(ritems.length+1)/(1-Math.exp(theta*(ritems.length+1)))+1/(1-Math.exp(theta));
			double negdLogLike=-v-GIDs.length*denom;
			return(negdLogLike);
		}
		else {
			double[] denoms=new double[ritems.length+1];
			denoms[0]=1;
			int i=1;
			double[] denoms2=new double[ritems.length+1];
			denoms2[0]=0;
			for(int r:ritems) {
				denoms[i]=denoms[i-1]+Math.exp(-theta*i);
				denoms2[i]=denoms2[i-1]+i*Math.exp(-theta*i);
				i++;
			}
			
			int M;
			int[] Mcounts=new int[ritems.length];
			for(int n=0;n<Ms.length;n++) {
				M=0;
				for(int r:ritems) {
					if(r<nitems) {
						M+=Ms[n][litem-nitems][r];
					}
				}
				Mcounts[M]++;
			}
			double denom=0;
			for(int m=0;m<Mcounts.length;m++) {
				denom+=Mcounts[m]*1/denoms[m]*denoms2[m];
			}
			double negdLogLike=-v+denom;//ONLY NEED ONE OF THESE
			return(negdLogLike);
		}
	}
	
	/*
	static double getThetaStandard(double[][] LRCount, double vlbar, CRIM model) {
		double[][] dlZ=new double[2][LRCount[0].length];
		double[][] Z=new double[2][LRCount[0].length];
		//HARD CODED FOR n=10000
		double epsilon=0.0001;
		if(model.dLogLikelihood(dlZ, LRCount, vlbar, 0.0001, 10000.0)<0) {
			return 0;
		}
		else {
			double thetal=0;
			double thetah=1;
			while(model.dLogLikelihood(dlZ, LRCount, vlbar, thetah, 10000.0)>0)
				thetah*=2;
			while(thetah-thetal>epsilon) {
				if(model.dLogLikelihood(dlZ, LRCount, vlbar, (thetah+thetal)/2, 10000.0)/2>0)
					thetal=(thetah+thetal)/2;
				else
					thetah=(thetah+thetal)/2;
			}
			return((thetah+thetal)/2);
		}
	}*/
	
	void setMs() {
		Ms=new int[GIDs.length][nlandmarks-1][nitems];
		for(int n=0;n<GIDs.length;n++) {
			for(int k=nlandmarks-2;k>=0;k--) {
				for(int m=0;m<nitems;m++) {
					Ms[n][k][m]=(GIDs[n][nitems+nlandmarks-1-k]>GIDs[n][m]) ? 1 : 0;
				}
			}
		}
	}
	
	int[] getMVec(int litem, int[] ritems) {
		int[] M=new int[Ms.length];
		for(int r=0;r<ritems.length;r++) {
			if(ritems[r]<nitems) {
				for(int n=0;n<Ms.length;n++) {
					M[n]+=Ms[n][litem-nitems][ritems[r]];
				}
			}
		}
		return(M);
	}
	
	void setGIDs(int[][] perms) {
		GIDs=new int[perms.length][perms[0].length];
		
		for(int n=0;n<perms.length;n++) {
			for(int k=0;k<perms[0].length;k++) {
				GIDs[n][perms[n][k]]=k;
			}
		}
		
		Q=new int[nitems+nlandmarks][nitems+nlandmarks];
		
		for (int n = 0; n < GIDs.length; n++) {
			for (int i = 0; i < nitems+nlandmarks - 1; i++) {
				for (int j = i + 1; j < nitems+nlandmarks; j++) {
					if (GIDs[n][i] < GIDs[n][j])
						Q[i][j]++;
					if (GIDs[n][i] > GIDs[n][j])
						Q[j][i]++;
				}
			}
		}
	}
	
	int[][] samplePerms(int samplesize, Random rand){
		int[][] retSample=new int[samplesize][nitems+nlandmarks];
		double[][] vcdf=getModelVCDF();
		//CDF for landmarks is the same, see detail below
		
		
		for(int n=0;n<samplesize;n++) {
			retSample[n][ntotal-1]=perm[ntotal-1];
			int M=1;
			if(retSample[n][ntotal-1]>=nitems) {
				M=0;
			}
			for(int k=ntotal-2;k>=0;k--) {
				double vrand=rand.nextDouble();
				// We can just multiply random double
				// by max index CDF for a landmark - caps the v max
				if(perm[k]>=nitems) {
					vrand*=vcdf[k][M];
				}
				int v;
				for(v=0;vrand>vcdf[k][v];v++);

				if(perm[k]>=nitems)
					M=v;
				else if(v<=M)
					M++;
				//System.out.println(vrand+" "+Arrays.toString(vcdf[k])+" "+v+" "+M);
				int ind=k;
				for(;v>0;v--) {
					retSample[n][ind++]=retSample[n][ind];
				}
				retSample[n][ind]=perm[k];
				//System.out.println(Arrays.toString(retSample[n]));
			}
		}
		return(retSample);
	}
	
	int[][] samplePerms(int samplesize){
		Random rand=new Random();
		return(samplePerms(samplesize,rand));
	}
	
	double[][] getModelVCDF(){
		double[][] modelVCDF=new double[nitems+nlandmarks-1][];
		for(int n=0;n<nitems+nlandmarks-1;n++)
			modelVCDF[n]=standardVlCDF(nitems+nlandmarks-1-n,thetas[n]);
		return(modelVCDF);
	}
	
	static double[] standardVlCDF(int ritems, double theta){
		double[] ret=new double[ritems+1];
		double total=0;
		for(int i=0;i<=ritems;i++) {
			ret[i]=Math.exp(-i*theta);
			total+=ret[i];
		}
		for(int i=0;i<=ritems;i++)
			ret[i]/=total;
		for(int j=1;j<=ritems;j++) 
			ret[j]+=ret[j-1];
		if(ret[ritems]!=1)
			ret[ritems]=1;
		return ret;
	}
}
