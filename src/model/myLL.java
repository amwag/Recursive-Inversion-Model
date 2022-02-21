package model;

import java.util.Arrays;

public class myLL{
	//Need a head node with no data
	// to properly append
	myLL prev;
	myLL next;
	double loglikelihood;
	double[] thetas;
	int[] items;
	
	myLL(int[] items_in, double[] thetas_in, double likelihood_in){
		this.loglikelihood=likelihood_in;
		this.thetas=thetas_in;
		this.items=items_in;
	}
	
	void appendTo(myLL node) {
		while(node.loglikelihood>this.loglikelihood) {
			if(node.next==null) {
				node.next=this;
				this.prev=node;
				return;
			}
			else
				node=node.next;
		}
		myLL last=node.prev;
		this.prev=last;
		this.next=node;
		last.next=this;
		node.prev=this;
	}
	
	public Boolean hasNext() {
		if(next!=null)
			return(true);
		else
			return(false);
	}
	
	public String toString() {
		return("LL: "+loglikelihood+" | Items: "+Arrays.toString(items))+" | Thetas: "+Arrays.toString(thetas);
	}
	
	public Boolean contains(int i) {
		for(int k:items)
			if(k==i)
				return(true);
		return(false);
	}
}
