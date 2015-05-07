#include <cuda.h>
#include <stdio.h>
#include "datadef.h"
#include "LCRNG.cuh"

__global__ void fission_kernel(unsigned N, unsigned starting_index, unsigned* remap, unsigned * rxn , unsigned * index, unsigned * yield , unsigned * rn_bank, unsigned* done, float** scatterdat){

	
	int tid_in = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid_in >= N){return;}       //return if out of bounds
	
	//remap to active
	int tid=remap[starting_index + tid_in];
	unsigned this_rxn = rxn[starting_index + tid_in];

	// print and return if wrong
	if (this_rxn < 811 | this_rxn > 845){printf("fission kernel accessing wrong reaction @ remapped dex %u sorted dex %u rxn %u\n",tid,starting_index + tid_in, this_rxn);return;} 

	//load rxn number, init values
	unsigned 	this_yield 	= 0;
	unsigned 	inu 		= 0;
	float 		nu 			= 0.0;
	unsigned	rn 			= rn_bank[ tid ];

	// determine integer yields
	if (this_rxn == 818 | this_rxn == 819 | this_rxn == 820 | this_rxn == 821){

		// load nu from arrays
		unsigned 	this_dex 	= index[tid];
	
		//load nu value, since e search has alrady been done, nu should be where the scatter array is (fission is always isotropic)
		memcpy(&nu, &scatterdat[this_dex], sizeof(float));
		if (nu==0.0){
			nu=2.8;
			printf("something is wrong with fission yields, nu = %6.4E, guessing %4.2f\n",0.0,nu); 
		}

		// get integer part
		inu = (unsigned) nu;
		
		// sample floor or ceil based on fractional part
		if((float)inu+get_rand(&rn) <= nu){
			this_yield = inu+1;
		}
		else{
			this_yield = inu;
		}

	}
	else if(this_rxn == 817 | this_rxn == 825 | this_rxn == 842){
		this_yield = 3;  
	}
	else if(this_rxn == 816 | this_rxn==824 | this_rxn == 811 | this_rxn == 824 | this_rxn == 829 | this_rxn == 830 | this_rxn == 841){
		this_yield = 2;  
	}
	else{
		this_yield = 1;
	}

	// write output and terminate history
	yield[tid] = this_yield;
	done[tid]  = 1;    // pop will re-activate this data slot on fixed-source runs (should be obsolete with active remapping, might need for fixed source?  I don't remember...)
 	rxn[starting_index + tid_in] = this_rxn+100;  //mark as done by putting it in the 900 block
	
	// write rn back if it was used for fission sampling
	if (this_rxn == 818 | this_rxn == 819 | this_rxn == 820 | this_rxn == 821){
		rn_bank[tid] = rn;  
	}  

}

void fission( cudaStream_t stream, unsigned NUM_THREADS, unsigned N, unsigned starting_index, unsigned* remap, unsigned * rxn , unsigned * index, unsigned * yield , unsigned * rn_bank, unsigned* done, float** scatterdat){

	if(N<1){return;}
	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	fission_kernel <<< blks, NUM_THREADS , 0 , stream >>> (   N,  starting_index, remap, rxn , index, yield , rn_bank, done, scatterdat);
	cudaThreadSynchronize();

}

