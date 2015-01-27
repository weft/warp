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
	//if(done[tid]){return;}

	// print and return if wrong
	if (this_rxn < 811 | this_rxn > 845){printf("fission kernel accessing wrong reaction @ remapped dex %u sorted dex %u rxn %u\n",tid,starting_index + tid_in, this_rxn);return;} 

	//load rxn number, init values
	unsigned 	this_yield 	= 0;
	unsigned 	inu 		= 0;
	float 		nu 			= 0.0;
	unsigned	rn 			= rn_bank[ tid ];

	//only do reactions with secondary neutrons
	//if (rxn[tid] == 18 | rxn[tid] == 16 | rxn[tid] == 17 | rxn[tid] == 37 | rxn[tid] == 24 | rxn[tid] == 41){}
	//else {return;} 

	//printf("in fission\n");

	if (this_rxn == 818){
		// load nu from arrays
		unsigned 	this_dex 	= index[tid];
	
		//load nu value, since e search has alrady been done!
		memcpy(&nu, &scatterdat[this_dex], sizeof(float));
		inu = (unsigned) nu;
	
		if((float)inu+get_rand(&rn) <= nu){
			this_yield = inu+1;
		}
		else{
			this_yield = inu;
		}
		//printf("nu %6.4E inu %u rn1 %6.4E yield %u\n",nu,inu,rn1,this_yield);
	}
	else if(this_rxn == 817){
		this_yield = 3;  
	}
	else if(this_rxn == 816 | this_rxn==824 | this_rxn == 841){
		this_yield = 2;
	}

	// write output and terminate history
	yield[tid] = this_yield;
	done[tid]  = 1;    // pop will re-activate this data slot on fixed-source runs
	rxn[starting_index + tid_in] = this_rxn+100;  //mark as done by putting it in the 900 block
	
	if (this_rxn == 818){
		rn_bank[tid] = rn;  //rn was used for fission sampling
	}  

}

void fission( cudaStream_t stream, unsigned NUM_THREADS, unsigned N, unsigned starting_index, unsigned* remap, unsigned * rxn , unsigned * index, unsigned * yield , unsigned * rn_bank, unsigned* done, float** scatterdat){

	if(N<1){return;}
	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	fission_kernel <<< blks, NUM_THREADS , 0 , stream >>> (   N,  starting_index, remap, rxn , index, yield , rn_bank, done, scatterdat);
	cudaThreadSynchronize();

}

