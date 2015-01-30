#include <cuda.h>
#include <stdio.h>
#include "datadef.h"
#include "LCRNG.cuh"

__global__ void microscopic_kernel(unsigned N, unsigned n_isotopes, unsigned n_columns, unsigned* remap, unsigned* isonum, unsigned * index, float * main_E_grid, unsigned * rn_bank, float * E, float * xs_data_MT , unsigned * xs_MT_numbers_total, unsigned * xs_MT_numbers,  float* xs_data_Q, unsigned * rxn, float* Q, unsigned* done){


	int tid_in = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid_in >= N){return;}

	unsigned 	this_rxn 		= rxn[tid_in];
	if(this_rxn>=800){return;} //return if flagged to resample or leaked (leak can be in here since set by macro and remap hasn't been done)

	//remap
	int tid=remap[tid_in];
	//printf("tid %u remapped_tid %u\n",tid_in,tid);

	// load from array
	unsigned 	this_tope 		= isonum[tid];
	unsigned 	dex 			= index[tid];
	unsigned 	tope_beginning;
	unsigned 	tope_ending;
	unsigned 	this_dex;
	float 		this_E  		= E[tid];
	unsigned	rn 				= rn_bank[tid];
	float 		rn1 			= get_rand(&rn);
	float 		cum_prob 		= 0.0;
	float 		this_Q 			= 0.0;
	unsigned 	k 				= 0;
	
	//printf("tid %u dex %u topes %u rxn %u this_tope %u\n",tid,dex,n_isotopes,this_rxn,this_tope);

	if (this_tope == 0){  //first isotope
		tope_beginning = n_isotopes + 0;
		tope_ending    = n_isotopes + xs_MT_numbers_total[0]-1;
	}
	else{  //interior space
		tope_beginning = n_isotopes + xs_MT_numbers_total[this_tope-1];
		tope_ending    = n_isotopes + xs_MT_numbers_total[this_tope]-1;
	}

	//printf("tope,begin,end = %u %u %u\n",this_tope,tope_beginning,tope_ending);

	float xs_total = 0.0;
	float e0 = main_E_grid[dex];
	float e1 = main_E_grid[dex+1];
	float t0,t1;

	// compute the total microscopic cross section for this material
	// linearly interpolate, dex is the row number
	t0 			= xs_data_MT[n_columns* dex    + this_tope];     
	t1 			= xs_data_MT[n_columns*(dex+1) + this_tope];
	xs_total 	= (t1-t0)/(e1-e0)*(this_E-e0) + t0 ;    

	// determine the reaction for this isotope
	for(k=tope_beginning; k<=tope_ending; k++){
		//linearly interpolate
		t0 = xs_data_MT[n_columns* dex    + k];     
		t1 = xs_data_MT[n_columns*(dex+1) + k];
		cum_prob += ( (t1-t0)/(e1-e0)*(this_E-e0) + t0 ) / xs_total;
		if(rn1 <= cum_prob){
			// reactions happen in reaction k
			this_rxn = xs_MT_numbers[k];
			//printf("tope %u beg/end %u %u rxn %u cum_prob %6.4E rn1 %6.4E this_E %6.4E (tot,es,91,abs) %6.4E %6.4E %6.4E %6.4E\n",this_tope,tope_beginning,tope_ending,this_rxn,cum_prob,rn1,this_E,xs_data_MT[n_columns* dex    + 0],xs_data_MT[n_columns* dex    + 1],xs_data_MT[n_columns* dex    + 46],xs_data_MT[n_columns* dex    + 47]);
			this_Q   = xs_data_Q[k];
			this_dex = n_columns* dex + k;
			break;
		}
	}

	if(this_rxn == 999999999){ // there is a gap in between the last MT and the total cross section, remap the rn to fit into the available data (effectively rescales the total cross section so everything adds up to it, if things aren't samples the first time around)
		printf("micro - REACTION NOT SAMPLED CORRECTLY! tope=%u E=%10.8E dex=%u rxn=%u cum_prob=%30.28E rn1=%30.28E\n",this_tope, this_E, dex, this_rxn, cum_prob,rn1); //most likely becasue rn1=1.0
		rn1 = rn1 * cum_prob;
		cum_prob = 0.0;
		for(k=tope_beginning; k<tope_ending; k++){
			//lienarly interpolate
			t0 = xs_data_MT[n_columns* dex    + k];     
			t1 = xs_data_MT[n_columns*(dex+1) + k];
			cum_prob += ( (t1-t0)/(e1-e0)*(this_E-e0) + t0 ) / xs_total;
			if(rn1 <= cum_prob){
				// reactions happen in reaction k
				this_rxn = xs_MT_numbers[k];
				this_Q   = xs_data_Q[k];
				this_dex = n_columns * dex + k;
				break;
			}
		}
	}


	// write results out
	//if( this_rxn >= 811 & this_rxn<850 & this_rxn!=818 ){printf("microscopic sampled tid %u rxn %d energy %6.4E\n",tid,this_rxn,this_E);}
	//printf("%u\n",this_rxn);
	if(this_rxn==4){printf("MT=4!!!\n");}
	if(this_rxn==56){printf("MT=56 at E %10.8E\n",this_E);}
	rxn[tid_in] = this_rxn;
	Q[tid] 	 = this_Q;
	rn_bank[tid] = rn;
	//also write MT array index to dex instead of energy vector index
	index[tid] = this_dex;


}

void microscopic( unsigned NUM_THREADS, unsigned N, unsigned n_isotopes, unsigned n_columns, unsigned* remap,unsigned* isonum, unsigned * index, float * main_E_grid, unsigned * rn_bank, float * E, float * xs_data_MT , unsigned * xs_MT_numbers_total, unsigned * xs_MT_numbers,  float* xs_data_Q, unsigned * rxn, float* Q, unsigned* done){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	microscopic_kernel <<< blks, NUM_THREADS >>> ( N,  n_isotopes, n_columns, remap, isonum, index, main_E_grid, rn_bank, E, xs_data_MT , xs_MT_numbers_total, xs_MT_numbers, xs_data_Q, rxn, Q, done);
	cudaThreadSynchronize();

}

