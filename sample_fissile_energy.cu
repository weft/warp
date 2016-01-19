#include <cuda.h>
#include <stdio.h>
#include "datadef.h"
#include "warp_device.cuh"

__global__ void sample_fissile_energy_kernel( unsigned N , float a , float b , unsigned* rn_bank , float* E ){

	int tid = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid >= N){return;}

	// declare
	unsigned	rn;
	float 		rn1, rn2, k, l, m, x, y, val1, val2;

	// get some random numbers
	rn		= rn_bank[ tid ];
	rn1		= get_rand(&rn);
	rn2		= get_rand(&rn);		

	// sample according to 3rd Monte Carlo sampler R12
	k		= 1.0 + (b/(8.0*a));
	l		= (k+sqrtf(k*k-1))/a;
	m		= a*l-1;
	x		= -logf(rn1);
	y		= -logf(rn2);
	val1	= (y-m*(x+1))*(y-m*(x+1));
	val2	= b*l*x;

	// reject until criterion satisfied
	while(val1>val2){
		rn1		= get_rand(&rn);
		rn2		= get_rand(&rn);
		x		= -logf(rn1);
		y		= -logf(rn2);
		val1	= (y-m*(x+1))*(y-m*(x+1));
		val2	= b*l*x;
	}

	// write sampled energy to array
	E[tid]	=	l*x;

	// write current rn back to bank
	rn_bank[tid]	=	rn;

	//printf("%10.8E\n",l*x);

}

void sample_fissile_energy( unsigned NUM_THREADS,  unsigned N , float a , float b , unsigned* rn_bank , float* E){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	sample_fissile_energy_kernel <<< blks, NUM_THREADS >>> (  N , a , b, rn_bank, E );
	cudaThreadSynchronize();

}

