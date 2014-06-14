#include <cuda.h>
#include <stdio.h>
#include "datadef.h"

__global__ void sample_fission_spectra_kernel(unsigned N, unsigned RNUM_PER_THREAD, unsigned* index, unsigned* done, unsigned* active, unsigned* rxn, float* rn_bank , float * E, source_point* space, float** energydata){

	int tid = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid >= N){return;}
	if (rxn[tid]!=18){return;}
	//if(done[tid]){return;}

	// external data
	unsigned 		dex  		= index[tid];
	float 			this_E 		= E[tid];
	float * 		this_array 	= energydata[dex];

	// internal data
	float 		Emin=1e-11;
	float 		Emax=20.0;
	unsigned 	n, offset, vlen, next_vlen, law;
	float 		sampled_E, phi, mu, rn1, rn2, last_E, next_E;
	float 		cdf0, cdf1, e0, e1, m, pdf0, pdf1, arg;
	const float  pi           =   3.14159265359 ;

	// sample spectrum, set data.  
	// reset self then write elsewhere

	//read in values
	rn1 = rn_bank[ tid*RNUM_PER_THREAD + 11 ];
	rn2 = rn_bank[ tid*RNUM_PER_THREAD + 12 ];
	offset = 5;
	//printf("rxn %u eptr %p\n",this_rxn,this_array);
	memcpy(&last_E,   	&this_array[0], sizeof(float));
	memcpy(&next_E,   	&this_array[1], sizeof(float));
	memcpy(&vlen,   	&this_array[2], sizeof(float));
	memcpy(&next_vlen,	&this_array[3], sizeof(float));
	memcpy(&law, 		&this_array[4], sizeof(float));
	//printf("dex=%u rxn=%u law=%u vlen/next= %u %u, E-last/this/next= %6.4E %6.4E %6.4E\n",dex,rxn[tid],law,vlen,next_vlen,last_E,this_E,next_E);
	//sample energy dist
	sampled_E = 0.0;
	if(  rn2 <= (next_E-this_E)/(next_E-last_E) ){   //sample last E
		for ( n=0 ; n<vlen-1 ; n++ ){
			cdf0 		= this_array[ (offset +   vlen ) + n+0];
			cdf1 		= this_array[ (offset +   vlen ) + n+1];
			pdf0		= this_array[ (offset + 2*vlen ) + n+0];
			pdf1		= this_array[ (offset + 2*vlen ) + n+1];
			e0  		= this_array[ (offset          ) + n+0];
			e1  		= this_array[ (offset          ) + n+1]; 
			if( rn1 >= cdf0 & rn1 < cdf1 ){
				break;
			}
		}
	}
	else{
		for ( n=0 ; n<next_vlen-1 ; n++ ){
			cdf0 		= this_array[ (offset + 3*vlen +   next_vlen ) + n+0];
			cdf1  		= this_array[ (offset + 3*vlen +   next_vlen ) + n+1];
			pdf0		= this_array[ (offset + 3*vlen + 2*next_vlen ) + n+0];
			pdf1		= this_array[ (offset + 3*vlen + 2*next_vlen ) + n+1];
			e0   		= this_array[ (offset + 3*vlen               ) + n+0];
			e1   		= this_array[ (offset + 3*vlen               ) + n+1];
			if( rn1 >= cdf0 & rn1 < cdf1 ){
				break;
			}
		}
	}

	// interpolate the values
	m 			= (pdf1 - pdf0)/(e1-e0);
	arg = pdf0*pdf0 + 2.0 * m * (rn1-cdf0);
	if(arg<0){arg=0.0;}
	sampled_E 	= e0 + (  sqrtf( arg ) - pdf0) / m ;
	//sampled_E = e0 + (rn1-cdf0)/pdf0;
	//printf("%u %u %u %u %u %p %6.4E %u %u %6.4E %6.4E %6.4E %6.4E %6.4E %6.4E %6.4E %6.4E %6.4E\n",tid,tid*RNUM_PER_THREAD + 12,fork,n,dex,this_array,rn1,next_vlen,vlen,this_E,e0,e1,cdf0,cdf1,pdf0,pdf1,m,sampled_E);


	//sample isotropic directions
	rn1 = rn_bank[ tid*RNUM_PER_THREAD + 13 ];
	rn2 = rn_bank[ tid*RNUM_PER_THREAD + 14 ];
	mu  = 2.0*rn1-1.0; 
	phi = 2.0*pi*rn2;
	
	//check limits
	if (sampled_E >= Emax){sampled_E = Emax * 0.9;}//printf("enforcing limits in pop data_dex=%u, sampled_E = %6.4E\n",data_dex,sampled_E);}
	if (sampled_E <= Emin){sampled_E = Emin * 1.1;}//printf("enforcing limits in pop data_dex=%u, sampled_E = %6.4E\n",data_dex,sampled_E);}

	// sync before writes
	__syncthreads();

	// set self data
	E    [ tid ] 		= sampled_E;
	space[ tid ].xhat 	= sqrtf(1.0-(mu*mu))*cosf(phi);
	space[ tid ].yhat 	= sqrtf(1.0-(mu*mu))*sinf(phi); 
	space[ tid ].zhat 	= mu;
	//done [ tid ] 		= 0;
	//rxn  [ tid ] 		= 0;

}

void sample_fission_spectra( unsigned NUM_THREADS, unsigned RNUM_PER_THREAD, unsigned N, unsigned* index, unsigned* done, unsigned* active, unsigned* rxn, float * rn_bank, float * E, source_point* space, float** energydata){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	sample_fission_spectra_kernel <<< blks, NUM_THREADS >>> (  N, RNUM_PER_THREAD, index, done, active, rxn, rn_bank, E, space, energydata );
	cudaThreadSynchronize();

}

