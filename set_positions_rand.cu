#include <cuda.h>
#include <stdio.h>
#include "datadef.h"
#include "warp_device.cuh"

__global__ void set_positions_rand_kernel(unsigned N , unsigned outer_cell_type, spatial_data * positions_ptr , unsigned * rn_bank , float x_min , float y_min , float z_min , float x_max , float y_max , float z_max ){

	int tid = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid>=N){return;}

	unsigned 	rn = rn_bank[tid];
	float rn1 = get_rand(&rn);
	float rn2 = get_rand(&rn);
	float rn3 = get_rand(&rn);
	float mu    = 2.0 * get_rand(&rn)-1.0;
	float theta = 6.28318530718 * get_rand(&rn);

	// set surface distance very far
	positions_ptr[tid].surf_dist =     500000; 

	// set directions isotropic
	positions_ptr[tid].xhat      =     sqrtf(1.0-mu*mu) * cosf( theta );
	positions_ptr[tid].yhat      =     sqrtf(1.0-mu*mu) * sinf( theta );
	positions_ptr[tid].zhat      =     mu;

	// set positions
	if(outer_cell_type==0){  // cube
		positions_ptr[tid].x         =     0.9 * ( ( x_max - x_min ) * rn1 + x_min );  
		positions_ptr[tid].y         =     0.9 * ( ( y_max - y_min ) * rn2 + y_min );  
		positions_ptr[tid].z         =     0.9 * ( ( z_max - z_min ) * rn3 + z_min ); 
	}
	else if(outer_cell_type==1){   //cylinder
    	float t = 6.28318530718 * rn1;
    	float r = x_max * sqrtf(rn2);
		positions_ptr[tid].x         =     0.9 * r * cosf(t);  
		positions_ptr[tid].y         =     0.9 * r * sinf(t);  
		positions_ptr[tid].z         =     0.9 * ( ( z_max - z_min ) * rn3 + z_min ); 
	}
	else if(outer_cell_type==2){  // hex
		positions_ptr[tid].x         =     1.4 * ( ( x_max - x_min ) * rn1 + x_min );  
		positions_ptr[tid].y         =     1.4 * ( ( y_max - y_min ) * rn2 + y_min );  
		positions_ptr[tid].z         =     1.4 * ( ( z_max - z_min ) * rn3 + z_min ); 
	}
	else if(outer_cell_type==3){  // sphere
		float t = 6.28318530718 * rn1;
    	float p = 2.0 * rn2 - 1.0;
    	float r = x_max * cbrtf(rn3);
		positions_ptr[tid].x         =     0.9 * r * sqrt(1.0-p*p) * cos(t);  
		positions_ptr[tid].y         =     0.9 * r * sqrt(1.0-p*p) * sin(t);
		positions_ptr[tid].z         =     0.9 * r * p; 
	}

	// update rn after using it
	rn_bank[tid]	=	rn;
  
	

}

void set_positions_rand( unsigned NUM_THREADS, unsigned N, unsigned outer_cell_type, spatial_data * d_space , unsigned * d_rn_bank, float * outer_cell_dims){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	set_positions_rand_kernel <<<  blks, NUM_THREADS >>> ( N , outer_cell_type, d_space , d_rn_bank, outer_cell_dims[0], outer_cell_dims[1], outer_cell_dims[2], outer_cell_dims[3], outer_cell_dims[4], outer_cell_dims[5]);
	cudaThreadSynchronize();

}