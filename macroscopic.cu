#include <cuda.h>
#include <stdio.h>
#include "datadef.h"
#include "LCRNG.cuh"

__global__ void macroscopic_kernel(unsigned N, unsigned n_isotopes, unsigned n_materials, unsigned n_columns, unsigned outer_cell, unsigned* remap, source_point * space, unsigned* isonum, unsigned* cellnum, unsigned * index, unsigned * matnum, unsigned* rxn, float * main_E_grid, unsigned * rn_bank, float * E, float * xs_data_MT , float* material_matrix, unsigned* done){


	int tid_in = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid_in >= N){return;}

	// return if terminated
	unsigned this_rxn=rxn[tid_in];
	if (this_rxn>900){return;}

	//remap
	int tid=remap[tid_in];

	// load from arrays
	unsigned 	this_mat 		= matnum[tid];
	unsigned 	dex 			= index[tid];   
	unsigned 	rn 				= rn_bank[tid];
	//unsigned 	cell 			= cellnum[tid];
	float 		this_E  		= E[tid];
	float		x 				= space[tid].x;
	float		y 				= space[tid].y;
	float		z 				= space[tid].z;
	float		xhat 			= space[tid].xhat;
	float		yhat 			= space[tid].yhat;
	float		zhat 			= space[tid].zhat;
	float		surf_dist 		= space[tid].surf_dist;
	unsigned 	enforce_BC 		= space[tid].enforce_BC;  
	float 		norm[3];
	memcpy(norm,space[tid].norm,3*sizeof(float));
	//norm[0] = 0.0;
	//norm[1] = 0.0;
	//norm[2] = 0.0;
	float 		samp_dist 		= 0.0;
	float 		cum_prob 		= 0.0;
	float 		diff			= 0.0;
	unsigned 	tope 			= 999999999;
	float 		epsilon 		= 5.0e-4;
	//unsigned    this_rxn		= 0;
	//unsigned    this_rxn = rxn[tid_in];
	unsigned 	isdone 			= 0;


	float macro_t_total = 0.0;
	float e0 = main_E_grid[dex];
	float e1 = main_E_grid[dex+1];
	float t0,t1;

	__syncthreads();

	// compute the total macroscopic cross section for this material
	for(int k=0; k<n_isotopes; k++){
		//lienarly interpolate
		//printf("val % 6.4E\n",s_material_matrix[n_isotopes*this_mat + k]);
		t0 = xs_data_MT[n_columns* dex    + k];     //dex is the row number
		t1 = xs_data_MT[n_columns*(dex+1) + k];
		macro_t_total += ( (t1-t0)/(e1-e0)*(this_E-e0) + t0 ) * material_matrix[n_isotopes*this_mat+k];    //interpolated micro times number density
		//printf("mat %u - density of tope %u = %6.3E\n",this_mat,k,material_matrix[n_isotopes*this_mat+k]);
	}

	// compute the interaction length
	samp_dist = -logf(get_rand(&rn))/macro_t_total;
	float rn1 = get_rand(&rn);

	// determine the isotope which the reaction occurs
	for(int k=0; k<n_isotopes; k++){
		//lienarly interpolate
		t0 = xs_data_MT[n_columns* dex    + k];     
		t1 = xs_data_MT[n_columns*(dex+1) + k];
		cum_prob += ( ( (t1-t0)/(e1-e0)*(this_E-e0) + t0 ) * material_matrix[n_isotopes*this_mat+k] ) / macro_t_total;
		if( k==n_isotopes-1 & cum_prob<1.0){cum_prob=1.0;}  //sometimes roundoff makes this a problem
		if( rn1 <= cum_prob){
			// reactions happen in isotope k
			tope = k;
			break;
		}
	}
	if(tope == 999999999){ 
		printf("macro - ISOTOPE NOT SAMPLED CORRECTLY! tope=%u E=%10.8E dex=%u mat=%u rn=%u cum_prob=%12.10E s_mm=%12.10E\n",tope, this_E, dex, this_mat, rn, cum_prob,material_matrix[n_isotopes*this_mat + tope]);
	}

	// do surf/samp compare
	//printf("norm %6.4E %6.4E %6.4E\n",norm[0],norm[1],norm[2]);
	//printf("% 6.4E % 6.4E % 6.4E % 6.4E % 6.4E % 6.4E % 6.4E % 6.4E %u\n",x,y,z,xhat,yhat,zhat,surf_dist,samp_dist,enforce_BC);
	diff = surf_dist - samp_dist;
	if( diff < 0.0 ){  //move to surface, set resample flag, push neutron epsilon away from surface so backscatter works right
		x += (surf_dist * xhat  +  1.2 * epsilon * norm[0]);
		y += (surf_dist * yhat  +  1.2 * epsilon * norm[1]);
		z += (surf_dist * zhat  +  1.2 * epsilon * norm[2]);
		this_rxn = 800;
		tope=999999998;  // make resampling a different isotope than mis-sampling
		// enforce BC
		if (enforce_BC){
			isdone = 1;
			this_rxn  = 999;
			tope=999999997;  // make leaking a different isotope than resampling
			//printf("leaked tid %u xyz % 6.4E % 6.4E % 6.4E dir % 6.4E % 6.4E % 6.4E\n",tid,x,y,z,xhat,yhat,zhat);
		}
	}
	else{  //move to sampled distance, null reaction
		x += samp_dist * xhat;
		y += samp_dist * yhat;
		z += samp_dist * zhat;
		this_rxn = 0;
	}

	//if(this_rxn==800){printf("resmapled tid %u xyz % 6.4E % 6.4E % 6.4E\n",tid,x,y,z);}
	//write outputs
	space[tid].x 			= x;
	space[tid].y			= y;
	space[tid].z			= z;
	space[tid].macro_t 		= macro_t_total;
	rxn[tid_in] 			= this_rxn;  // rxn is sorted WITH the remapping vector, i.e. its index does not need to be remapped
	isonum[tid] 			= tope;
	rn_bank[tid] 			= rn;
	done[tid] 				= isdone;


}

void macroscopic( unsigned NUM_THREADS,  unsigned N, unsigned Ntopes, unsigned n_materials, unsigned n_col , unsigned outer_cell, unsigned* remap, source_point * space, unsigned* isonum, unsigned* cellnum, unsigned * index, unsigned * matnum, unsigned* rxn, float * main_E_grid, unsigned * rn_bank, float * E, float * xs_data_MT , float* material_matrix, unsigned* done){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	macroscopic_kernel <<< blks, NUM_THREADS >>> ( N, Ntopes, n_materials, n_col, outer_cell, remap, space, isonum, cellnum, index, matnum, rxn, main_E_grid, rn_bank, E, xs_data_MT , material_matrix, done);
	cudaThreadSynchronize();

}

