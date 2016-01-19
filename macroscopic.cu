#include <cuda.h>
#include <stdio.h>
#include "datadef.h"
#include "warp_device.cuh"

__global__ void macroscopic_kernel(unsigned N, unsigned n_materials, cross_section_data* d_xsdata, particle_data* d_particles, unsigned* d_remap, float* d_number_density_matrix){


	int tid_in = threadIdx.x+blockIdx.x*blockDim.x; 
	if (tid_in >= N){return;}

	// declare shared variables
	__shared__ 	unsigned			n_isotopes;
	__shared__ 	unsigned			total_reaction_channels;
	__shared__ 	float*				energy_grid;
	__shared__ 	float*				xs;
	__shared__	spatial_data*		space;
	__shared__	unsigned*			rxn;
	__shared__	float*				E;
	__shared__	unsigned*			rn_bank;
	__shared__	unsigned*			matnum;
	__shared__	unsigned*			isonum;
	__shared__	unsigned*			index;

	// have thread 0 of block copy all pointers and static info into shared memory
	if (threadIdx.x == 0){
		n_isotopes					= d_xsdata[0].n_isotopes;
		total_reaction_channels		= d_xsdata[0].total_reaction_channels;
		energy_grid 				= d_xsdata[0].energy_grid
		xs 							= d_xsdata[0].xs;
		space						= d_particles[0].space;
		rxn							= d_particles[0].rxn;
		E							= d_particles[0].E;
		rn_bank						= d_particles[0].rn_bank;
		matnum						= d_particles[0].matnum;
		isonum						= d_particles[0].isonum;
		index						= d_particles[0].index;
	}

	// make sure shared loads happen before anything else
	__syncthreads();

	// return if terminated
	unsigned this_rxn=rxn[tid_in];
	if (this_rxn>=900){return;}

	//remap
	int tid=d_remap[tid_in];

	// declare
	float 		norm[3];
	float 		samp_dist 		= 0.0;
	float 		diff			= 0.0;
	unsigned 	tope 			= 999999999;
	float 		epsilon 		= 2.0e-5;
	float 		dotp 			= 0.0;
	float 		macro_t_total 	= 0.0;
	int   		flag 			= 0;
	float surf_minimum, xhat_new, yhat_new, zhat_new;

	// load from arrays
	unsigned 	this_mat 		=  matnum[tid];
	unsigned 	dex 			=   index[tid];  
	unsigned 	rn 				= rn_bank[tid];
	float 		this_E  		=       E[tid];
	float		x 				=   space[tid].x;
	float		y 				=   space[tid].y;
	float		z 				=   space[tid].z;
	float		xhat 			=   space[tid].xhat;
	float		yhat 			=   space[tid].yhat;
	float		zhat 			=   space[tid].zhat;
	float		surf_dist 		=   space[tid].surf_dist;
	unsigned 	enforce_BC 		=   space[tid].enforce_BC;  
	memcpy(norm,space[tid].norm,3*sizeof(float));

	// compute some things
	unsigned 	n_columns 		= n_isotopes + total_reaction_channels;
	float 		e0 				= energy_grid[dex];
	float 		e1 				= energy_grid[dex+1];

	if(this_mat>=n_materials){
		printf("MACRO - this_mat %u > n_materials %u!!!!!\n",this_mat,n_materials);
        rxn[tid_in]   = 1001;  
        isonum[tid]   = 0;
		return;
	}

	if (this_rxn>801)printf("multiplicity %u entered macro at E %10.8E\n",this_rxn,this_E);

	// compute the total macroscopic cross section for this material
	macro_t_total = compute_macro_t(	n_isotopes,
										e0, e1, this_E,
										&d_number_density_matrix[this_mat],  
										&xs[dex*n_columns],  
										&xs[(dex+1)*n_columns] 				);

	// compute the interaction length
	samp_dist = -logf(get_rand(&rn))/macro_t_total;

	// determine the isotope in the material for this cell
	tope = sample_isotope(	n_isotopes, macro_t_total, get_rand(&rn),
							e0, e1, this_E,
							&d_number_density_matrix[this_mat],  
							&xs[dex*n_columns],  
							&xs[(dex+1)*n_columns]					);


	// do surf/samp compare, calculate epsilon projection onto neutron trajectory
	// dotp positive = neutron is inside the cell (normal points out, trajectory must be coming from inside)
	// dotp negative = neutorn is outside the cell
	diff = surf_dist - samp_dist;
	dotp = norm[0]*xhat + norm[1]*yhat + norm[2]*zhat;
	surf_minimum = epsilon / fabsf(dotp);
	
	// surface logic
	if( diff < surf_minimum ){  // need to make some decisions so epsilon is handled correctly
		// if not outer cell, neutron placement is too close to surface.  risk of next interaction not seeing the surface due to epsilon.
		// preserve if in this cell or next, but make sure neutron is at least an epsilon away from the surface.
		if (diff < 0.0){ // next cell, enforce BC or push through
			if (enforce_BC == 1){  // black BC
				x += (surf_dist + 2.1*surf_minimum) * xhat;
				y += (surf_dist + 2.1*surf_minimum) * yhat;
				z += (surf_dist + 2.1*surf_minimum) * zhat;
				this_rxn  = 999;  // leaking is 999
				tope=999999997;  
			}
			else if(enforce_BC == 2){  // specular reflection BC
				// move epsilon off of surface
				x += ((surf_dist*xhat) + 1.2*epsilon*norm[0]); // do not need copysign since dotp should always be positive (neutron always inside the outer cell)
				y += ((surf_dist*yhat) + 1.2*epsilon*norm[1]);
				z += ((surf_dist*zhat) + 1.2*epsilon*norm[2]);
				// calculate reflection
				xhat_new = -(2.0 * dotp * norm[0]) + xhat; 
				yhat_new = -(2.0 * dotp * norm[1]) + yhat; 
				zhat_new = -(2.0 * dotp * norm[2]) + zhat; 
				// flags
				this_rxn = 801;  // reflection is 801 
				tope=999999996;  
			}
			else{   // next cell, move to intersection point, then move *out* epsilon along surface normal
				x += surf_dist*xhat + copysignf(1.0,dotp)*1.2*epsilon*norm[0];
				y += surf_dist*yhat + copysignf(1.0,dotp)*1.2*epsilon*norm[1];
				z += surf_dist*zhat + copysignf(1.0,dotp)*1.2*epsilon*norm[2];
				this_rxn = 800;  // resampling is 800
				tope=999999998;  
			}
		else{   // this cell, move to intersection point, then move *in* epsilon along surface normal 
			x += surf_dist*xhat - copysignf(1.0,dotp)*1.2*epsilon*norm[0];
			y += surf_dist*yhat - copysignf(1.0,dotp)*1.2*epsilon*norm[1];
			z += surf_dist*zhat - copysignf(1.0,dotp)*1.2*epsilon*norm[2];
			this_rxn = 0;
		}
	}
	else{  // near side of minimum, can simply move the neutron
			x += samp_dist * xhat;
			y += samp_dist * yhat;
			z += samp_dist * zhat;
			this_rxn = 0;
	}

	//write outputs
	space[tid].x 			= x;
	space[tid].y			= y;
	space[tid].z			= z;
	space[tid].macro_t 		= macro_t_total;
	if(enforce_BC==2){
		space[tid].xhat = xhat_new;  // write reflected directions for specular BC
		space[tid].yhat = yhat_new;
		space[tid].zhat = zhat_new;
	}
	rxn[tid_in] 			= this_rxn;  // rxn is sorted WITH the remapping vector, i.e. its index does not need to be remapped
	isonum[tid] 			= tope;
	rn_bank[tid] 			= rn;

	//printf("tope[%d]=%u\n", tid,tope);

}

void macroscopic(unsigned NUM_THREADS, unsigned N, unsigned n_materials, cross_section_data* d_xsdata, particle_data* d_particles, unsigned* d_remap, float* d_number_density_matrix ){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	macroscopic_kernel <<< blks, NUM_THREADS >>> ( N, n_materials, d_xsdata, d_particles, d_remap, d_number_density_matrix);
	cudaThreadSynchronize();

}

