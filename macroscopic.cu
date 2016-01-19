#include <cuda.h>
#include <stdio.h>
#include "datadef.h"
#include "LCRNG.cuh"

__global__ void macroscopic_kernel(unsigned N, unsigned n_materials, cross_section_data* d_xsdata, particle_data* d_particles, unsigned* d_remap, float* d_number_density_matrix){


	int tid_in = threadIdx.x+blockIdx.x*blockDim.x; 
	if (tid_in >= N){return;}

	// declare shared variables
	__shared__ 	unsigned			n_isotopes;				
	//__shared__ 	unsigned			energy_grid_len;		
	__shared__ 	unsigned			total_reaction_channels;
	//__shared__ 	unsigned*			rxn_numbers;			
	//__shared__ 	unsigned*			rxn_numbers_total;		
	__shared__ 	float*				energy_grid;			
	//__shared__ 	float*				rxn_Q;						
	__shared__ 	float*				xs;						
	//__shared__ 	float*				awr;					
	//__shared__ 	float*				temp;					
	//__shared__ 	dist_container*		dist_scatter;			
	//__shared__ 	dist_container*		dist_energy; 
	__shared__	spatial_data*		space;	
	__shared__	unsigned*			rxn;	
	__shared__	float*				E;		
	//__shared__	float*				Q;		
	__shared__	unsigned*			rn_bank;
	//__shared__	unsigned*			cellnum;
	__shared__	unsigned*			matnum;	
	__shared__	unsigned*			isonum;	
	//__shared__	unsigned*			yield;	
	//__shared__	float*				weight;	
	__shared__	unsigned*			index;	

	// have thread 0 of block copy all pointers and static info into shared memory
	if (threadIdx.x == 0){
		n_isotopes					= d_xsdata[0].n_isotopes;								
		//energy_grid_len				= d_xsdata[0].energy_grid_len;				
		total_reaction_channels		= d_xsdata[0].total_reaction_channels;
		//rxn_numbers 				= d_xsdata[0].rxn_numbers;						
		//rxn_numbers_total			= d_xsdata[0].rxn_numbers_total;					
		energy_grid 				= d_xsdata[0].energy_grid;						
		//rxn_Q 						= d_xsdata[0].Q;												
		xs 							= d_xsdata[0].xs;												
		//awr 						= d_xsdata[0].awr;										
		//temp 						= d_xsdata[0].temp;										
		//dist_scatter 				= d_xsdata[0].dist_scatter;						
		//dist_energy 				= d_xsdata[0].dist_energy; 
		space						= d_particles[0].space;
		rxn							= d_particles[0].rxn;
		E							= d_particles[0].E;
		//Q							= d_particles[0].Q;	
		rn_bank						= d_particles[0].rn_bank;
		//cellnum						= d_particles[0].cellnum;
		matnum						= d_particles[0].matnum;
		isonum						= d_particles[0].isonum;
		//yield						= d_particles[0].yield;
		//weight						= d_particles[0].weight;
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
	float 		cum_prob 		= 0.0;
	float 		diff			= 0.0;
	unsigned 	tope 			= 999999999;
	float 		epsilon 		= 2.0e-5;
	float 		dotp 			= 0.0;
	float 		macro_t_total 	= 0.0;
	unsigned 	n_columns 		= n_isotopes + total_reaction_channels;

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

	float e0 = energy_grid[dex];
	float e1 = energy_grid[dex+1];
	float t0,t1,number_density,surf_minimum, xhat_new, yhat_new, zhat_new;

	if(this_mat>=n_materials){
		printf("MACRO - this_mat %u > n_materials %u!!!!!\n",this_mat,n_materials);
        rxn[tid_in]   = 1001;  
        isonum[tid]   = 0;
		return;
	}

	if (this_rxn>801)printf("multiplicity %u entered macro at E %10.8E\n",this_rxn,this_E);

	// compute the total macroscopic cross section for this material
	for(int k=0; k<n_isotopes; k++){
		number_density = d_number_density_matrix[n_isotopes*this_mat+k];
		if(number_density > 0.0){
			//lienarly interpolate
			t0 = xs[n_columns* dex    + k];     //dex is the row number
			t1 = xs[n_columns*(dex+1) + k];
			macro_t_total += ( (t1-t0)/(e1-e0)*(this_E-e0) + t0 ) * number_density;    //interpolated micro times number density
		}
	}

	// compute the interaction length
	samp_dist = -logf(get_rand(&rn))/macro_t_total;
	float rn1 = get_rand(&rn);
	int flag = 0;
	// determine the isotope which the reaction occurs
	for(int k=0; k<n_isotopes; k++){
		number_density = d_number_density_matrix[n_isotopes*this_mat+k];
		if(number_density > 0.0){
			flag=1;
			//lienarly interpolate
			t0 = xs[n_columns* dex    + k];     
			t1 = xs[n_columns*(dex+1) + k];
			cum_prob += ( ( (t1-t0)/(e1-e0)*(this_E-e0) + t0 ) * number_density ) / macro_t_total;
			if( k==n_isotopes-1 & cum_prob<1.0){cum_prob=1.0;}  //sometimes roundoff makes this a problem
			if( rn1 <= cum_prob){
				// reactions happen in isotope k
				tope = k;
				break;
			}
		}
	}
	if(tope == 999999999){ 
		printf("macro - ISOTOPE NOT SAMPLED CORRECTLY!  Resampling with scaled probability... tid %u rn1 %10.8E cum_prob %10.8E flag %d this_mat %u rxn %u\n",tid,rn1,cum_prob,flag,this_mat,this_rxn);
		rn1 = rn1 * cum_prob;
		for(int k=0; k<n_isotopes; k++){
               		number_density = d_number_density_matrix[n_isotopes*this_mat+k];
                	if(number_density > 0.0){
                	        //lienarly interpolate
                	        t0 = xs[n_columns* dex    + k];
                	        t1 = xs[n_columns*(dex+1) + k];
                	        cum_prob += ( ( (t1-t0)/(e1-e0)*(this_E-e0) + t0 ) * number_density ) / macro_t_total;
                	        if( k==n_isotopes-1 & cum_prob<1.0){cum_prob=1.0;}  //sometimes roundoff makes this a problem
                	        if( rn1 <= cum_prob){
                	                // reactions happen in isotope k
                	                tope = k;
                	                break;
                	        }
                	}
        	}
	}
	if(tope == 999999999){
                printf("macro - ISOTOPE  NOT SAMPLED CORRECTLY AFTER RESAMPLING WITH SCALED PROBABILITY! tid=%u rn1 %10.8E cum_prob %10.8E\n",tid,rn1,cum_prob);
	}

	// do surf/samp compare, calculate epsilon projection onto neutron trajectory
	diff = surf_dist - samp_dist;
	dotp = norm[0]*xhat + norm[1]*yhat + norm[2]*zhat;
	surf_minimum = -epsilon / dotp;   // dotp *should* never be zero since optix won't report parallel lines
	if (surf_minimum>surf_dist){surf_minimum = surf_dist;}

	// complain if the dot product is positive.  normal should always be in the reflective sense, and normal should always be negative then.
	if (dotp > 0.0 | dotp < -1.0){
		//printf("norms(%u,:)=[%6.4E,%6.4E,%6.4E,%6.4E,%6.4E,%6.4E];\n",tid_in+1,x+surf_dist*xhat,y+surf_dist*yhat,z+surf_dist*zhat,norm[0],norm[1],norm[2]);
		printf("!!! (surface normal) dot (neutron direction) is positive or < -1!   dotp = %10.8E\n",dotp);
	}
	
	// surface logic
	if( diff < surf_minimum ){  //move to surface, set resample flag, push neutron epsilon away from surface so backscatter works right
		// enforce BC
		if (enforce_BC == 1){  // black BC
			x += (surf_dist + 2.1*surf_minimum) * xhat;
			y += (surf_dist + 2.1*surf_minimum) * yhat;
			z += (surf_dist + 2.1*surf_minimum) * zhat;
			this_rxn  = 999;
			tope=999999997;  // make leaking a different isotope than resampling
		}
		else if(enforce_BC == 2){  // specular reflection BC
			// move epsilon off of surface
			x += ((surf_dist*xhat) + 1.2*epsilon*norm[0]);
			y += ((surf_dist*yhat) + 1.2*epsilon*norm[1]);
			z += ((surf_dist*zhat) + 1.2*epsilon*norm[2]);
			// calculate reflection
			xhat_new = -(2.0 * dotp * norm[0]) + xhat; 
			yhat_new = -(2.0 * dotp * norm[1]) + yhat; 
			zhat_new = -(2.0 * dotp * norm[2]) + zhat; 
			// flags
			this_rxn = 801;  // reflection is 801 
			tope=999999996;  // make reflection a different isotope 
		}
		else{ // if not outer cell, push through surface, set resample
			x += (surf_dist*xhat - 1.2*epsilon*norm[0]);
			y += (surf_dist*yhat - 1.2*epsilon*norm[1]);
			z += (surf_dist*zhat - 1.2*epsilon*norm[2]);
			this_rxn = 800;
			tope=999999998;  // make resampling a different isotope than mis-sampling
		}
	}
	else{  //move to sampled distance (adjust if within epsilon of boundary), null reaction
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

	printf("tope[%d]=%u\n", tid,tope);

}

void macroscopic(unsigned NUM_THREADS, unsigned N, unsigned n_materials, cross_section_data* d_xsdata, particle_data* d_particles, unsigned* d_remap, float* d_number_density_matrix ){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	macroscopic_kernel <<< blks, NUM_THREADS >>> ( N, n_materials, d_xsdata, d_particles, d_remap, d_number_density_matrix);
	cudaThreadSynchronize();

}

