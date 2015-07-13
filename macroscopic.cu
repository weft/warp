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

	// declare
	float 		norm[3];
	float 		samp_dist 		= 0.0;
	float 		cum_prob 		= 0.0;
	float 		diff			= 0.0;
	unsigned 	tope 			= 999999999;
	float 		epsilon 		= 2.0e-5;
	unsigned 	isdone 			= 0;
	float 		dotp 			= 0.0;

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
	memcpy(norm,space[tid].norm,3*sizeof(float));
//	norm[0]=0;
//	norm[1]=0;
//	norm[2]=0;

	float macro_t_total = 0.0;
	float e0 = main_E_grid[dex];
	float e1 = main_E_grid[dex+1];
	float t0,t1,number_density,surf_minimum, xhat_new, yhat_new, zhat_new;

	if(this_mat>=n_materials){
		printf("MACRO - this_mat %u > n_materials %u!!!!!\n",this_mat,n_materials);
        	rxn[tid_in]                     = 1001;  
        	isonum[tid]                     = 0;
		return;
	}

	if(dex>178889){printf("tid_in %u -> tid %u, dex %u\n",tid_in, tid, dex);}
	if(this_mat<0 | this_mat>2){printf("MATERIAL INVALID, this_mat = %u\n",this_mat);}
	if(n_isotopes<0 | n_isotopes>15){printf("N_ISOTOPES INVALID, n_isotopes = %u\n",n_isotopes);}
	if(n_columns<0 | n_columns>394){printf("N_COLUMNS INVALID, n_columns = %u\n",n_columns);}
	//if(tid_in >= 370899 & tid_in <=370901){printf("this_mat %u n_isotopes %u n_columns %u dex %u\n", this_mat, n_isotopes, n_columns, dex);}
	//if(tid_in >= 4635554 | tid_in <= 4635557){printf("n_isotopes %u this_mat %u\n",n_isotopes,this_mat);}

	// compute the total macroscopic cross section for this material
	for(int k=0; k<n_isotopes; k++){
		//printf("index %u\n", (n_isotopes*this_mat+k));
		number_density = material_matrix[n_isotopes*this_mat+k];
		if(number_density > 0.0){
			//lienarly interpolate
			//printf("val % 6.4E\n",material_matrix[n_isotopes*this_mat + k]);
			t0 = xs_data_MT[n_columns* dex    + k];     //dex is the row number
			t1 = xs_data_MT[n_columns*(dex+1) + k];
			macro_t_total += ( (t1-t0)/(e1-e0)*(this_E-e0) + t0 ) * number_density;    //interpolated micro times number density
			//printf("mat %u - density of tope %u = %6.3E\n",this_mat,k,material_matrix[n_isotopes*this_mat+k]);
		}
	}

	//if(N==1228932){printf("tid_in %u -> tid %u, dex %u\n",tid_in, tid, dex);}

	// compute the interaction length
	samp_dist = -logf(get_rand(&rn))/macro_t_total;
	float rn1 = get_rand(&rn);
	int flag = 0;
	// determine the isotope which the reaction occurs
	for(int k=0; k<n_isotopes; k++){
		number_density = material_matrix[n_isotopes*this_mat+k];
		if(number_density > 0.0){
			flag=1;
			//lienarly interpolate
			t0 = xs_data_MT[n_columns* dex    + k];     
			t1 = xs_data_MT[n_columns*(dex+1) + k];
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
               		number_density = material_matrix[n_isotopes*this_mat+k];
                	if(number_density > 0.0){
                	        //lienarly interpolate
                	        t0 = xs_data_MT[n_columns* dex    + k];
                	        t1 = xs_data_MT[n_columns*(dex+1) + k];
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
			isdone = 1;
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
			isdone = 0;
			tope=999999996;  // make reflection a different isotope 
		}
		else{ // if not outer cell, push through surface, set resample
			x += (surf_dist*xhat - 1.2*epsilon*norm[0]);
			y += (surf_dist*yhat - 1.2*epsilon*norm[1]);
			z += (surf_dist*zhat - 1.2*epsilon*norm[2]);
			this_rxn = 800;
			isdone = 0;
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
	done[tid] 			= isdone;


}

void macroscopic( unsigned NUM_THREADS,  unsigned N, unsigned Ntopes, unsigned n_materials, unsigned n_col , unsigned outer_cell, unsigned* remap, source_point * space, unsigned* isonum, unsigned* cellnum, unsigned * index, unsigned * matnum, unsigned* rxn, float * main_E_grid, unsigned * rn_bank, float * E, float * xs_data_MT , float* material_matrix, unsigned* done){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	macroscopic_kernel <<< blks, NUM_THREADS >>> ( N, Ntopes, n_materials, n_col, outer_cell, remap, space, isonum, cellnum, index, matnum, rxn, main_E_grid, rn_bank, E, xs_data_MT , material_matrix, done);
	cudaThreadSynchronize();

}

