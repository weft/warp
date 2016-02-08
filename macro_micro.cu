#include <cuda.h>
#include <stdio.h>
#include "datadef.h"
#include "warp_device.cuh"
#include "check_cuda.h"

__global__ void macro_micro_kernel(unsigned N, unsigned converged, unsigned n_materials, unsigned n_tallies, cross_section_data* d_xsdata, particle_data* d_particles, tally_data* d_tally, unsigned* d_remap, float* d_number_density_matrix){
/*
This kernel does a lot.  It does all the total interaction processes:  
	samples the distance to the next interaction, 
	determines if the interaction needs to be resampled if it is past the nearest surface
	scores a flux tally if the current cell is flagged for one, and
	determines the reaction type (if interaction distance doesn't need to be resampled)
All neutrons need these things done, so these routines all live in the same routine.
*/

	// declare shared variables
	__shared__ 	unsigned			n_isotopes;				
	__shared__ 	unsigned			energy_grid_len;		
	__shared__ 	unsigned			total_reaction_channels;
	__shared__ 	unsigned*			rxn_numbers;			
	__shared__ 	unsigned*			rxn_numbers_total;		
	__shared__ 	float*				energy_grid;			
	__shared__ 	float*				rxn_Q;						
	__shared__ 	float*				xs;						
	//__shared__ 	float*				awr;					
	//__shared__ 	float*				temp;					
	//__shared__ 	dist_container*		dist_scatter;			
	//__shared__ 	dist_container*		dist_energy; 
	__shared__	spatial_data*		space;	
	__shared__	unsigned*			rxn;	
	__shared__	float*				E;		
	__shared__	float*				Q;		
	__shared__	unsigned*			rn_bank;
	__shared__	unsigned*			cellnum;
	__shared__	unsigned*			matnum;	
	__shared__	unsigned*			isonum;	
	__shared__	int*				talnum;
	//__shared__	unsigned*			yield;	
	__shared__	float*				weight;	
	__shared__	unsigned*			index;	
	// dynamically allocated shared array for material matrix
	extern __shared__	float		s_number_density_matrix[];


	// have thread 0 of block copy all pointers and static info into shared memory
	if (threadIdx.x == 0){
		n_isotopes					= d_xsdata[0].n_isotopes;								
		energy_grid_len				= d_xsdata[0].energy_grid_len;				
		total_reaction_channels		= d_xsdata[0].total_reaction_channels;
		rxn_numbers 				= d_xsdata[0].rxn_numbers;						
		rxn_numbers_total			= d_xsdata[0].rxn_numbers_total;					
		energy_grid 				= d_xsdata[0].energy_grid;						
		rxn_Q 						= d_xsdata[0].Q;												
		xs 							= d_xsdata[0].xs;												
		//awr 						= d_xsdata[0].awr;										
		//temp 						= d_xsdata[0].temp;										
		//dist_scatter 				= d_xsdata[0].dist_scatter;						
		//dist_energy 				= d_xsdata[0].dist_energy; 
		space						= d_particles[0].space;
		rxn							= d_particles[0].rxn;
		E							= d_particles[0].E;
		Q							= d_particles[0].Q;	
		rn_bank						= d_particles[0].rn_bank;
		cellnum						= d_particles[0].cellnum;
		matnum						= d_particles[0].matnum;
		isonum						= d_particles[0].isonum;
		talnum 						= d_particles[0].talnum;
		//yield						= d_particles[0].yield;
		weight						= d_particles[0].weight;
		index						= d_particles[0].index;
		// copy material matrix, hopefully small enough to fit!
		memcpy(s_number_density_matrix,d_number_density_matrix,n_isotopes*n_materials*sizeof(float));
	}

	// make sure shared loads happen before anything else
	__syncthreads();

	// return immediately if out of bounds
	int tid_in = threadIdx.x+blockIdx.x*blockDim.x; 
	if (tid_in >= N){return;}

	// return if terminated
	unsigned this_rxn = rxn[tid_in];
	if (this_rxn>=900){return;}
	else              {this_rxn = 0;}

	// remap
	int tid=d_remap[tid_in];

	// declare
	float		norm[3];
	float		samp_dist		= 0.0;
	float		diff			= 0.0;
	unsigned	this_tope		= 999999999;
	unsigned	array_dex		= 0;
	float		dotp			= 0.0;
	float		macro_t_total	= 0.0;
	const float	epsilon			= 2.0e-5;
	const float	push_value		= 2.0;
	float surf_minimum, xhat_new, yhat_new, zhat_new, this_Q;

	// load from arrays
	unsigned	this_mat		=  matnum[tid];
	int 		tally_index 	=  talnum[tid];
	unsigned	dex				=   index[tid];  
	unsigned	rn				= rn_bank[tid];
	float		this_E			=       E[tid];
	float		x				=   space[tid].x;
	float		y				=   space[tid].y;
	float		z				=   space[tid].z;
	float		xhat			=   space[tid].xhat;
	float		yhat			=   space[tid].yhat;
	float		zhat			=   space[tid].zhat;
	float		surf_dist		=   space[tid].surf_dist;
	unsigned	enforce_BC		=   space[tid].enforce_BC;  
	memcpy(     norm,               space[tid].norm,  3*sizeof(float) );

	//
	//
	//
	//  MACROSCOPIC SECTION
	//  -- find interacting isotope
	//
	//

	// compute some things
	unsigned 	n_columns 		= n_isotopes + total_reaction_channels;
	float 		e0, e1;

	// check
	if(this_mat>=n_materials){
		printf("MACRO - tid %u this_mat %u > n_materials %u!!!!!\n",tid,this_mat,n_materials);
		return;
	}

	if (dex>=4294967294){

		if(dex==4294967294){
			// out of bounds below data
			dex = 0;
			e0	= energy_grid[dex];
			e1	= 0.0;

		}
		else{
			// out of bounds above data
			dex = energy_grid_len-1;
			e0	= energy_grid[dex];
			e1	= 0.0;
		}

		// outside data, pass one array
		// compute the total macroscopic cross section for this material
		macro_t_total = sum_cross_section(	n_isotopes,
											e0, this_E,
											&s_number_density_matrix[this_mat],  
											&xs[ dex   *n_columns]				);
	
		// determine the isotope in the material for this cell
		this_tope = sample_cross_section(	n_isotopes, macro_t_total, get_rand(&rn),
											e0, this_E,
											&s_number_density_matrix[this_mat],  
											&xs[ dex   *n_columns]	 			);

	}
	else{

		// energy edges
		e0 = energy_grid[dex];
		e1 = energy_grid[dex+1];

		// inside the data, pass two arrays
		// compute the total macroscopic cross section for this material
		macro_t_total = sum_cross_section(	n_isotopes,
											e0, e1, this_E,  
											&s_number_density_matrix[this_mat],
											&xs[ dex   *n_columns],  
											&xs[(dex+1)*n_columns] 				);
	
		// determine the isotope in the material for this cell
		this_tope = sample_cross_section(	n_isotopes, macro_t_total, get_rand(&rn),
											e0, e1, this_E,
											&s_number_density_matrix[this_mat],  
											&xs[ dex   *n_columns],  
											&xs[(dex+1)*n_columns]					);

	}

	if (this_tope==n_isotopes) {printf("this_tope==n_isotopes, E %6.4E macro_t %6.4E\n",this_E,macro_t_total);}

	// compute the interaction length
	samp_dist = -logf(get_rand(&rn))/macro_t_total;

	// do surf/samp compare
	diff = surf_dist - samp_dist;

	// calculate epsilon projection onto neutron trajectory
	// dotp positive = neutron is inside the cell (normal points out, trajectory must be coming from inside)
	// dotp negative = neutron is outside the cell
	dotp = 	norm[0]*xhat + 
			norm[1]*yhat + 
			norm[2]*zhat;
	surf_minimum = 1.2 * epsilon / fabsf(dotp);
	
	// surface logic
	if( diff < surf_minimum ){  // need to make some decisions so epsilon is handled correctly
		// if not outer cell, neutron placement is too close to surface.  risk of next interaction not seeing the surface due to epsilon.
		// preserve if in this cell or next, but make sure neutron is at least an epsilon away from the surface.
		if (diff < 0.0){ // next cell, enforce BC or push through
			if (enforce_BC == 1){  // black BC
				//printf("Killing at black BC...\n");
				x += (surf_dist + 3.0*push_value*surf_minimum) * xhat;
				y += (surf_dist + 3.0*push_value*surf_minimum) * yhat;
				z += (surf_dist + 3.0*push_value*surf_minimum) * zhat;
				this_rxn  = 999;  // leaking is 999
				this_tope = 999999997;  
			}
			else if(enforce_BC == 2){  // specular reflection BC
				printf("Reflecting at specular BC NOT IMPLEMENTED YET!\n");
				// move epsilon off of surface
				x += ((surf_dist*xhat) + copysignf(1.0,dotp)*push_value*epsilon*norm[0]); 
				y += ((surf_dist*yhat) + copysignf(1.0,dotp)*push_value*epsilon*norm[1]);
				z += ((surf_dist*zhat) + copysignf(1.0,dotp)*push_value*epsilon*norm[2]);
				// calculate reflection
				xhat_new = -(2.0 * dotp * norm[0]) + xhat; 
				yhat_new = -(2.0 * dotp * norm[1]) + yhat; 
				zhat_new = -(2.0 * dotp * norm[2]) + zhat; 
				// flags
				this_rxn  = 801;  // reflection is 801 
				this_tope = 999999996;  
			}
			else{   // next cell, move to intersection point, then move *out* epsilon along surface normal
				//printf("Moving to next cell...\n");
				x += surf_dist*xhat + copysignf(1.0,dotp)*push_value*epsilon*norm[0];
				y += surf_dist*yhat + copysignf(1.0,dotp)*push_value*epsilon*norm[1];
				z += surf_dist*zhat + copysignf(1.0,dotp)*push_value*epsilon*norm[2];
				this_rxn  = 800;  // resampling is 800
				this_tope = 999999998;  
				}
			}
		else{   // this cell, move to intersection point, then move *in* epsilon along surface normal 
			//printf("In this cell...\n");
			x += surf_dist*xhat - copysignf(1.0,dotp)*push_value*epsilon*norm[0];
			y += surf_dist*yhat - copysignf(1.0,dotp)*push_value*epsilon*norm[1];
			z += surf_dist*zhat - copysignf(1.0,dotp)*push_value*epsilon*norm[2];
			this_rxn = 0;
		}
	}
	else{  // near side of minimum, can simply move the neutron
			x += samp_dist * xhat;
			y += samp_dist * yhat;
			z += samp_dist * zhat;
			this_rxn = 0;
	}

	// skip rest if leaked or resampled
	if (this_rxn==0){
	
		//
		//
		//
		//  TALLY SECTION
		//  -- score tally if valid
		//  
		//
	
		// only score if converged and if cell flagged for a tally. 
		if( tally_index>=0 & converged==1){
	
			// check
			tally_data	this_tally	= d_tally[tally_index];
			if( this_tally.cell != cellnum[tid]){
				printf("TALLY CELL AND REPORTED INDEX MISMATCH!!!  tid %u tally_index %i tally_cell %u cellnum[tid] %u\n",tid,tally_index,this_tally.cell,cellnum[tid]);
			}
	
			// load  
			float 		this_weight		= weight[tid];
			unsigned 	bin_index 		= 0;
			float 		E_max			= this_tally.E_max;
			float 		E_min			= this_tally.E_min;
			unsigned 	tally_length	= this_tally.length;
	
			// determine bin number
			bin_index = logf(this_E/E_min)/logf(E_max/E_min)*(tally_length);
	
			//score the bins atomicly, could be bad if many neutrons are in a single bin since this will serialize their operations
			atomicAdd(&this_tally.score[ bin_index], this_weight/macro_t_total);
			atomicAdd(&this_tally.square[bin_index], this_weight/(macro_t_total * macro_t_total));
			atomicInc(&this_tally.count[ bin_index], 4294967295);
	
		}
	
		//
		//
		//
		//  MICROSCOPIC SECTION
		//  -- find reaction type
		//  -- only reachable if collided!
		//
			
		unsigned	col_start	=	0;
		unsigned	col_end		=	0;
		unsigned	this_col	=	0;
		float		micro_t		=	0.0;
		
		// compute the index ranges 
		if(this_tope>=n_isotopes){
			printf("micro - ISOTOPE NUMBER FROM MACRO > NUMBER OF ISOTOPES!  n_isotopes %u tope %u\n",n_isotopes,this_tope);
			return;
		}
		else{
			col_start	=	n_isotopes + rxn_numbers_total[this_tope];
			col_end		=	n_isotopes + rxn_numbers_total[this_tope+1];
		}

		if (dex>=4294967294){

			if(dex==4294967294){
				// out of bounds below data
				dex = 0;
				e0	= energy_grid[dex];
				e1	= 0.0;
			}
			else{
				// out of bounds above data
				dex = energy_grid_len-1;
				e0	= energy_grid[dex];
				e1	= 0.0;
			}

			// compute the interpolated total microscopic cross section for this isotope.  Use non-multiplier function overload.  Remember that total xs is stored in the first n_isotopes of columns, then come the individual reaction cross sections...
			micro_t = sum_cross_section(	1,
											e0, this_E,
											&xs[ dex   *n_columns + this_tope]	);

			// determine the reaction/Q for this isotope, use non-multiplier function overload.  Returns index from col_start!
			this_col = col_start + sample_cross_section(	(col_end-col_start), micro_t, get_rand(&rn),
															e0, this_E,
															&xs[ dex   *n_columns + col_start]			);
	
		}
		else{
		
			// compute the interpolated total microscopic cross section for this isotope.  Use non-multiplier function overload.  Remember that total xs is stored in the first n_isotopes of columns, then come the individual reaction cross sections...
			micro_t = sum_cross_section(	1,
											e0, e1, this_E,
											&xs[ dex   *n_columns + this_tope],  
											&xs[(dex+1)*n_columns + this_tope] );
			
			// determine the reaction/Q for this isotope, use non-multiplier function overload.  Returns index from col_start!
			this_col = col_start + sample_cross_section(	(col_end-col_start), micro_t, get_rand(&rn),
															e0, e1, this_E,
															&xs[ dex   *n_columns + col_start],  
															&xs[(dex+1)*n_columns + col_start]			);
		}

		// the the numbers for this column
		this_rxn	=	rxn_numbers[this_col];
		this_Q		=	rxn_Q[      this_col];
		array_dex	=	dex*n_columns + this_col; 
		
		//printf("this_tope %u this_E %6.4E micro_t %6.4E col_start %u col_end %u \n",this_tope,this_E,micro_t,col_start,col_end);
		
		// errors?
		if(this_rxn == 999999999){ 
			printf("micro - REACTION NOT SAMPLED CORRECTLY! tope=%u E=%10.8E dex=%u rxn=%u\n",this_tope, this_E, dex, this_rxn); //most likely becasue rn1=1.0
		}
		if(this_rxn == 3 | this_rxn==4 | this_rxn ==5 | this_rxn ==10 | this_rxn ==27){
			printf("MT=%u!!!, changing to 1102...\n",this_rxn);
			this_rxn = 1102;
		}
	
	}
	
	//
	//
	//
	//  OUTPUT
	//  -- write output to arrays
	//
	//

	if(this_rxn==0){printf("rxn for tid_in %d / tid %d is still ZERO at end of macro_micro!\n", tid_in, tid);}

	rxn[    tid_in]			=	this_rxn;			// rxn is sorted WITH the remapping vector, i.e. its index does not need to be remapped
	Q[      tid]			=	this_Q; 			// put reaction Q value into particle Q value
	rn_bank[tid]			=	rn;
	index[  tid]			=	array_dex;			// write MT array index to dex instead of energy vector index
	isonum[ tid]			=	this_tope;
	space[  tid].x			=	x;
	space[  tid].y			=	y;
	space[  tid].z			=	z;
	if( enforce_BC==2 ){
		space[tid].xhat		=	xhat_new;			// write reflected directions for specular BC
		space[tid].yhat		=	yhat_new;
		space[tid].zhat		=	zhat_new;
	}

}

void macro_micro(unsigned NUM_THREADS, unsigned N, unsigned converged, unsigned n_materials, unsigned n_isotopes, unsigned n_tallies, cross_section_data* d_xsdata, particle_data* d_particles, tally_data* d_tally, unsigned* d_remap, float* d_number_density_matrix ){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	macro_micro_kernel <<< blks, NUM_THREADS, n_materials*n_isotopes*sizeof(float) >>> ( N, converged, n_materials, n_tallies, d_xsdata, d_particles, d_tally, d_remap, d_number_density_matrix);
	check_cuda(cudaThreadSynchronize());

}

