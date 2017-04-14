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
		n_isotopes			= d_xsdata[0].n_isotopes;			
		energy_grid_len			= d_xsdata[0].energy_grid_len;				
		total_reaction_channels		= d_xsdata[0].total_reaction_channels;
		rxn_numbers 			= d_xsdata[0].rxn_numbers;			
		rxn_numbers_total		= d_xsdata[0].rxn_numbers_total;	
		energy_grid 			= d_xsdata[0].energy_grid;			
		rxn_Q 				= d_xsdata[0].Q;			
		xs 				= d_xsdata[0].xs;					
		//awr 				= d_xsdata[0].awr;					
		//temp 				= d_xsdata[0].temp;				
		//dist_scatter 			= d_xsdata[0].dist_scatter;				
		//dist_energy 			= d_xsdata[0].dist_energy; 
		space				= d_particles[0].space;
		rxn				= d_particles[0].rxn;
		E				= d_particles[0].E;
		Q				= d_particles[0].Q;	
		rn_bank				= d_particles[0].rn_bank;
		cellnum				= d_particles[0].cellnum;
		matnum				= d_particles[0].matnum;
		isonum				= d_particles[0].isonum;
		talnum 				= d_particles[0].talnum;
		//yield				= d_particles[0].yield;
		weight				= d_particles[0].weight;
		index				= d_particles[0].index;
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
	float		plane_vec[3];
	float		samp_dist	= 0.0;
	float		diff		= 0.0;
	unsigned	this_tope	= 999999999;
	unsigned	array_dex	= 0;
	unsigned	adj_dex		= 0;
	float		dotp		= 0.0;
	float		macro_t_total	= 0.0;
	const float	epsilon		= 5.0e-6;
	const float	push_value	= 2.0;
	float surf_minimum, this_Q;
	//float xhat_new, yhat_new, zhat_new;
	unsigned	cell_local[10];
	unsigned	mat_local[10];
	float		dist_local[10];
	float		macro_maj	= 0.0;
	float		x_new = 0.0;
	float		y_new = 0.0;
	float		z_new = 0.0;

	// load from arrays
	unsigned	this_mat		=  matnum[tid];
	unsigned	this_cell		=  cellnum[tid];
	int 		tally_index 		=  talnum[tid];
	unsigned	dex			=   index[tid];  
	unsigned	rn			= rn_bank[tid];
	float		this_E			=       E[tid];
	float		x			=   space[tid].x;
	float		y			=   space[tid].y;
	float		z			=   space[tid].z;
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

	for(int i = 0; i < 10; i++)
	{
		cell_local[i] = space[tid].cell[i];
		mat_local[i] = space[tid].mat[i];
		dist_local[i] = space[tid].dist[i];
	}

	for(int i = 1; i < 10; i++)
	{
		if(dist_local[i] > 0){dist_local[i] += dist_local[i-1];}
	}

	for(int i = 0; i < 9; i++)
	{
		if(cell_local[i] == cell_local[i+1] && cell_local[i] != -1)
		{
			for(int j = i+1; j > 0; j--)
			{
				cell_local[j] = cell_local[j-1];
				mat_local[j] = mat_local[j-1];
			}
		}
	}

	cell_local[0] = this_cell;
	mat_local[0] = this_mat;

	// compute some things
	unsigned	n_columns 		= n_isotopes + total_reaction_channels;
	float 		e0, e1, rnd, rnd2, rnd3;
	int found = 0;

	// check
	if(this_mat>=n_materials){
		printf("MACRO - tid %u this_mat %u > n_materials %u!!!!!\n",tid,this_mat,n_materials);
		return;
	}

	rnd = get_rand(&rn);
	rnd2 = get_rand(&rn);

	if (dex>=4294967294){

		if(dex==4294967294){
			// out of bounds below data
			adj_dex = 0;
			e0	= energy_grid[adj_dex];

		}
		else{
			// out of bounds above data
			adj_dex = energy_grid_len-1;
			e0	= energy_grid[adj_dex];
		}

		for(int i = 0; i < n_materials; i++)
		{
			macro_t_total = sum_cross_section(n_isotopes, e0, this_E,
				&s_number_density_matrix[i*n_isotopes],	&xs[adj_dex*n_columns]);
			if(macro_t_total > macro_maj){macro_maj = macro_t_total;}
		}

		// compute the interaction length
		samp_dist = -logf(rnd)/macro_maj;
		x_new = x + (samp_dist * xhat);
		y_new = y + (samp_dist * yhat);
		z_new = z + (samp_dist * zhat);
		found = 0; 

		for(int i = 0; i < 10; i++)
		{
			if(samp_dist < dist_local[i] && cell_local[i] != -1)
			{
				this_cell = cell_local[i];
				this_mat = mat_local[i];
				found = 1;
				break;
			}
		}

		// outside data, pass one array
		// compute the total macroscopic cross section for this material
		macro_t_total = sum_cross_section(n_isotopes, e0, this_E,
				&s_number_density_matrix[this_mat*n_isotopes],	&xs[adj_dex*n_columns]);

		if(!isfinite(macro_t_total) | macro_t_total<=0.0){printf("1 macro_t_total is wrong:  % 6.4E e0 % 6.4E e1 % 6.4E \n",macro_t_total,e0,e1);}
	
		// determine the isotope in the material for this cell
		this_tope = sample_cross_section(n_isotopes, macro_t_total, rnd2, e0, this_E,
			    &s_number_density_matrix[this_mat*n_isotopes], &xs[adj_dex*n_columns]);

	}
	else{

		// energy edges
		e0 = energy_grid[dex];
		e1 = energy_grid[dex+1];

		for(int i = 0; i < n_materials; i++)
		{
			macro_t_total = sum_cross_section(n_isotopes, e0, e1, this_E,  
				&s_number_density_matrix[i*n_isotopes], &xs[dex*n_columns],
				&xs[(dex+1)*n_columns]);
			if(macro_t_total > macro_maj){macro_maj = macro_t_total;}
		}

		// compute the interaction length
		samp_dist = -logf(rnd)/macro_maj;
		x_new = x + (samp_dist * xhat);
		y_new = y + (samp_dist * yhat);
		z_new = z + (samp_dist * zhat);
		found = 0; 

		for(int i = 0; i < 10; i++)
		{
			if(samp_dist < dist_local[i] && cell_local[i] != -1)
			{
				this_cell = cell_local[i];
				this_mat = mat_local[i];
				found = 1;
				break;
			}
		}

		// inside the data, pass two arrays
		// compute the total macroscopic cross section for this material
		macro_t_total = sum_cross_section(n_isotopes, e0, e1, this_E,  
				&s_number_density_matrix[this_mat*n_isotopes], &xs[dex*n_columns],
				&xs[(dex+1)*n_columns]);

		if(!isfinite(macro_t_total) | macro_t_total<=0.0){printf("2 macro_t_total is wrong:  % 6.4E e0 % 6.4E e1 % 6.4E \n",macro_t_total,e0,e1);}
	
		// determine the isotope in the material for this cell
		this_tope = sample_cross_section(n_isotopes, macro_t_total, rnd2, e0, e1, this_E,
			    &s_number_density_matrix[this_mat*n_isotopes], &xs[dex*n_columns],  
			    &xs[(dex+1)*n_columns]);

	}

	if (this_tope==n_isotopes) {printf("this_tope==n_isotopes, tid %d E %6.4E macro_t_total %6.4E rn %12.10E dex %u\n",tid,this_E,macro_t_total,rnd2,dex);}

	// calculate epsilon projection onto neutron trajectory
	// dotp positive = neutron is inside the cell (normal points out, trajectory must be coming from inside)
	// dotp negative = neutron is outside the cell
	dotp = 	norm[0]*xhat + norm[1]*yhat + norm[2]*zhat;
	surf_minimum = 1.2 * epsilon / fabsf(dotp);

	rnd3 = get_rand(&rn);
	
	// surface logic
	if(found){ 
		if(rnd3 > macro_t_total/macro_maj){this_rxn = 800; this_tope = 999999998;}
		else{this_rxn = 0;}
	}
	else{
//	if( diff < surf_minimum ){  // need to make some decisions so epsilon is handled correctly
		// if not outer cell, neutron placement is too close to surface.  risk of next interaction not seeing the surface due to epsilon.
		// preserve if in this cell or next, but make sure neutron is at least an epsilon away from the surface.
		if (enforce_BC == 1){  // black BC
			//printf("Killing at black BC...\n");
			this_rxn  = 999;  // leaking is 999
			this_tope = 999999997;  
		}
		else if(enforce_BC == 2){  // specular reflection BC
			// calculate plane vector and normalize to one
			plane_vec[0] = xhat - dotp*norm[0];
			plane_vec[1] = yhat - dotp*norm[1];
			plane_vec[2] = zhat - dotp*norm[2];
			plane_vec[0] = plane_vec[0] / sqrtf(plane_vec[0]*plane_vec[0]+plane_vec[1]*plane_vec[1]+plane_vec[2]*plane_vec[2]);
			plane_vec[1] = plane_vec[1] / sqrtf(plane_vec[0]*plane_vec[0]+plane_vec[1]*plane_vec[1]+plane_vec[2]*plane_vec[2]);
			plane_vec[2] = plane_vec[2] / sqrtf(plane_vec[0]*plane_vec[0]+plane_vec[1]*plane_vec[1]+plane_vec[2]*plane_vec[2]);
			// first move intersection point back in the incoming direction in case close to a wall
			x += surf_dist*xhat - push_value*epsilon*plane_vec[0];
			y += surf_dist*yhat - push_value*epsilon*plane_vec[1];
			z += surf_dist*zhat - push_value*epsilon*plane_vec[2];
			// move epsilon off of surface
			x += - copysignf(1.0,dotp)*push_value*epsilon*norm[0]; 
			y += - copysignf(1.0,dotp)*push_value*epsilon*norm[1];
			z += - copysignf(1.0,dotp)*push_value*epsilon*norm[2];
			// calculate reflection
			xhat += -(2.0 * dotp * norm[0]);
			yhat += -(2.0 * dotp * norm[1]);   
			zhat += -(2.0 * dotp * norm[2]);   
			// flags
			this_rxn  = 801;  // reflection is 801 
			this_tope = 999999996;  
		}
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
			atomicAdd(&this_tally.square[bin_index], this_weight*this_weight/(macro_t_total * macro_t_total));
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
		//float		micro_rxn	=	0.0;
		
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
				adj_dex = 0;
				e0		= energy_grid[adj_dex];
			}
			else{
				// out of bounds above data
				adj_dex = energy_grid_len-1;
				e0		= energy_grid[adj_dex];
			}

			// compute the interpolated total microscopic cross section for this isotope.  Use non-multiplier function overload.  Remember that total xs is stored in the first n_isotopes of columns, then come the individual reaction cross sections...
			micro_t = sum_cross_section(1,e0, this_E,&xs[ adj_dex *n_columns + this_tope]);

			if(!isfinite(micro_t) | micro_t<=0.0){printf("1 micro_t is wrong:  % 6.4E e0 % 6.4E e1 % 6.4E \n!",micro_t,e0,e1);}

			// determine the reaction/Q for this isotope, use non-multiplier function overload.  Returns index from col_start!
			this_col = col_start + sample_cross_section((col_end-col_start), micro_t, 
					get_rand(&rn),e0, this_E,&xs[ adj_dex   *n_columns + col_start]);

			// compute array index
			array_dex	=	adj_dex*n_columns + this_col;

			// get interpolated cross section for this rxn (debug)
			//micro_rxn = sum_cross_section(	1,
			//								e0, e1, this_E,
			//								&xs[ adj_dex   *n_columns + this_col],  
			//								&xs[(adj_dex+1)*n_columns + this_col] );

	
		}
		else{

			// energy edges
			e0 = energy_grid[dex];
			e1 = energy_grid[dex+1];
		
			// compute the interpolated total microscopic cross section for this isotope.  Use non-multiplier function overload.  Remember that total xs is stored in the first n_isotopes of columns, then come the individual reaction cross sections...
			micro_t = sum_cross_section(1,e0, e1, this_E,&xs[ dex   *n_columns + this_tope],  
								&xs[(dex+1)*n_columns + this_tope] );

			if(!isfinite(micro_t) | micro_t<=0.0){printf("2 micro_t is wrong:  % 6.4E e0 % 6.4E e1 % 6.4E \n",micro_t,e0,e1);}
			
			// determine the reaction/Q for this isotope, use non-multiplier function overload.  Returns index from col_start!
			this_col = col_start + sample_cross_section((col_end-col_start), micro_t, 
				get_rand(&rn),e0, e1, this_E,&xs[ dex   *n_columns + col_start],  
				&xs[(dex+1)*n_columns + col_start]);

			// compute array index
			array_dex	=	dex*n_columns + this_col;

			// get interpolated cross section for this rxn (debug)
			//micro_rxn = sum_cross_section(	1,
			//								e0, e1, this_E,
			//								&xs[ dex   *n_columns + this_col],  
			//								&xs[(dex+1)*n_columns + this_col] );

		}

		// the the numbers for this column
		this_rxn	=	rxn_numbers[this_col];
		this_Q		=	rxn_Q[      this_col]; 
		
		//printf("this_tope %u this_E %6.4E micro_t %6.4E col_start %u col_end %u \n",this_tope,this_E,micro_t,col_start,col_end);
		//printf("(rxn,E,micro_t,micro_rxn) %4u % 6.4E % 6.4E % 6.4E \n",this_rxn,this_E,micro_t,micro_rxn);

		// errors?
		if(this_rxn == 999999999){ 
			printf("micro - REACTION NOT SAMPLED CORRECTLY! tope=%u E=%10.8E dex=%u rxn=%u\n",this_tope, this_E, dex, this_rxn); //most likely becasue rnd=1.0
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
	if(this_rxn==4){printf("rxn for tid_in %d / tid %d is 4 at end of macro_micro!\n", tid_in, tid);}
	//if(this_rxn==999){printf("(macro_t_total,surf_dist,samp_dist,surf_minimum,dotp) % 6.4E % 6.4E % 6.4E % 6.4E % 6.4E\n",macro_t_total,surf_dist,samp_dist,surf_minimum,dotp);}

	rxn[    tid_in]		=	this_rxn;			// rxn is sorted WITH the remapping vector, i.e. its index does not need to be remapped
	Q[      tid]		=	this_Q; 			// put reaction Q value into particle Q value
	rn_bank[tid]		=	rn;
	index[  tid]		=	array_dex;			// write MT array index to dex instead of energy vector index
	isonum[ tid]		=	this_tope;
	space[  tid].x		=	x_new;
	space[  tid].y		=	y_new;
	space[  tid].z		=	z_new;
	space[  tid].xhat	=	xhat;
	space[  tid].yhat	=	yhat;
	space[  tid].zhat	=	zhat;
	cellnum[tid] = this_cell;
	matnum[tid] = this_mat;

	//if( enforce_BC==2 ){
	//	printf("(norm) % 6.4E % 6.4E % 6.4E (xyz-hat) % 6.4E % 6.4E % 6.4E (xyz-hat_new) % 6.4E % 6.4E % 6.4E\n",norm[0],norm[1],norm[2],xhat,yhat,zhat,xhat_new,yhat_new,zhat_new);
	//	space[tid].xhat		=	xhat_new;			// write reflected directions for specular BC
	//	space[tid].yhat		=	yhat_new;
	//	space[tid].zhat		=	zhat_new;
	//}

}


/**
 * \brief a
 * \details b
 *
 * @param[in]    NUM_THREADS              - the number of threads to run per thread block
 * @param[in]    N                        - the total number of threads to launch on the grid
 * @param[in]    converged                - flag for tally scoring
 * @param[in]    n_materials              - number of materials
 * @param[in]    n_isotopes               - number of isotopes
 * @param[in]    n_tallies                - number of tallies
 * @param[in]    d_xsdata                 - device pointer to cross section data pointer array 
 * @param[in]    d_particles              - device pointer to particle data pointer array 
 * @param[in]    d_tally                  - device pointer to tally array
 * @param[in]    d_remap                  - device pointer to data remapping vector
 * @param[in]    d_number_density_matrix  - device pointer to material number density array
 */ 
void macro_micro(unsigned NUM_THREADS, unsigned N, unsigned converged, unsigned n_materials, unsigned n_isotopes, unsigned n_tallies, cross_section_data* d_xsdata, particle_data* d_particles, tally_data* d_tally, unsigned* d_remap, float* d_number_density_matrix ){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	macro_micro_kernel <<< blks, NUM_THREADS, n_materials*n_isotopes*sizeof(float) >>> ( N, converged, n_materials, n_tallies, d_xsdata, d_particles, d_tally, d_remap, d_number_density_matrix);
	check_cuda(cudaThreadSynchronize());

}

