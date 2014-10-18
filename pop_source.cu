#include <cuda.h>
#include <stdio.h>
#include "datadef.h"
#include "wfloat3.h"
#include "LCRNG.cuh"

__global__ void pop_source_kernel(unsigned N, unsigned* isonum, unsigned* completed, unsigned* scanned, unsigned* remap, unsigned* yield, unsigned* done, unsigned* index, unsigned* rxn, source_point* space, float* E , unsigned* rn_bank, float**  energydata, source_point* space_out, float* E_out, float * awr_list){

	int tid = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid >= N){return;}

	// return if no yield
	if (yield[tid]==0){return;}

	//if(done[tid]){return;}

	// external data
	unsigned 		position 	= scanned[tid];
	unsigned 		this_tope 	= isonum[tid];
	unsigned 		this_yield 	= yield[tid];
	unsigned 		dex  		= index[tid];
	float 			this_E 		= E[tid];
	unsigned 		this_rxn 	= rxn[tid]; 
	unsigned 		rn 			= rn_bank[tid];
	float * 		this_array 	= energydata[dex];
	unsigned 		data_dex 	= 0;
	source_point 	this_space 	= space[tid];
	wfloat3 		hats_old(this_space.xhat,this_space.yhat,this_space.zhat);
	//float 			this_awr	= awr_list[this_tope];
	//printf("tid %drxn %u yield %u space % 6.4E % 6.4E % 6.4E\n",tid,this_rxn,this_yield,this_space.x,this_space.y,this_space.z);

	// internal data
	float 		Emin=1e-11;
	float 		Emax=20.0;
	unsigned 	k, n, offset, vlen, next_vlen, law;
	float 		sampled_E, phi, mu, rn1, rn2, last_E, next_E, e_start, E0, E1, Ek, next_e_start, next_e_end, last_e_start, last_e_end, diff;;
	float 		cdf0, cdf1, e0, e1, m, pdf0, pdf1, arg,x,y,z;
	const float  pi           =   3.14159265359 ;
	const float  m_n          =   1.00866491600 ; // u
	float 		speed_n 	=   sqrtf(2.0*this_E/m_n);

	//read in values
	offset = 5;
	//printf("rxn %u eptr %p\n",this_rxn,this_array);
	memcpy(&last_E,   	&this_array[0], sizeof(float));
	memcpy(&next_E,   	&this_array[1], sizeof(float));
	memcpy(&vlen,   	&this_array[2], sizeof(float));
	memcpy(&next_vlen,	&this_array[3], sizeof(float));
	memcpy(&law, 		&this_array[4], sizeof(float)); 
	float r = (this_E-last_E)/(next_E-last_E);
	last_e_start = this_array[ offset ];
	last_e_end   = this_array[ offset + vlen - 1 ];
	next_e_start = this_array[ offset + 3*vlen ];
	next_e_end   = this_array[ offset + 3*vlen + next_vlen - 1];
	for(k=0 ; k < this_yield ; k++ ){
		//get proper data index
		data_dex = completed[ position+k ];
		//printf("tid %u position %u k %u data_dex %u done %u (xyz) % 6.4E % 6.4E % 6.4E\n",tid,position,k,data_dex,done[data_dex],this_space.x,this_space.y,this_space.z);
		//make sure data is done
		//if(!done[data_dex]){printf("overwriting into active data!\n");}
		//copy in values
		rn1 = get_rand(&rn);
		rn2 = get_rand(&rn);
		//sample energy dist
		sampled_E = 0.0;
		if(  rn2 >= r ){   //sample last E
			diff = next_e_end - next_e_start;
			e_start = next_e_start;
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
			diff = next_e_end - next_e_start;
			e_start = next_e_start;
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
	
		// lin-lin interpolation
		m 	= (pdf1 - pdf0)/(e1-e0);
		arg = pdf0*pdf0 + 2.0 * m * (rn1-cdf0);
		if(arg<0){
			E0 = e0 + (e1-e0)/(cdf1-cdf0)*(rn1-cdf0);
		}
		else{
			E0 	= e0 + (  sqrtf( arg ) - pdf0) / m ;
		}
		// histogram interpolation
		//E0 = e0 + (rn1-cdf0)/pdf0;
		
		//scale it
		E1 = last_e_start + r*( next_e_start - last_e_start );
		Ek = last_e_end   + r*( next_e_end   - last_e_end   );
		sampled_E = E1 +(E0-e_start)*(Ek-E1)/diff;
		//sampled_E = E0;
		//printf("%6.4E\n",sampled_E);

		if(this_rxn==918){
			//sample isotropic directions
			rn1 = get_rand(&rn);
			rn2 = get_rand(&rn);
			mu  = 2.0*rn1-1.0; 
			phi = 2.0*pi*rn2;
			x = sqrtf(1.0-(mu*mu))*cosf(phi);
			y = sqrtf(1.0-(mu*mu))*sinf(phi);
			z = mu;
			//printf("% 6.4E % 6.4E % 6.4E\n",x,y,z);
		}
		else{  
			// pass old values from (n,2/3/4n) 
			x 	= this_space.xhat;
			y 	= this_space.yhat;
			z 	= this_space.zhat;	
		}
	
		//printf("tid %u k %u mu % 6.4E phi % 6.4E rn1 % 6.4E rn2 % 6.4E compactpos %u realpos %u\n",tid,k,mu,phi,rn1,rn2, position, completed[k+position]);
	
		//check data
		//printf("done? %u\n",done[ data_dex ]);
	
		//check limits
		if (sampled_E >= Emax){sampled_E = Emax * 0.99;}//printf("enforcing limits in pop data_dex=%u, sampled_E = %6.4E\n",data_dex,sampled_E);}
		if (sampled_E <= Emin){sampled_E = Emin * 1.01;}//printf("enforcing limits in pop data_dex=%u, sampled_E = %6.4E\n",data_dex,sampled_E);}


		// sync before writes
		__syncthreads();

		// set data
		space_out[ data_dex ].x 			= this_space.x;
		space_out[ data_dex ].y 			= this_space.y;
		space_out[ data_dex ].z 			= this_space.z;
		space_out[ data_dex ].xhat 			= x;
		space_out[ data_dex ].yhat 			= y;
		space_out[ data_dex ].zhat 			= z;
		space_out[ data_dex ].enforce_BC 	= 0;
		space_out[ data_dex ].surf_dist 	= 99999.0;
		space_out[ data_dex ].macro_t 		= 8.675309;
		rxn  	 [ data_dex ]				= this_rxn;
		index 	 [ data_dex ] 				= dex;
		isonum   [ data_dex ]  				= this_tope;
		rn_bank  [tid] = rn;

		//printf("rxn %u in pop dex %u\n",this_rxn,data_dex);
		if(this_rxn==918){
			E_out 	 [ data_dex ] 	= sampled_E;
		}
		else{ // pass (n,2/3/4n) to cscatter
			E_out 	 [ data_dex ] 	= this_E;
		}
		//printf("rxn %u dex %u E_out %6.4E data_dex %u\n",this_rxn,dex,E_out[data_dex],data_dex);
		//done [ data_dex ] 		= 0;
		//yield[ data_dex ] 		= 0;
		//printf("%u % 6.4E % 6.4E % 6.4E % 6.4E\n",data_dex,sampled_E,space_out[ data_dex ].x,space_out[ data_dex ].y,space_out[ data_dex ].z); 

	}


}

void pop_source( unsigned NUM_THREADS,  unsigned N, unsigned* isonum, unsigned* d_completed, unsigned* d_scanned, unsigned* d_remap, unsigned* d_yield, unsigned* d_done, unsigned* d_index, unsigned* d_rxn, source_point* d_space, float* d_E , unsigned* d_rn_bank, float ** energydata, source_point* space_out, float* E_out, float * awr_list){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	pop_source_kernel <<< blks, NUM_THREADS >>> ( N, isonum, d_completed, d_scanned, d_remap, d_yield, d_done, d_index, d_rxn, d_space, d_E , d_rn_bank, energydata, space_out, E_out, awr_list);
	cudaThreadSynchronize();

}

