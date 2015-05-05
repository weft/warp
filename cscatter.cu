#include <cuda.h>
#include <stdio.h>
#include "datadef.h"
#include "wfloat3.h"
#include "binary_search.h"
#include "LCRNG.cuh"

__global__ void cscatter_kernel(unsigned N, unsigned run_mode, unsigned starting_index, unsigned* remap, unsigned* isonum, unsigned * index, unsigned * rn_bank, float * E, source_point * space, unsigned * rxn, float * awr_list, float * Q, unsigned * done, float** scatterdat, float** energydat){


	int tid_in = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid_in >= N){return;}       //return if out of bounds
	
	//remap
	int tid, this_rxn;

	// check
	if(run_mode){ // regular transport
		tid=remap[starting_index + tid_in];
		this_rxn = rxn[starting_index + tid_in];
		//if (this_rxn != 91 ){printf("cscatter kernel accessing wrong reaction @ dex %u rxn %u\n",tid, this_rxn);return;}  //print and return if not continuum scatter
	}
	else{  // pop mode
		tid=tid_in;
		this_rxn = rxn[tid];
		if ( this_rxn != 916 & this_rxn != 917 & this_rxn != 937 & this_rxn != 924 & this_rxn != 922 & this_rxn != 928 & this_rxn != 924 & this_rxn !=932 & this_rxn != 933 & this_rxn != 941){return;}  // return if not multiplicity reaction
		
	}

	//constants
	const float  pi           =   3.14159265359 ;
	const float  m_n          =   1.00866491600 ; // u
	//const float  temp         =   0.025865214e-6;    // MeV
	const float  E_cutoff     =   1e-11;
	const float  E_max        =   20.0; //MeV
	// load history data
	unsigned 	this_tope 	= isonum[tid];
	unsigned 	this_dex	= index[tid];
	float 		this_E 		= E[tid];
	wfloat3 	hats_old(space[tid].xhat,space[tid].yhat,space[tid].zhat);
	float 		this_awr	= awr_list[this_tope];
	float * 	this_Sarray = scatterdat[this_dex];
	float * 	this_Earray =  energydat[this_dex];
	unsigned	rn			= rn_bank[ tid ];

	// check pointers
	if(this_Sarray == 0x0){
		printf("null pointer in cscatter!,dex %u rxn %u tope %u E %6.4E run mode %u\n",this_dex,this_rxn,this_tope,this_E,run_mode);
		return;
	}
	if(this_Earray == 0x0){
		printf("null pointer in cscatter!,dex %u rxn %u tope %u E %6.4E run mode %u\n",this_dex,this_rxn,this_tope,this_E,run_mode);
		return;
	}


	// internal kernel variables
	float 		mu, phi, next_E, last_E, sampled_E, e_start, E0, E1, Ek, next_e_start, next_e_end, last_e_start, last_e_end, diff;
    unsigned 	vlen, next_vlen, offset, n, law, intt; 
    unsigned  	isdone = 0;
	float  		speed_n          	=   sqrtf(2.0*this_E/m_n);
	float 		E_new				=   0.0;
	//float 		a 					= 	this_awr/(this_awr+1.0);
	wfloat3 	v_n_cm,v_t_cm,v_n_lf,v_t_lf,v_cm, hats_new, hats_target;
	float 		cdf0,e0,A,R,pdf0,rn1,cdf1,pdf1,e1;
	//float		cdf1,pdf1,e1;
	//float 		v_rel,E_rel;


	// make target isotropic (continuum is high energy, good approx)
	//mu  = (2.0*get_rand(&rn)) - 1.0;
	//phi = 2.0*pi*get_rand(&rn);
	//hats_target.x = sqrtf(1.0-(mu*mu))*cosf(phi);
	//hats_target.y = sqrtf(1.0-(mu*mu))*sinf(phi); 
	//hats_target.z = mu;
	
	// make speed vectors, assume high enough energy to approximate target as stationary
	v_n_lf = hats_old    * speed_n;
	v_t_lf = hats_target * 0.0;

	// calculate  v_cm
	v_cm = (v_n_lf + (v_t_lf*this_awr))/(1.0+this_awr);

	//transform neutron velocity into CM frame
	v_n_cm = v_n_lf - v_cm;
	v_t_cm = v_t_lf - v_cm;

	//
	//sample energy
	//
	//read in values
	offset = 6;
	memcpy(&last_E,   	&this_Earray[0], sizeof(float));
	memcpy(&next_E,   	&this_Earray[1], sizeof(float));
	memcpy(&vlen,   	&this_Earray[2], sizeof(float));
	memcpy(&next_vlen,	&this_Earray[3], sizeof(float));
	memcpy(&law, 		&this_Earray[4], sizeof(float));
	memcpy(&intt, 		&this_Earray[5], sizeof(float));


	//LAW = 44:
	if (law==44){
		if(run_mode==0 & this_E<last_E){
			this_E = last_E;
		}

		float r = (this_E-last_E)/(next_E-last_E);
		if(r<0){
			printf("DATA NOT WITHIN ENERGY INTERVAL tid %u r % 10.8E rxn %u isotope %u this_E % 10.8E last_E % 10.8E next_E % 10.8E dex %u\n",tid,r,this_rxn,this_tope,this_E,last_E,next_E,this_dex);
		}	

		last_e_start = this_Earray[ offset ];
		last_e_end   = this_Earray[ offset + vlen - 1 ];
		next_e_start = this_Earray[ offset + 3*vlen ];
		next_e_end   = this_Earray[ offset + 3*vlen + next_vlen - 1];
		//unsigned svlen,next_svlen,len;
		//printf("law=%u vlen/next= %u %u, E-last/this/next= %6.4E %6.4E %6.4E\n",law,vlen,next_vlen,last_E,this_E,next_E);
		sampled_E = 0.0;
		rn1 = get_rand(&rn);
		if(  get_rand(&rn) >= r ){   //sample last E
			diff = last_e_end - last_e_start;
			e_start = last_e_start;
			//n = binary_search( &this_Earray[ offset + vlen ] , rn1, vlen);
			for ( n=0 ; n<vlen-1 ; n++ ){
				cdf0 		= this_Earray[ (offset +   vlen ) + n+0];
				cdf1 		= this_Earray[ (offset +   vlen ) + n+1];
				if( rn1 >= cdf0 & rn1 < cdf1 ){
					break;
				}
			}
			pdf0		= this_Earray[ (offset + 2*vlen ) + n+0];
			pdf1		= this_Earray[ (offset + 2*vlen ) + n+1];
			e0  		= this_Earray[ (offset          ) + n+0];
			e1  		= this_Earray[ (offset          ) + n+1]; 
			offset = 4;
			A = this_Sarray[ (offset)      + n ];
			R = this_Sarray[ (offset+vlen) + n ];
		}
		else{
			diff = next_e_end - next_e_start;
			e_start = next_e_start;
			//n = binary_search( &this_Earray[ offset + 3*vlen + next_vlen] , rn1, next_vlen);
			for ( n=0 ; n<next_vlen-1 ; n++ ){
				cdf0 		= this_Earray[ (offset + 3*vlen +   next_vlen ) + n+0];
				cdf1  		= this_Earray[ (offset + 3*vlen +   next_vlen ) + n+1];
				if( rn1 >= cdf0 & rn1 < cdf1 ){
					break;
				}
			}
			pdf0		= this_Earray[ (offset + 3*vlen + 2*next_vlen ) + n+0];
			pdf1		= this_Earray[ (offset + 3*vlen + 2*next_vlen ) + n+1];
			e0   		= this_Earray[ (offset + 3*vlen               ) + n+0];
			e1   		= this_Earray[ (offset + 3*vlen               ) + n+1];
			offset = 4;
			A = this_Sarray[ (offset+2*vlen)           +n  ] ;
			R = this_Sarray[ (offset+2*vlen+next_vlen) +n  ];
		}
	
		if (intt==1){
		// histogram interpolation
			E0 = e0 + (rn1-cdf0)/pdf0;
		}
		else if(intt==2){
		////lin-lin interpolation
		//	float m   = (pdf1 - pdf0)/(e1-e0);
		//	float arg = pdf0*pdf0 + 2.0 * m * (rn1-cdf0);
		//	if(arg<0){
		//		E0 = e0 + (e1-e0)/(cdf1-cdf0)*(rn1-cdf0);
		//	}
		}
		else{
		//	E0 	= e0 + (  sqrtf( arg ) - pdf0) / m ;
		}
	
		//scale it
		E1 = last_e_start + r*( next_e_start - last_e_start );
		Ek = last_e_end   + r*( next_e_end   - last_e_end   );
		sampled_E = E1 +(E0-e_start)*(Ek-E1)/diff;
		//printf("sampled_E % 10.8E E1 % 10.8E Ek % 10.8E E0 % 10.8E e_start % 10.8E diff % 10.8E last_e_end % 10.8E r % 10.8E next_e_end % 10.8E last_e_end % 10.8E \n",sampled_E,E1,Ek,E0,e_start,diff,last_e_end,r,next_e_end,last_e_end);
		//sampled_E = E0;
	
		// find mu
		rn1 = get_rand(&rn);
		if(get_rand(&rn)>R){
			float T = (2.0*rn1-1.0)*sinhf(A);
			mu = logf(T+sqrtf(T*T+1.0))/A;
		}
		else{
			mu = logf(rn1*expf(A)+(1.0-rn1)*expf(-A))/A;
		}
	}
	else if (law==61){

	}

	// rotate direction vector
	hats_old = v_n_cm / v_n_cm.norm2();
	hats_old = hats_old.rotate(mu, get_rand(&rn));

	//  scale to sampled energy
	v_n_cm = hats_old * sqrtf(2.0*sampled_E/m_n);
	
	// transform back to L
	v_n_lf = v_n_cm + v_cm;
	hats_new = v_n_lf / v_n_lf.norm2();
	hats_new = hats_new / hats_new.norm2(); // get higher precision, make SURE vector is length one
	
	// calculate energy in lab frame
	E_new = 0.5 * m_n * v_n_lf.dot(v_n_lf);

	// enforce limits
	if ( E_new <= E_cutoff | E_new > E_max ){
		isdone=1;
		this_rxn = 998;  // ecutoff code
	}
	
	//if(this_rxn==91){printf("%u % 6.4E %6.4E %6.4E %6.4E %u %u\n",this_rxn,mu,sampled_E,this_E,E_new, vlen, next_vlen);}
	//if(this_rxn==91){printf("%6.4E %6.4E %6.4E\n",E_new,this_E,E_new/this_E);}
	//printf("n,vlen %u %u S,Eptrs %p %p Enew,samp %6.4E %6.4E A,R %6.4E %6.4E\n",n,vlen,this_Sarray,this_Earray,E_new,sampled_E,A,R);
	//printf("%u dex %u sampled_E % 10.8E norm2_lf % 10.8E mu % 10.8E\n",tid,this_dex,sampled_E,v_n_lf.norm2(),mu);

	// write reaction results if in trasport mode
	if(run_mode){
		done[tid]       			= isdone;
		rxn[starting_index+tid_in] 	= this_rxn;
	}
	
	// write universal results
	E[tid]          			= E_new;
	space[tid].xhat 			= hats_new.x;
	space[tid].yhat 			= hats_new.y;
	space[tid].zhat 			= hats_new.z;
	rn_bank[tid] 				= rn;	

}

void cscatter( cudaStream_t stream, unsigned NUM_THREADS, unsigned run_mode, unsigned N, unsigned starting_index, unsigned* remap, unsigned* isonum, unsigned * index, unsigned * rn_bank, float * E, source_point * space ,unsigned * rxn, float* awr_list, float * Q, unsigned* done, float** scatterdat, float** energydat){

	if(N<1){return;}
	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;
	
	cscatter_kernel <<< blks, NUM_THREADS , 0 , stream >>> (  N, run_mode, starting_index, remap, isonum, index, rn_bank, E, space, rxn, awr_list, Q, done, scatterdat, energydat);
	cudaThreadSynchronize();

}

