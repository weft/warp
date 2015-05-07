#include <cuda.h>
#include <stdio.h>
#include "datadef.h"
#include "wfloat3.h"
#include "LCRNG.cuh"
#include "binary_search.h"

__device__ void process_fission(unsigned this_yield, unsigned* rn, unsigned position, unsigned this_tope, unsigned this_awr, float this_E, source_point this_space, float* this_Earray, float* this_Sarray, source_point* space_out, float* E_out){

	// internal data
	unsigned 	k, n, offset, vlen, next_vlen, law, data_dex, intt;
	float 		sampled_E, phi, mu, rn1, rn2, last_E, next_E, e_start, E0, E1, Ek, next_e_start, next_e_end, last_e_start, last_e_end, diff, r;
	float 		cdf0, cdf1, e0, e1, m, pdf0, pdf1, arg,x,y,z,A, R;
	const float pi 		= 3.14159265359;
	const float Emin 	= 1e-11;
	const float Emax 	= 20.0;

	//read in values
	offset = 6;
	memcpy(&last_E,   	&this_Earray[0], sizeof(float));
	memcpy(&next_E,   	&this_Earray[1], sizeof(float));
	memcpy(&vlen,   	&this_Earray[2], sizeof(float));
	memcpy(&next_vlen,	&this_Earray[3], sizeof(float));
	memcpy(&law, 		&this_Earray[4], sizeof(float)); 
	memcpy(&intt, 		&this_Earray[5], sizeof(float)); 


	if (law == 4){
		r = (this_E-last_E)/(next_E-last_E);
		last_e_start = this_Earray[ offset ];
		last_e_end   = this_Earray[ offset + vlen - 1 ];
		next_e_start = this_Earray[ offset + 3*vlen ];
		next_e_end   = this_Earray[ offset + 3*vlen + next_vlen - 1];
	}

	// loop over the for the number of (rebased) yielded particles
	for(k=0 ; k < this_yield ; k++ ){

		//get proper data index
		data_dex = position+k ;
		
		if (law==4){ // tabular
			rn1 = get_rand(rn);
			rn2 = get_rand(rn);
	
			//sample energy dist
			sampled_E = 0.0;
			if(  rn2 >= r ){   //sample last E
				diff = next_e_end - next_e_start;
				e_start = next_e_start;
				for ( n=0 ; n<vlen-1 ; n++ ){
					cdf0 		= this_Earray[ (offset +   vlen ) + n+0];
					cdf1 		= this_Earray[ (offset +   vlen ) + n+1];
					pdf0		= this_Earray[ (offset + 2*vlen ) + n+0];
					pdf1		= this_Earray[ (offset + 2*vlen ) + n+1];
					e0  		= this_Earray[ (offset          ) + n+0];
					e1  		= this_Earray[ (offset          ) + n+1]; 
					if( rn1 >= cdf0 & rn1 < cdf1 ){
						break;
					}
				}
			}
			else{
				diff = next_e_end - next_e_start;
				e_start = next_e_start;
				for ( n=0 ; n<next_vlen-1 ; n++ ){
					cdf0 		= this_Earray[ (offset + 3*vlen +   next_vlen ) + n+0];
					cdf1  		= this_Earray[ (offset + 3*vlen +   next_vlen ) + n+1];
					pdf0		= this_Earray[ (offset + 3*vlen + 2*next_vlen ) + n+0];
					pdf1		= this_Earray[ (offset + 3*vlen + 2*next_vlen ) + n+1];
					e0   		= this_Earray[ (offset + 3*vlen               ) + n+0];
					e1   		= this_Earray[ (offset + 3*vlen               ) + n+1];
					if( rn1 >= cdf0 & rn1 < cdf1 ){
						break;
					}
				}
			}
		
			if (intt==2){// lin-lin interpolation
				m 	= (pdf1 - pdf0)/(e1-e0);
				arg = pdf0*pdf0 + 2.0 * m * (rn1-cdf0);
				if(arg<0){
					E0 = e0 + (e1-e0)/(cdf1-cdf0)*(rn1-cdf0);
				}
				else{
					E0 	= e0 + (  sqrtf( arg ) - pdf0) / m ;
				}
			}
			else if(intt==1){// histogram interpolation
				E0 = e0 + (rn1-cdf0)/pdf0;
			}
			
			//scale it
			E1 = last_e_start + r*( next_e_start - last_e_start );
			Ek = last_e_end   + r*( next_e_end   - last_e_end   );
			sampled_E = E1 +(E0-e_start)*(Ek-E1)/diff;

		}
		else if (law==44){
			if(this_Sarray == 0x0){  //Sarray value is nu for 918, might have to put check in here
				printf("null pointer in pop multiplicity Sarray!, tope %u E %6.4E\n",this_tope,this_E);
				return;
			}
			//printf("%6.4E %6.4E %u %u %u ... %6.4E %6.4E ... %6.4E %6.4E %6.4E\n",this_Earray[0],this_Earray[1],vlen,next_vlen,law,this_Earray[5],this_Earray[6],this_Earray[offset],this_Earray[offset+vlen],this_Earray[offset+vlen+1]);
			if(this_E<last_E){this_E=last_E;}
			float r = (this_E-last_E)/(next_E-last_E);
			if(r<0){
				printf("r less than zero in source pop for multiplicity, r % 10.8E isotope %u this_E % 10.8E last_E % 10.8E next_E % 10.8E\n",r,this_tope,this_E,last_E,next_E);
			}
			last_e_start = this_Earray[ offset ];
			last_e_end   = this_Earray[ offset + vlen - 1 ];
			next_e_start = this_Earray[ offset + 3*vlen ];
			next_e_end   = this_Earray[ offset + 3*vlen + next_vlen - 1];
	
			//sample energy dist
			sampled_E = 0.0;
			rn1 = get_rand(rn);
			if(  get_rand(rn) >= r ){   //sample last E
				diff = last_e_end - last_e_start;
				e_start = last_e_start;
				//n = binary_search( &this_Earray[ offset + vlen ] , rn1, vlen);
				for ( n=0 ; n<vlen-1 ; n++ ){
					cdf0 		= this_Earray[ (offset +   vlen ) + n+0];
					cdf1 		= this_Earray[ (offset +   vlen ) + n+1];
					pdf0		= this_Earray[ (offset + 2*vlen ) + n+0];
					pdf1		= this_Earray[ (offset + 2*vlen ) + n+1];
					e0  		= this_Earray[ (offset          ) + n+0];
					e1  		= this_Earray[ (offset          ) + n+1]; 
					if( rn1 >= cdf0 & rn1 < cdf1 ){
						break;
					}
				}
				offset = 6;
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
					pdf0		= this_Earray[ (offset + 3*vlen + 2*next_vlen ) + n+0];
					pdf1		= this_Earray[ (offset + 3*vlen + 2*next_vlen ) + n+1];
					e0   		= this_Earray[ (offset + 3*vlen               ) + n+0];
					e1   		= this_Earray[ (offset + 3*vlen               ) + n+1];
					if( rn1 >= cdf0 & rn1 < cdf1 ){
						break;
					}
				}
				offset = 6;
				A = this_Sarray[ (offset+3*vlen)           +n  ] ;
				R = this_Sarray[ (offset+3*vlen+next_vlen) +n  ];
			}
			
			// histogram interpolation, intt=1
			E0 = e0 + (rn1-cdf0)/pdf0;
			
			//scale it
			E1 = last_e_start + r*( next_e_start - last_e_start );
			Ek = last_e_end   + r*( next_e_end   - last_e_end   );
			sampled_E = E1 +(E0-e_start)*(Ek-E1)/diff;
			
			// find mu
			rn1 = get_rand(rn);
			if(get_rand(rn)>R){
				float T = (2.0*rn1-1.0)*sinhf(A);
				mu = logf(T+sqrtf(T*T+1.0))/A;
			}
			else{
					mu = logf(rn1*expf(A)+(1.0-rn1)*expf(-A))/A;
				}
		}
		else if(law==7){   // maxwellian fission

			// get tabulated temperature
			float t0 = this_Earray[ offset     ];
			float t1 = this_Earray[ offset + 1 ];
			float U  = this_Earray[ offset + vlen       ];
			      e0 = this_Earray[ offset + vlen*2     ];
			      e1 = this_Earray[ offset + vlen*2 + 1 ];
			float  T = 0.0;

			// interpolate T
			if (e1==e0){  // in top bin, both values are the same
				T = t0;
			}
			else if (intt==2){// lin-lin interpolation
				m = (this_E - e0)/(e1 - e0);
                T = (1.0 - m)*t0 + m*t1;
			}
			else if(intt==1){// histogram interpolation
				T  = (t1 - t0)/(e1 - e0) * this_E + t0;
			}
			else{
				printf("something is wrong in pop for law 7");
			}

			// rejection sample
			sampled_E = 9999999.9;
			while ( sampled_E > (this_E - U)){
				m = cosf(pi*get_rand(rn)/2.0);
				sampled_E = -T * ( m*m*logf(get_rand(rn))  +   logf(get_rand(rn)) );
			}

		}
		else if (law==9){   //evaporation spectrum

			// get tabulated temperature  
			float t0 = this_Earray[ offset              ];
			float t1 = this_Earray[ offset + 1          ];
			float U  = this_Earray[ offset + vlen       ];
			      e0 = this_Earray[ offset + vlen*2     ];
			      e1 = this_Earray[ offset + vlen*2 + 1 ];
			float  T = 0.0;

			// interpolate T
			if (e1==e0){  // in top bin, both values are the same
				T = t0;
			}
			else if (intt==2){// lin-lin interpolation
				m = (this_E - e0)/(e1 - e0);
                T = (1.0 - m)*t0 + m*t1;
			}
			else if(intt==1){// histogram interpolation
				T  = (t1 - t0)/(e1 - e0) * this_E + t0;
			}

			// rejection sample
			m  = (this_E - U)/T;
			e0 = 1.0-expf(-m);
			x  = -logf(1.0-e0*get_rand(rn)) - logf(1.0-e0*get_rand(rn));
			while (  x>m ) {
				x  = -logf(1.0-e0*get_rand(rn)) - logf(1.0-e0*get_rand(rn));
			}

			// mcnp5 volIII pg 2-43
			sampled_E = T * x;

			//isotropic mu
			mu  = 2.0*get_rand(rn)-1.0;

		}
		else{
			printf("LAW %u NOT HANDLED IN POP!\n",law);
		}

		//sample isotropic directions
		rn1 = get_rand(rn);
		rn2 = get_rand(rn);
		mu  = 2.0*rn1-1.0; 
		phi = 2.0*pi*rn2;
		x = sqrtf(1.0-(mu*mu))*cosf(phi);
		y = sqrtf(1.0-(mu*mu))*sinf(phi);
		z = mu;
		//printf("% 6.4E % 6.4E % 6.4E\n",x,y,z);
	
		//check limits
		if (sampled_E >= Emax){sampled_E = Emax * 0.99;}//printf("enforcing limits in pop data_dex=%u, sampled_E = %6.4E\n",data_dex,sampled_E);}
		if (sampled_E <= Emin){sampled_E = Emin * 1.01;}//printf("enforcing limits in pop data_dex=%u, sampled_E = %6.4E\n",data_dex,sampled_E);}

		// sync before writes
		__syncthreads();

		// set data
		//printf("(xyz) %6.4E %6.4E %6.4E (dir) %6.4E %6.4E %6.4E E %6.4E\n",this_space.x,this_space.y,this_space.z,this_space.xhat,this_space.yhat,this_space.zhat,sampled_E);
		space_out[ data_dex ].x 			= this_space.x;
		space_out[ data_dex ].y 			= this_space.y;
		space_out[ data_dex ].z 			= this_space.z;
		space_out[ data_dex ].xhat 			= x;
		space_out[ data_dex ].yhat 			= y;
		space_out[ data_dex ].zhat 			= z;
		space_out[ data_dex ].enforce_BC 	= 0;
		space_out[ data_dex ].surf_dist 	= 99999.0;
		space_out[ data_dex ].macro_t 		= 8.675309;
		E_out 	 [ data_dex ] 				= sampled_E;
	}

}
//__device__ void process_multiplicity(unsigned this_yield, unsigned* rn, unsigned position, unsigned this_tope, unsigned this_awr, float this_E, source_point this_space, float* this_Earray, float* this_Sarray, source_point* space_out, float* E_out){
//
//	//constants
//	//const float  pi           =   3.14159265359 ;
//	const float  m_n  =   1.00866491600 ; // u
//	const float  Emin =   1e-11;
//	const float  Emax =   20.0; //MeV
//
//	// internal kernel variables
//	wfloat3 	hats_old(this_space.xhat,this_space.yhat,this_space.zhat);
//	float 		mu, next_E, last_E, sampled_E, e_start, E0, E1, Ek, next_e_start, next_e_end, last_e_start, last_e_end, diff;
//    unsigned 	k, vlen, next_vlen, offset, n, law, data_dex, intt; 
//	float  		speed_n          	=   sqrtf(2.0*this_E/m_n);
//	wfloat3 	v_n_cm,v_t_cm,v_n_lf,v_t_lf,v_cm, hats_new, hats_target;
//	float 		cdf0,e0,A,R,pdf0,rn1,cdf1,pdf1,e1;
//
//	//get proper data index
//	data_dex = position+k ;
//	
//	// make speed vectors
//	v_n_lf = hats_old    * speed_n;
//	v_t_lf = hats_target * 0.0;
//
//	// calculate  v_cm
//	v_cm = (v_n_lf + (v_t_lf*this_awr))/(1.0+this_awr);
//
//	//transform neutron velocity into CM frame
//	v_n_cm = v_n_lf - v_cm;
//	v_t_cm = v_t_lf - v_cm;
//
//	//
//	//sample energy
//	//
//	//read in values
//	offset = 6;
//	memcpy(&last_E,   	&this_Earray[0], sizeof(float));
//	memcpy(&next_E,   	&this_Earray[1], sizeof(float));
//	memcpy(&vlen,   	&this_Earray[2], sizeof(float));
//	memcpy(&next_vlen,	&this_Earray[3], sizeof(float));
//	memcpy(&law, 		&this_Earray[4], sizeof(float));
//	memcpy(&intt, 		&this_Earray[5], sizeof(float));
//
//	if (law==4){
//		float r = (this_E-last_E)/(next_E-last_E);
//		last_e_start = this_Earray[ offset ];
//		last_e_end   = this_Earray[ offset + vlen - 1 ];
//		next_e_start = this_Earray[ offset + 3*vlen ];
//		next_e_end   = this_Earray[ offset + 3*vlen + next_vlen - 1];
//		rn1 = get_rand(rn);
//		float rn2 = get_rand(rn);
//	
//		//sample energy dist
//		sampled_E = 0.0;
//		if(  rn2 >= r ){   //sample last E
//			diff = next_e_end - next_e_start;
//			e_start = next_e_start;
//			for ( n=0 ; n<vlen-1 ; n++ ){
//				cdf0 		= this_Earray[ (offset +   vlen ) + n+0];
//				cdf1 		= this_Earray[ (offset +   vlen ) + n+1];
//				pdf0		= this_Earray[ (offset + 2*vlen ) + n+0];
//				pdf1		= this_Earray[ (offset + 2*vlen ) + n+1];
//				e0  		= this_Earray[ (offset          ) + n+0];
//				e1  		= this_Earray[ (offset          ) + n+1]; 
//				if( rn1 >= cdf0 & rn1 < cdf1 ){
//					break;
//				}
//			}
//		}
//		else{
//			diff = next_e_end - next_e_start;
//			e_start = next_e_start;
//			for ( n=0 ; n<next_vlen-1 ; n++ ){
//				cdf0 		= this_Earray[ (offset + 3*vlen +   next_vlen ) + n+0];
//				cdf1  		= this_Earray[ (offset + 3*vlen +   next_vlen ) + n+1];
//				pdf0		= this_Earray[ (offset + 3*vlen + 2*next_vlen ) + n+0];
//				pdf1		= this_Earray[ (offset + 3*vlen + 2*next_vlen ) + n+1];
//				e0   		= this_Earray[ (offset + 3*vlen               ) + n+0];
//				e1   		= this_Earray[ (offset + 3*vlen               ) + n+1];
//				if( rn1 >= cdf0 & rn1 < cdf1 ){
//					break;
//				}
//			}
//		}
//		
//		if (intt==2){// lin-lin interpolation
//			float m 	= (pdf1 - pdf0)/(e1-e0);
//			float arg = pdf0*pdf0 + 2.0 * m * (rn1-cdf0);
//			if(arg<0){
//				E0 = e0 + (e1-e0)/(cdf1-cdf0)*(rn1-cdf0);
//			}
//			else{
//				E0 	= e0 + (  sqrtf( arg ) - pdf0) / m ;
//			}
//		}
//		else if(intt==1){// histogram interpolation
//			E0 = e0 + (rn1-cdf0)/pdf0;
//		}
//		
//		//scale it
//		E1 = last_e_start + r*( next_e_start - last_e_start );
//		Ek = last_e_end   + r*( next_e_end   - last_e_end   );
//		sampled_E = E1 +(E0-e_start)*(Ek-E1)/diff;
//
//		//isotropic mu
//		mu  = 2.0*get_rand(rn)-1.0;
//
//	}	
//
//	else if (law==9){   //evaopration spectrum
//
//		// get tabulated temperature
//		float t0 = this_Earray[ offset              ];
//		float t1 = this_Earray[ offset + 1          ];
//		float U  = this_Earray[ offset + vlen       ];
//		      e0 = this_Earray[ offset + vlen*2     ];
//		      e1 = this_Earray[ offset + vlen*2 + 1 ];
//		float  T = 0.0;
//		float  m = 0.0;
//
//		// interpolate T
//		if (e1==e0){  // in top bin, both values are the same
//				T = t0;
//			}
//		else if (intt==2){// lin-lin interpolation
//			m = (this_E - e0)/(e1 - e0);
//            T = (1.0 - m)*t0 + m*t1;
//		}
//		else if(intt==1){// histogram interpolation
//			T  = (t1 - t0)/(e1 - e0) * this_E + t0;
//		}
//
//		// rejection sample
//		m  = (this_E - U)/T;
//		e0 = 1.0-expf(-m);
//		float x  = -logf(1.0-e0*get_rand(rn)) - logf(1.0-e0*get_rand(rn));
//		while (  x>m ) {
//			x  = -logf(1.0-e0*get_rand(rn)) - logf(1.0-e0*get_rand(rn));
//		}
//
//		// mcnp5 volIII pg 2-43
//		sampled_E = T * x;
//
//		//isotropic mu
//		mu  = 2.0*get_rand(rn)-1.0;
//
//	}
//	else{
//		printf("LAW %u NOT HANDLED IN multiplicity POP!\n",law);
//	}
//		
//	// rotate direction vector
//	hats_old = v_n_cm / v_n_cm.norm2();
//	hats_old = hats_old.rotate(mu, get_rand(rn));
//	
//	//  scale to sampled energy
//	v_n_cm = hats_old * sqrtf(2.0*sampled_E/m_n);
//	
//	// transform back to L
//	v_n_lf = v_n_cm + v_cm;
//	hats_new = v_n_lf / v_n_lf.norm2();
//	hats_new = hats_new / hats_new.norm2(); // get higher precision, make SURE vector is length one
//	
//	// calculate energy in lab frame
//	//sampled_E = 0.5 * m_n * v_n_lf.dot(v_n_lf);
//
//	//check limits
//	if (sampled_E >= Emax){sampled_E = Emax * 0.99;}//printf("enforcing limits in pop data_dex=%u, sampled_E = %6.4E\n",data_dex,sampled_E);}
//	if (sampled_E <= Emin){sampled_E = Emin * 1.01;}//printf("enforcing limits in pop data_dex=%u, sampled_E = %6.4E\n",data_dex,sampled_E);}
//
//	// sync before writes
//	__syncthreads();
//
//	// write results
//	space_out[ data_dex ].x 			= this_space.x;
//	space_out[ data_dex ].y 			= this_space.y;
//	space_out[ data_dex ].z 			= this_space.z;
//	space_out[ data_dex ].xhat 			= hats_new.x;
//	space_out[ data_dex ].yhat 			= hats_new.y;
//	space_out[ data_dex ].zhat 			= hats_new.z;
//	space_out[ data_dex ].enforce_BC 	= 0;
//	space_out[ data_dex ].surf_dist 	= 99999.0;
//	space_out[ data_dex ].macro_t 		= 8.675309;
//	E_out 	 [ data_dex ] 				= sampled_E;
//
//}
__global__ void pop_source_kernel(unsigned N, unsigned* isonum, unsigned* completed, unsigned* scanned, unsigned* remap, unsigned* yield, unsigned* done, unsigned* index, unsigned* rxn, source_point* space, float* E , unsigned* rn_bank, float**  energydata, float**  scatterdata, source_point* space_out, float* E_out, float * awr_list){

	int tid = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid >= N){return;}

	// return if no yield
	unsigned 		this_yield 	= yield[tid];
	if (this_yield==0){return;}

	// load in external data
	source_point 	this_space 	= space  [tid];
	unsigned 		position 	= scanned[tid];
	unsigned 		this_tope 	= isonum [tid];
	unsigned 		dex  		= index  [tid];
	unsigned 		rn 			= rn_bank[tid];
	unsigned 		this_rxn 	= rxn    [tid];
	float 			this_E 		= E      [tid]; 

	__syncthreads();

	float*	 		this_Sarray = scatterdata[dex];
	float*			this_Earray = energydata [dex];

	__syncthreads();

	// check data array pointers
	if(this_Earray == 0x0){
		printf("null pointer in pop Earray,tid %u dex %u rxn %u tope %u E %6.4E\n",tid,dex,this_rxn,this_tope,this_E);
		return;
	}
	//if(this_Sarray == 0x0){  //Sarray value is nu for 918, might have to put check in here
	//	printf("null pointer in pop Sarray!,tid %u dex %u rxn %u tope %u E %6.4E\n",tid,dex,this_rxn,this_tope,this_E);
	//	return;
	//}

	// sampled based on reaction type
	if(this_rxn>=916 & this_rxn<=945 ){
		     process_fission(this_yield, &rn, position, this_tope, awr_list[this_tope], this_E, this_space, this_Earray, this_Sarray,  space_out, E_out);
	//}
	//else if(this_rxn == 916 | this_rxn==924 | this_rxn == 911 | this_rxn == 924 | this_rxn == 929 | this_rxn == 930 | this_rxn == 941 | this_rxn == 917 | this_rxn == 925 | this_rxn == 942 ){
	//else if (){
		//printf("processing multiplicity rxn %d\n",this_rxn);
	//	process_multiplicity(this_yield, &rn, position, this_tope, awr_list[this_tope], this_E, this_space, this_Earray, this_Sarray, space_out, E_out);
	}
	else{
		printf("tid %u REACTION %u HAS NONZERO YIELD IN SOURCE POP!\n",tid,this_rxn);
		return;
	}

	// update rn bank
	rn_bank[tid] = rn;

}

void pop_source( unsigned NUM_THREADS,  unsigned N, unsigned* isonum, unsigned* d_completed, unsigned* d_scanned, unsigned* d_remap, unsigned* d_yield, unsigned* d_done, unsigned* d_index, unsigned* d_rxn, source_point* d_space, float* d_E , unsigned* d_rn_bank, float ** energydata, float** scatterdata, source_point* space_out, float* E_out, float * awr_list){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	pop_source_kernel <<< blks, NUM_THREADS >>> ( N, isonum, d_completed, d_scanned, d_remap, d_yield, d_done, d_index, d_rxn, d_space, d_E , d_rn_bank, energydata, scatterdata, space_out, E_out, awr_list);
	cudaThreadSynchronize();

}

