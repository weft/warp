inline __device__ float get_rand(unsigned* in)
{
/*
increments the random number with LCRNG 
adapated from OpenMC again
values from http://www.ams.org/journals/mcom/1999-68-225/S0025-5718-99-00996-5/S0025-5718-99-00996-5.pdf
since 32-bit math is being used, 30 bits are used here
*/
	const unsigned a   		= 116646453;		 		// multiplier
	const unsigned c   		= 7;						// constant add, must be odd
	const unsigned mask   	= 1073741823; 				// 2^30-1
	const float norm   		= 9.31322574615478515625E-10;	// 2^-30
	unsigned nextint = (a * in[0] +  c) & mask; 			// mod by truncation
	float randout = nextint*norm;
	if(randout>=1.0){
		randout=0.9999999;
		//printf("RN=1.0  %u %u %10.8E\n",in[0],nextint,randout);
	}
	in[0]=nextint;
	return randout;   						// return normalized float
}

inline __device__ float sum_cross_section( unsigned length , float energy0, float energy1, float this_E, float* multiplier, float* array0, float* array1){
/*
Calculates the sum of a cross section range.  This routine HAS a multiplier array.  Returns sum.
*/
	float macro_t_total = 0.0;

	for( int k=0; k<length; k++ ){
		//linearly interpolate and accumulate
		macro_t_total += ( (array1[k]-array0[k])/(energy1-energy0)*(this_E-energy0) + array0[k] ) * multiplier[k];
	}

	return macro_t_total;

}

inline __device__ float sum_cross_section( unsigned length , float energy0, float energy1, float this_E, float* array0, float* array1){
/*
Calculates the sum of a cross section range.  Returns sum.
*/
	float macro_t_total = 0.0;

	for( int k=0; k<length; k++ ){
		//linearly interpolate and accumulate
		macro_t_total += ( (array1[k]-array0[k])/(energy1-energy0)*(this_E-energy0) + array0[k] );
	}

	return macro_t_total;

}

inline __device__ unsigned sample_cross_section( unsigned length , float normalize, float rn, float energy0, float energy1, float this_E, float* multiplier, float* array0, float* array1){
/*
Samples the isotope/reaction once a normalization factor is known (material/isotope total macroscopic cross section).  This routine HAS a multiplier array.  Returns array index.
*/
	unsigned	index				= 0;
	float		cumulative_value	= 0.0;

	for( index=0; index<length; index++ ){
		//linearly interpolate and accumulate
		cumulative_value += ( (array1[index]-array0[index])/(energy1-energy0)*(this_E-energy0) + array0[index] ) * multiplier[index] / normalize;
		if ( rn <= cumulative_value ){
			break;
		}
	}

	return index;

}

inline __device__ unsigned sample_cross_section( unsigned length , float normalize, float rn, float energy0, float energy1, float this_E, float* array0, float* array1){
/*
Samples the isotope/reaction once a normalization factor is known (material/isotope total macroscopic cross section).  Returns array index.
*/
	unsigned	index				= 0;
	float		cumulative_value	= 0.0;

	for( index=0; index<length; index++ ){
		//linearly interpolate and accumulate
		cumulative_value += ( (array1[index]-array0[index])/(energy1-energy0)*(this_E-energy0) + array0[index] ) / normalize;
		if ( rn <= cumulative_value ){
			break;
		}
	}

	return index;

}

inline __device__ unsigned sample_probability_distribution( unsigned length , unsigned intt , float rn , float* var , float* pdf, float* cdf ){
/*
Samples a probability distribution with historgram or lin-lin interpolation.  Returns sampled index.
*/
	unsigned	index				= 0;
	float		cumulative_value	= 0.0;

	// scan the CDF,
	for(index=0;index<length;index++){

		cumulative_value += index;

	}

	if(intt==1){        // histogram interpolation

	}
	else if(intt==2){  // lin-lin interpolation

	}
	else{
		printf("INTT=%u NOT HANDLED!\n",intt);
		index = 4294967295;   // return -1
	}

	return index;

}