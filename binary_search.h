#ifndef BINARY_SEARCH_H
#define BINARY_SEARCH_H

__forceinline__ __device__ unsigned binary_search( float * array , float value, unsigned len ){

	// load data
	unsigned donesearching = 0;
	unsigned cnt  = 1;
	unsigned powtwo = 2;
	int dex  = (len) / 2;  //N_energies starts at 1, duh

	// edge check
	if(value < array[0] | value > array[len-1]){
		//printf("device binary search value outside array range! %p %d val % 10.8f ends % 10.8f % 10.8f\n",array,len,value,array[0],array[len-1]);
		//printf("val %6.4E len %u outside %6.4E %6.4E %6.4E %6.4E %6.4E %6.4E ... %6.4E %6.4E\n",value,len,array[0],array[1],array[2],array[3],array[4],array[5],array[len-1],array[len]);
		//return 0;
	}

	// search
	while(!donesearching){

		powtwo = powtwo * 2;
		if      ( 	array[dex]   <= value && 
					array[dex+1] >  value ) { donesearching = 1; }
		else if ( 	array[dex]   >  value ) { dex  = dex - (( len / powtwo) + 1) ; cnt++; }  // +1's are to do a ceiling instead of a floor on integer division
		else if ( 	array[dex]   <  value ) { dex  = dex + (( len / powtwo) + 1) ; cnt++; }

		if(cnt>30){
			donesearching=1;
			printf("device binary search iteration overflow! dex %d ptr %p %d val % 10.8f ends % 10.8f % 10.8f\n",dex,array,len,value,array[0],array[len-1]);
			dex=0;
		}

		// edge checks... fix later???
		if(dex<0){
			dex=0;
		}
		if(dex>=len){
			dex=len-1;
		}
	}

	// output index
	return dex;

}

#endif
