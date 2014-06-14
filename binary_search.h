__forceinline__ __device__ unsigned binary_search( float * array , float value, unsigned len ){

	// load data
	unsigned donesearching = 0;
	unsigned cnt  = 1;
	unsigned powtwo = 2;
	int dex  = (len) / 2;  //N_energiesgth starts at 1, duh

	while(!donesearching){

		powtwo = powtwo * 2;
		if      ( 	array[dex]   <= value && 
					array[dex+1] >  value ) { donesearching = 1; }
		else if ( 	array[dex]   >  value ) { dex  = dex - (( len / powtwo) + 1) ; cnt++; }  // +1's are to do a ceiling instead of a floor on integer division
		else if ( 	array[dex]   <  value ) { dex  = dex + (( len / powtwo) + 1) ; cnt++; }

		if(cnt>30){
			donesearching=1;
			printf("device binary search iteration overflow! %p %d val % 10.8f ends % 10.8f % 10.8f\n",array,len,value,array[0],array[len-1]);
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