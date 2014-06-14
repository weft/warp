#include <iostream>
#include <sstream>
#include <stdio.h>
#include <cuda.h>

void write_to_file(unsigned* array_in , unsigned N , std::string filename){

	FILE* f  = fopen(filename.c_str(),"w");
	unsigned * hostdata = new unsigned [N];
	cudaMemcpy(hostdata,array_in,N*sizeof(unsigned),cudaMemcpyDeviceToHost);

	for(unsigned k = 0;  k<N ;k++){
		fprintf(f,"%u\n",hostdata[k]);
	}

	delete hostdata;
	fclose(f);

}

