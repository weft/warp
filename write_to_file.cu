#include <iostream>
#include <sstream>
#include <stdio.h>
#include <cuda.h>

/**
 * \brief writes CUDA array to a text file
 * \details copies the cuda array to a local buffer, writes to buffer to a new file, then frees the local memory
 * @param[in] array_in   - device pointer to array to write
 * @param[in] N          - number of elements to write
 * @param[in] filename   - name for the file
 */ 
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

