#include <cuda.h>
#include <stdio.h>
#include "datadef.h"

void copy_to_device(void* dest,void* source ,unsigned bytes){

	cudaMemcpy(dest,source,bytes,cudaMemcpyHostToDevice);

}

void copy_from_device(void* dest,void* source ,unsigned bytes){

	cudaMemcpy(dest,source,bytes,cudaMemcpyDeviceToHost);

}