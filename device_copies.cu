#include <cuda.h>
#include <stdio.h>
#include "datadef.h"

/**
 * \brief function to do a host-to-device copy
 * \details something
 *
 * @param[in] dest 	  - cuda pointer on device
 * @param[in] source  - pointer on host
 * @param[in] bytes   - number of bytes to copy
 */
void copy_to_device(void* dest,void* source ,unsigned bytes){

	cudaMemcpy(dest,source,bytes,cudaMemcpyHostToDevice);

}

/**
 * \brief function to do a device-to-host copy
 *
 * @param[in] dest 	  - pointer on host
 * @param[in] source  - cuda pointer on device
 * @param[in] bytes   - number of bytes to copy
 */
void copy_from_device(void* dest,void* source ,unsigned bytes){

	cudaMemcpy(dest,source,bytes,cudaMemcpyDeviceToHost);

}

/**
 * \brief function to do a device memory allocation
 *
 * @param[in] dest  - pointer on host
 * @param[in] bytes - number of bytes to copy
 */
void allocate_on_device(void** dest,unsigned bytes){

	cudaMalloc(dest,bytes);

}

/**
 * \brief function to do a device memory allocation
 *
 * @param[in] dest  - pointer on device to deallocate
 */
void deallocate_on_device(void* dest){

	cudaFree(dest);

}