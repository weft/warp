#ifndef CHECK_CUDA_H
#define CHECK_CUDA_H

/**
 * \brief CUDA error check wrapper, host only
 * \details this inline function prints detailed information about the return code of host-side
 * CUDA functions and where they occur in the code
 * @param[in] code 	- CUDA error code output from some host-side CUDA function call
 * @param[in] file 	- file where erroring fuction is, written inline by preprocessor
 * @param[in] line 	- line of file where erroring fuction is, written inline by preprocessor
 * @param[in] abort - flag to exit on error, default is true
 */
__host__ inline void check_cuda(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}
#define check_cuda(ans) { check_cuda((ans), __FILE__, __LINE__); }

#endif