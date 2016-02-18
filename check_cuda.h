#ifndef CHECK_CUDA_H
#define CHECK_CUDA_H

// CUDA error check wrapper
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