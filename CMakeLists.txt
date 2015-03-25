#set project name
project(WARP)

#define cmake commands/policies
cmake_minimum_required(VERSION 2.8)
find_package(CUDA QUIET REQUIRED)

#pass options to NVCC
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -arch sm_20 -m64  -use_fast_math --compiler-options '-fPIC' --compiler-bindir '/usr/bin/g++-4.4')

#build CUDA executables
cuda_add_executable(set_positions_rand set_positions_rand.cu)
