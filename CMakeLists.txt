#run command
#cmake -DCMAKE_C_COMPILER:FILEPATH=/usr/bin/gcc-4.4 ..

#set project name
project(WARP)

#define cmake commands/policies
#need >= v2.8.10, see http://redmine.gromacs.org/issues/1051
cmake_minimum_required(VERSION 2.8.10)
find_package(CUDA QUIET REQUIRED)
set(CMAKE_MODULE_PATH "/usr/local/OptiX-3.0.1/SDK/CMake/;${CMAKE_MODULE_PATH}")
include_directories("/usr/local/OptiX-3.0.1/include/")
find_package(OptiX REQUIRED)

set(CUDA_64_BIT_DEVICE_CODE ON)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -arch=sm_20; --use_fast_math; --compiler-options '-fPIC')

cuda_compile_ptx(camera camera.cu optix.h)

add_custom_target(ptx ALL
    DEPENDS ${camera} camera.cu
    SOURCES camera.cu)

#build CUDA executables
#cuda_add_executable(set_positions_rand set_positions_rand.cu)
