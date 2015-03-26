#run command
#cmake -DCMAKE_C_COMPILER:FILEPATH=/usr/bin/gcc-4.4 ..

#set project name
project(WARP)

#define cmake commands/policies
#need >= v2.8.10, see http://redmine.gromacs.org/issues/1051
cmake_minimum_required(VERSION 2.8.10)

#find CUDA; we need that
find_package(CUDA QUIET REQUIRED)

#find OptiX; also required
#see http://bit.ly/19nExJN for inspiration
#note: had to append "/usr/local/OptiX-3.0.1/lib64" to $PATH in order for optix.h & friends to be found
set(OptiX_INSTALL_DIR "/usr/local/OptiX-3.0.1")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${OptiX_INSTALL_DIR}/SDK/CMake")
find_package(OptiX REQUIRED)
include_directories("${OptiX_INCLUDE}")

set(CUDA_64_BIT_DEVICE_CODE ON)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -arch=sm_20; --use_fast_math; --compiler-options '-fPIC')

cuda_compile_ptx(camera camera.cu)
cuda_compile_ptx(hits hits.cu)
cuda_compile_ptx(miss miss.cu)
cuda_compile_ptx(box box.cu)
cuda_compile_ptx(cylinder cylinder.cu)
cuda_compile_ptx(hex hex.cu)
cuda_compile_ptx(hits_mesh hits_mesh.cu)
cuda_compile_ptx(box_mesh box_mesh.cu)
cuda_compile_ptx(cylinder_mesh cylinder_mesh.cu)
cuda_compile_ptx(hex_mesh hex_mesh.cu)
cuda_compile_ptx(sphere_mesh sphere_mesh.cu)
#http://stackoverflow.com/questions/26198294/nvcc-compile-to-ptx-using-cmakes-cuda-compile-ptx
add_custom_command(OUTPUT camera.ptx COMMAND ${CMAKE_COMMAND} -E rename ${camera} camera.ptx DEPENDS ${camera}) 
add_custom_target(ptx ALL
    DEPENDS camera.ptx camera.cu
    SOURCES camera.cu)

#build CUDA executables
#cuda_add_executable(set_positions_rand set_positions_rand.cu)
