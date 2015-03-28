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
add_custom_command(OUTPUT hits.ptx COMMAND ${CMAKE_COMMAND} -E rename ${hits} hits.ptx DEPENDS ${hits}) 
add_custom_command(OUTPUT miss.ptx COMMAND ${CMAKE_COMMAND} -E rename ${miss} miss.ptx DEPENDS ${miss}) 
add_custom_command(OUTPUT box.ptx COMMAND ${CMAKE_COMMAND} -E rename ${box} box.ptx DEPENDS ${box}) 
add_custom_command(OUTPUT cylinder.ptx COMMAND ${CMAKE_COMMAND} -E rename ${cylinder} cylinder.ptx DEPENDS ${cylinder}) 
add_custom_command(OUTPUT hex.ptx COMMAND ${CMAKE_COMMAND} -E rename ${hex} hex.ptx DEPENDS ${hex}) 
add_custom_command(OUTPUT hits_mesh.ptx COMMAND ${CMAKE_COMMAND} -E rename ${hits_mesh} hits_mesh.ptx DEPENDS ${hits_mesh}) 
add_custom_command(OUTPUT box_mesh.ptx COMMAND ${CMAKE_COMMAND} -E rename ${box_mesh} box_mesh.ptx DEPENDS ${box_mesh}) 
add_custom_command(OUTPUT cylinder_mesh.ptx COMMAND ${CMAKE_COMMAND} -E rename ${cylinder_mesh} cylinder_mesh.ptx DEPENDS ${cylinder_mesh}) 
add_custom_command(OUTPUT hex_mesh.ptx COMMAND ${CMAKE_COMMAND} -E rename ${hex_mesh} hex_mesh.ptx DEPENDS ${hex_mesh}) 
add_custom_command(OUTPUT sphere_mesh.ptx COMMAND ${CMAKE_COMMAND} -E rename ${sphere_mesh} sphere_mesh.ptx DEPENDS ${sphere_mesh}) 
add_custom_target(camera ALL
    DEPENDS camera.ptx camera.cu
    SOURCES camera.cu)
add_custom_target(hits ALL
    DEPENDS hits.ptx hits.cu
    SOURCES hits.cu)
add_custom_target(miss ALL
    DEPENDS miss.ptx miss.cu
    SOURCES miss.cu)
add_custom_target(box ALL
    DEPENDS box.ptx box.cu
    SOURCES box.cu)
add_custom_target(cylinder ALL
    DEPENDS cylinder.ptx cylinder.cu
    SOURCES cylinder.cu)
add_custom_target(hex ALL
    DEPENDS hex.ptx hex.cu
    SOURCES hex.cu)
add_custom_target(hits_mesh ALL
    DEPENDS hits_mesh.ptx hits_mesh.cu
    SOURCES hits_mesh.cu)
add_custom_target(box_mesh ALL
    DEPENDS box_mesh.ptx box_mesh.cu
    SOURCES box_mesh.cu)
add_custom_target(cylinder_mesh ALL
    DEPENDS cylinder_mesh.ptx cylinder_mesh.cu
    SOURCES cylinder_mesh.cu)
add_custom_target(hex_mesh ALL
    DEPENDS hex_mesh.ptx hex_mesh.cu
    SOURCES hex_mesh.cu)
add_custom_target(sphere_mesh ALL
    DEPENDS sphere_mesh.ptx sphere_mesh.cu
    SOURCES sphere_mesh.cu)

#build CUDA executables
cuda_add_executable(set_positions_rand set_positions_rand.cu)
