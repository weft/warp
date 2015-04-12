#run command
#cmake -DCMAKE_C_COMPILER:FILEPATH=/usr/bin/gcc-4.4 ..

#set project name
project(WARP)

#define cmake commands/policies
#need >= v2.8.10, see http://redmine.gromacs.org/issues/1051
cmake_minimum_required(VERSION 2.8.10)

set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} LD_LIBRARY_PATH)

#find CUDA; we need that
find_package(CUDA QUIET REQUIRED)

#find OptiX; also required
#see http://bit.ly/19nExJN for inspiration
#note: had to append "/usr/local/OptiX-3.0.1/lib64" to $PATH in order for optix.h & friends to be found
#i'm having to export that at the front of $PATH each time? strange...
#set(OptiX_INSTALL_DIR "/usr/local/OptiX-3.0.1")
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${OptiX_INSTALL_DIR}/SDK/CMake")
#set(OptiX_INSTALL_DIR "${CMAKE_SOURCE_DIR}/../" CACHE PATH "Path to OptiX installed location.")
#set(OptiX_DIR "${OptiX_INSTALL_DIR}/SDK/CMake")
#find_package(OptiX REQUIRED)
#include_directories("${OptiX_INCLUDE}")

set(OptiX_ROOT_DIR "/usr/local/OptiX-3.0.1")

FIND_PATH(OptiX_INCLUDE_DIR optix.h
   ${OptiX_ROOT_DIR}/include
)

FIND_LIBRARY(OptiX_LIBRARY
  NAMES optix optix${CMAKE_SIZEOF_VOID_P}
  PATHS
     ${OptiX_ROOT_DIR}/lib64
)

IF(OptiX_LIBRARY)
  IF (OptiX_INCLUDE_DIR)
    # OK, found all we need
    SET(OptiX_FOUND TRUE)		
    GET_FILENAME_COMPONENT(OptiX_LINK_DIRECTORIES ${OptiX_LIBRARY} PATH)

  ELSE (OptiX_INCLUDE_DIR)
    MESSAGE("OptiX include dir not found. Set OptiX_ROOT_DIR to find it.")
  ENDIF(OptiX_INCLUDE_DIR)
ELSE(OptiX_LIBRARY)
  MESSAGE("OptiX lib not found. Set OptiX_ROOT_DIR to find it.")
ENDIF(OptiX_LIBRARY)

MARK_AS_ADVANCED(
  OptiX_INCLUDE_DIR
  OptiX_LIBRARY
  OptiX_LINK_DIRECTORIES		
)

include_directories(${OptiX_INCLUDE_DIR})

set(CUDA_64_BIT_DEVICE_CODE ON)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -arch=sm_20; --use_fast_math; --compiler-options '-fPIC')

include_directories("${CUDA_TOOLKIT_INCLUDE}")
include_directories("${CUDA_INCLUDE_DIRS}")
set (CUDPP_ROOT_DIR "/home/krowland/cudpp-2.1")
#http://suren.me/webbzr/normxcorr/trunk/annotate/27/dict_hw/cmake/FindCUDPP.cmake?remember=7
FIND_PATH(CUDPP_INCLUDE_DIR cudpp.h
   ${CUDPP_ROOT_DIR}/cudpp/include
   ${CUDPP_ROOT_DIR}/include
)

FIND_LIBRARY(CUDPP_LIBRARY
  NAMES cudpp cudpp${CMAKE_SIZEOF_VOID_P}
  PATHS
     ${CUDPP_ROOT_DIR}/lib
)

IF(CUDPP_LIBRARY)
  IF (CUDPP_INCLUDE_DIR)
    # OK, found all we need
    SET(CUDPP_FOUND TRUE)		
    GET_FILENAME_COMPONENT(CUDPP_LINK_DIRECTORIES ${CUDPP_LIBRARY} PATH)

  ELSE (CUDPP_INCLUDE_DIR)
    MESSAGE("CUDPP include dir not found. Set CUDPP_ROOT_DIR to find it.")
  ENDIF(CUDPP_INCLUDE_DIR)
ELSE(CUDPP_LIBRARY)
  MESSAGE("CUDPP lib not found. Set CUDPP_ROOT_DIR to find it.")
ENDIF(CUDPP_LIBRARY)

MARK_AS_ADVANCED(
  CUDPP_INCLUDE_DIR	
  CUDPP_LIBRARY
  CUDPP_LINK_DIRECTORIES		
)

SET(libs ${libs} ${CUDPP_LIBRARY})
SET(libs ${libs} ${OptiX_LIBRARY})
SET(libs ${libs} ${CUDA_LIBRARY})
SET(libs ${libs} ${CUDA_CUDART_LIBRARY})
SET(libs ${libs} ${CUDA_CUBLAS_LIBRARY})

include_directories(${CUDPP_INCLUDE_DIR})

#http://stackoverflow.com/questions/11041299/python-h-no-such-file-or-directory
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

find_package(PNG REQUIRED)

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
cuda_compile(set_positions_rand set_positions_rand.cu)
cuda_compile(find_E_grid_index find_E_grid_index.cu)
cuda_compile(find_E_grid_index_quad find_E_grid_index_quad.cu)
cuda_compile(sample_fission_spectra sample_fission_spectra.cu)
cuda_compile(sample_fixed_source sample_fixed_source.cu)
cuda_compile(macrosocpic macroscopic.cu)
cuda_compile(microscopic microscopic.cu)
cuda_compile(copy_points copy_points.cu)
cuda_compile(tally_spec tally_spec.cu)
cuda_compile(escatter escatter.cu)
cuda_compile(iscatter iscatter.cu)
cuda_compile(cscatter cscatter.cu)
cuda_compile(fission fission.cu)
cuda_compile(absorb absorb.cu)
cuda_compile(print_histories print_histories.cu)
cuda_compile(pop_secondaries pop_secondaries.cu)
cuda_compile(pop_source pop_source.cu)
cuda_compile(rebase_yield rebase_yield.cu)
cuda_compile(flip_done flip_done.cu)
cuda_compile(device_copies device_copies.cu)
cuda_compile(reaction_edges3 reaction_edges3.cu) 
cuda_compile(write_to_file write_to_file.cu)

# Sources:
SET(warp_executable_srcs
    print_banner.cpp
    warp.pyc
    primitive.cpp
    unionize.pyc
    main.cpp
#    primitive_test.cpp
    mt19937ar.cpp
    optix_stuff.cpp
    whistory.cpp
#    optixtest.cpp
#    optix_stuff_test.cpp
#    reaction_edges_test.cpp
#    wgeometry_test.cpp
#    whistory_test.cpp
    wgeometry.cpp
)

SET(warp_cuda_srcs
    ${set_positions_rand}
    ${find_E_grid_index}
    ${find_E_grid_index_quad}
    ${sample_fission_spectra}
    ${sample_fixed_source}
    ${macroscopic}
    ${microscopic}
    ${copy_points}
    ${tally_spec}
    ${escatter}
    ${iscatter}
    ${cscatter}
    ${fission}
    ${absorb}
    ${print_histories}
    ${pop_secondaries}
    ${pop_source}
    ${rebase_yield}
    ${flip_done}
    ${device_copies}
    ${reaction_edges3}
    ${write_to_file}
)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

# actual target:
ADD_EXECUTABLE(warp ${warp_executable_srcs} ${warp_cuda_srcs})

# link libs
TARGET_LINK_LIBRARIES(warp ${libs})

# add install target:
INSTALL(TARGETS warp DESTINATION bin)

#note: how to do "block comments"
#if(FALSE)
#lines to be commented out
#endif()

#look into:
#add target library
#link 
