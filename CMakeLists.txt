#set project name
project(WARP)

#define cmake commands/policies
#need >= v2.8.10, see http://redmine.gromacs.org/issues/1051
cmake_minimum_required(VERSION 2.8.10)

set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} LD_LIBRARY_PATH)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake") 
include_directories(${CMAKE_SOURCE_DIR})

#find CUDA; we need that
find_package(CUDA QUIET REQUIRED)
include_directories("${CUDA_TOOLKIT_INCLUDE}")
include_directories("${CUDA_INCLUDE_DIRS}")

set(OptiX_ROOT_DIR CACHE PATH "Root directory of OptiX")

find_path(OptiX_INCLUDE_DIR optix.h
   ${OptiX_ROOT_DIR}/include
)

find_library(OptiX_LIBRARY
  NAMES optix optix${CMAKE_SIZEOF_VOID_P}
  PATHS
     ${OptiX_ROOT_DIR}/lib64
)

if(OptiX_LIBRARY)
  if (OptiX_INCLUDE_DIR)
    # OK, found all we need
    set(OptiX_FOUND TRUE)		
    get_filename_component(OptiX_LINK_DIRECTORIES ${OptiX_LIBRARY} PATH)

  else (OptiX_INCLUDE_DIR)
    message("OptiX include dir not found. Set OptiX_ROOT_DIR to find it.")
  endif(OptiX_INCLUDE_DIR)
else(OptiX_LIBRARY)
  message("OptiX lib not found. Set OptiX_ROOT_DIR to find it.")
endif(OptiX_LIBRARY)

mark_as_advanced(
  OptiX_INCLUDE_DIR
  OptiX_LIBRARY
  OptiX_LINK_DIRECTORIES		
)

include_directories(${OptiX_INCLUDE_DIR})

set(CUDA_64_BIT_DEVICE_CODE ON)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -arch=sm_30; --use_fast_math; --compiler-options '-fPIC')

set (CUDPP_ROOT_DIR CACHE PATH "Root directory for CUDPP")
#http://suren.me/webbzr/normxcorr/trunk/annotate/27/dict_hw/cmake/FindCUDPP.cmake?remember=7
find_path(CUDPP_INCLUDE_DIR cudpp.h
   ${CUDPP_ROOT_DIR}/cudpp/include
   ${CUDPP_ROOT_DIR}/include
)

find_library(CUDPP_LIBRARY
  NAMES cudpp cudpp${CMAKE_SIZEOF_VOID_P}
  PATHS
     ${CUDPP_ROOT_DIR}/lib
)

if(CUDPP_LIBRARY)
  if (CUDPP_INCLUDE_DIR)
    # OK, found all we need
    set(CUDPP_FOUND TRUE)
    get_filename_component(CUDPP_LINK_DIRECTORIES ${CUDPP_LIBRARY} PATH)

  else (CUDPP_INCLUDE_DIR)
    message("CUDPP include dir not found. Set CUDPP_ROOT_DIR to find it.")
  endif(CUDPP_INCLUDE_DIR)
else(CUDPP_LIBRARY)
  message("CUDPP lib not found. Set CUDPP_ROOT_DIR to find it.")
endif(CUDPP_LIBRARY)

mark_as_advanced(
  CUDPP_INCLUDE_DIR     
  CUDPP_LIBRARY
  CUDPP_LINK_DIRECTORIES                
)

# find numpy and include the numpy headers
find_package(Numpy REQUIRED)
include_directories("${NUMPY_INCLUDE_DIR}")

set(libs ${libs} ${CUDPP_LIBRARY})
set(libs ${libs} ${OptiX_LIBRARY})
set(libs ${libs} ${CUDA_LIBRARIES})
set(libs ${libs} ${CUDA_CUDART_LIBRARY})
set(libs ${libs} ${CUDA_CUBLAS_LIBRARY})
set(libs ${libs} ${CUDA_curand_LIBRARY})

include_directories(${CUDPP_INCLUDE_DIR})

#http://stackoverflow.com/questions/11041299/python-h-no-such-file-or-directory
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
set(libs ${libs} ${PYTHON_LIBRARIES})

find_package(PNG REQUIRED)
set(libs ${libs} ${PNG_LIBRARY})

#mac stuff
if(APPLE)
 set(CMAKE_OSX_VERSION "10.11" CACHE STRING "OSX Version number")
 include_directories("/usr/local/include")
 set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mmacosx-version-min=${CMAKE_OSX_VERSION}")
 set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mmacosx-version-min=${CMAKE_OSX_VERSION}")
 set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} --compiler-options '-mmacosx-version-min=${CMAKE_OSX_VERSION}'")
endif(APPLE)

include( UseCython )

execute_process(COMMAND "${CYTHON_EXECUTABLE}" "-V"
    ERROR_VARIABLE CYTHON_VERSION
    ERROR_STRIP_TRAILING_WHITESPACE)
message( STATUS "Cython Version: " ${CYTHON_VERSION} )



cuda_compile_ptx(camera camera.cu)
cuda_compile_ptx(miss miss.cu)
cuda_compile_ptx(hits_mesh hits_mesh.cu)
cuda_compile_ptx(box_mesh box_mesh.cu)
cuda_compile_ptx(cylinder_mesh cylinder_mesh.cu)
cuda_compile_ptx(hex_mesh hex_mesh.cu)
cuda_compile_ptx(sphere_mesh sphere_mesh.cu)
#http://stackoverflow.com/questions/26198294/nvcc-compile-to-ptx-using-cmakes-cuda-compile-ptx
add_custom_command(OUTPUT camera.ptx COMMAND ${CMAKE_COMMAND} -E rename ${camera} camera.ptx DEPENDS ${camera}) 
add_custom_command(OUTPUT miss.ptx COMMAND ${CMAKE_COMMAND} -E rename ${miss} miss.ptx DEPENDS ${miss}) 
add_custom_command(OUTPUT hits_mesh.ptx COMMAND ${CMAKE_COMMAND} -E rename ${hits_mesh} hits_mesh.ptx DEPENDS ${hits_mesh}) 
add_custom_command(OUTPUT box_mesh.ptx COMMAND ${CMAKE_COMMAND} -E rename ${box_mesh} box_mesh.ptx DEPENDS ${box_mesh}) 
add_custom_command(OUTPUT cylinder_mesh.ptx COMMAND ${CMAKE_COMMAND} -E rename ${cylinder_mesh} cylinder_mesh.ptx DEPENDS ${cylinder_mesh}) 
add_custom_command(OUTPUT hex_mesh.ptx COMMAND ${CMAKE_COMMAND} -E rename ${hex_mesh} hex_mesh.ptx DEPENDS ${hex_mesh}) 
add_custom_command(OUTPUT sphere_mesh.ptx COMMAND ${CMAKE_COMMAND} -E rename ${sphere_mesh} sphere_mesh.ptx DEPENDS ${sphere_mesh}) 
add_custom_target(camera ALL
    DEPENDS camera.ptx camera.cu
    SOURCES camera.cu)
add_custom_target(miss ALL
    DEPENDS miss.ptx miss.cu
    SOURCES miss.cu)
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
cuda_compile(sample_fissile_energy sample_fissile_energy.cu)
cuda_compile(macro_micro macro_micro.cu)
cuda_compile(copy_points copy_points.cu)
cuda_compile(scatter_level scatter_level.cu)
cuda_compile(scatter_conti scatter_conti.cu)
cuda_compile(scatter_multi scatter_multi.cu)
cuda_compile(fission fission.cu)
cuda_compile(pop_fission pop_fission.cu)
cuda_compile(rebase_yield rebase_yield.cu)
cuda_compile(reaction_edges3 reaction_edges3.cu) 
cuda_compile(write_to_file write_to_file.cu)
cuda_compile(safety_check safety_check.cu)
cuda_compile(check_pointers check_pointers.cu)

# Sources:
set(warp_executable_srcs
    wprimitive.cpp
    mt19937ar.cpp
    optix_stuff.cpp
    whistory.cpp
    wgeometry.cpp
)

set(warp_cuda_srcs
    ${set_positions_rand}
    ${find_E_grid_index}
    ${sample_fissile_energy}
    ${macro_micro}
    ${copy_points}
    ${scatter_level}
    ${scatter_conti}
    ${scatter_multi}
    ${fission}
    ${pop_fission}
    ${rebase_yield}
    ${reaction_edges3}
    ${write_to_file}
    ${safety_check}
    ${check_pointers}
)

ADD_LIBRARY(warp SHARED ${warp_executable_srcs} ${warp_cuda_srcs})

# link libs
TARGET_LINK_LIBRARIES(warp ${libs})

# actual target:
add_executable(dissertation ${warp_executable_srcs} ${warp_cuda_srcs} main.cpp)
add_executable(benchmarks   ${warp_executable_srcs} ${warp_cuda_srcs} benchmarks.cpp)

# standalone ace
set_source_files_properties("${PROJECT_SOURCE_DIR}/ace.pyx"
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(ace ace.pyx)
#target_link_libraries(ace pyne)

# link libs
target_link_libraries(dissertation ${libs})
target_link_libraries(benchmarks ${libs})

#enable_testing()
#find_package(GTest REQUIRED)
#add_subdirectory(test)

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})
include_directories(${PYTHON_INCLUDE_PATH})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set_source_files_properties(warp.i PROPERTIES CPLUSPLUS ON)
set_source_files_properties(warp.i PROPERTIES SWIG_FLAGS "-c++")
swig_add_module(warp python warp.i)
swig_link_libraries(warp warp ${PYTHON_LIBRARIES})
