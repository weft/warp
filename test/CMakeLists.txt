find_package(Threads)
if(CMAKE_USE_PTHREADS_INIT)  # The pthreads library is available.
    set(CXX_FLAGS "${CXX_FLAGS} -DGTEST_HAS_PTHREAD=1")
endif()

include_directories(${GTEST_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR})

set(warp_test_srcs
#	optix_stuff_test.cpp
	primitive_test.cpp
#	reaction_edges_test.cpp
#	wgeometry_test.cpp
#	whistory_test.cpp
)

add_executable(warp_tests ${warp_test_srcs})
target_link_libraries(warp_tests libwarp ${GTEST_BOTH_LIBRARIES})

add_test(warp_test warp_tests)
