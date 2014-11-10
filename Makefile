#
CC =  /usr/bin/gcc-4.4
CXX = /usr/bin/g++-4.4
OPTIX = /usr/local/OptiX-3.0.1
NVCC = nvcc
ARCH = -arch sm_20
C_FLAGS = -O3 -m64 -fPIC
NVCC_FLAGS = -m64  -use_fast_math --compiler-options '-fPIC' --compiler-bindir '/usr/bin/gcc-4.4'
CURAND_LIBS = -lcurand
OPTIX_FLAGS = -I$(OPTIX)/include -L$(OPTIX)/lib64 
OPTIX_LIBS = -loptix 
CUDA_FLAGS = -I/usr/local/cuda/include -L/usr/local/cuda/lib64
CUDPP_PATH = /home/krowland/cudpp-2.1
CUDPP_FLAGS = -I$(CUDPP_PATH)/include -L$(CUDPP_PATH)/lib
CUDPP_LIBS = -lcudpp_hash -lcudpp
PYTHON_FLAGS = -I/home/krowland/anaconda/include/python2.7 -L/home/krowland/anaconda/lib/python2.7
PYTHON_LIBS = -lpython2.7
PNG_FLAGS = -L/home/krowland/anaconda/lib
PNG_LIBS = -lpng15 
LIBS =
#google test stuff
GTEST_DIR = /home/krowland/gtest-1.7.0
USER_DIR = /home/krowland/warp/warp
CPPFLAGS += -I$(GTEST_DIR)/include
CXXFLAGS += -g -Wall -Wextra
TESTS = optix_stuff_test primitive_test wgeometry_test whistory_test
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
	        $(GTEST_DIR)/include/gtest/internal/*.h
#/google test stuff

COBJS =	mt19937ar.o \
		print_banner.o \
		set_positions_rand.o \
		copy_points.o \
		macroscopic.o \
		microscopic.o \
		find_E_grid_index.o \
		find_E_grid_index_quad.o \
		sample_fission_spectra.o \
		sample_fixed_source.o \
		tally_spec.o \
		escatter.o \
		iscatter.o \
		cscatter.o \
		fission.o \
		absorb.o \
		make_mask.o \
		print_histories.o \
		pop_secondaries.o \
		pop_source.o \
		rebase_yield.o \
		flip_done.o \
		whistory.o \
		wgeometry.o \
		optix_stuff.o \
		primitive.o \
		write_to_file.o \
		reaction_edges2.o \
		device_copies.o \

ptx_objects = 	camera.ptx \
				hits.ptx \
				miss.ptx \
				box.ptx \
				cylinder.ptx \
				hex.ptx\
				hits_mesh.ptx \
				box_mesh.ptx \
				cylinder_mesh.ptx \
				hex_mesh.ptx\
				sphere_mesh.ptx\


all:  	$(ptx_objects) \
		$(COBJS) \
		libwarp.so \
		python \ $(TESTS)

clean:
	rm -f *.ptx *.o *.so gpu debug optixtest warp_wrap.cxx warp.py

tests: optix_stuff_test primitive_test whistory_test wgeometry_test

GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

camera.ptx:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) $(OPTIX_FLAGS) $(OPTIX_LIBS) -ptx camera.cu

hits.ptx:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) $(OPTIX_FLAGS) $(OPTIX_LIBS) -ptx hits.cu

miss.ptx:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) $(OPTIX_FLAGS) $(OPTIX_LIBS) -ptx miss.cu

box.ptx:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) $(OPTIX_FLAGS) $(OPTIX_LIBS) -ptx box.cu

cylinder.ptx:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) $(OPTIX_FLAGS) $(OPTIX_LIBS) -ptx cylinder.cu

hex.ptx:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) $(OPTIX_FLAGS) $(OPTIX_LIBS) -ptx hex.cu

hits_mesh.ptx:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) $(OPTIX_FLAGS) $(OPTIX_LIBS) -ptx hits_mesh.cu

box_mesh.ptx:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) $(OPTIX_FLAGS) $(OPTIX_LIBS) -ptx box_mesh.cu

cylinder_mesh.ptx:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) $(OPTIX_FLAGS) $(OPTIX_LIBS) -ptx cylinder_mesh.cu

hex_mesh.ptx:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) $(OPTIX_FLAGS) $(OPTIX_LIBS) -ptx hex_mesh.cu

sphere_mesh.ptx:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) $(OPTIX_FLAGS) $(OPTIX_LIBS) -ptx sphere_mesh.cu

mt19937ar.o:
	$(CXX) $(C_FLAGS) -c -O mt19937ar.cpp

print_banner.o:
	$(NVCC) $(NVCC_FLAGS) -c -O print_banner.cpp

set_positions_rand.o:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) -c set_positions_rand.cu

find_E_grid_index.o:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) -c find_E_grid_index.cu

find_E_grid_index_quad.o:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) -c find_E_grid_index_quad.cu

sample_fission_spectra.o:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) -c sample_fission_spectra.cu

sample_fixed_source.o:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) -c sample_fixed_source.cu

macroscopic.o:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) -c macroscopic.cu

microscopic.o:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) -c microscopic.cu

copy_points.o:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) -c copy_points.cu

tally_spec.o:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) -c tally_spec.cu

escatter.o: 
	$(NVCC) $(ARCH) $(NVCC_FLAGS) -c escatter.cu

iscatter.o: 
	$(NVCC) $(ARCH) $(NVCC_FLAGS) -c iscatter.cu

cscatter.o: 
	$(NVCC) $(ARCH) $(NVCC_FLAGS) -c cscatter.cu

fission.o:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) -c fission.cu

absorb.o:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) -c absorb.cu

make_mask.o:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) -c make_mask.cu

print_histories.o:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) -c print_histories.cu

pop_secondaries.o:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) -c pop_secondaries.cu

pop_source.o:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) -c pop_source.cu

rebase_yield.o:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) -c rebase_yield.cu

flip_done.o:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) -c flip_done.cu

device_copies.o:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) -c device_copies.cu

reaction_edges2.o:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) -c reaction_edges2.cu

write_to_file.o:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) -c write_to_file.cu

whistory.o: $(USER_DIR)/whistory.cpp $(USER_DIR)/whistory.h $(GTEST_HEADERS)
	$(CXX) $(C_FLAGS) $(OPTIX_FLAGS) $(CUDPP_FLAGS) $(PNG_FLAGS) $(PYTHON_FLAGS) $(CUDA_FLAGS) -c whistory.cpp

wgeometry.o: $(USER_DIR)/wgeometry.cpp $(USER_DIR)/wgeometry.h $(GTEST_HEADERS)
	$(CXX) $(C_FLAGS) $(CPPFLAGS) $(CXXFLAGS) -c wgeometry.cpp

primitive.o: $(USER_DIR)/primitive.cpp $(USER_DIR)/primitive.h $(GTEST_HEADERS)
	$(CXX) $(C_FLAGS) $(CPPFLAGS) $(CXXFLAGS) -c primitive.cpp

optix_stuff.o: device_copies.o $(USER_DIR)/optix_stuff.cpp $(USER_DIR)/optix_stuff.h $(GTEST_HEADERS)
	$(CXX) $(C_FLAGS) $(CPPFLAGS) $(CXXFLAGS) $(OPTIX_FLAGS) $(OPTIX_LIBS) $(CUDA_FLAGS) $(PNG_FLAGS) -lm -c optix_stuff.cpp

#google test
optix_stuff_test.o : $(USER_DIR)/optix_stuff_test.cpp \
	$(USER_DIR)/optix_stuff.h $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPTIX_FLAGS) $(CUDPP_FLAGS) $(PYTHON_FLAGS) $(CUDA_FLAGS) $(PNG_FLAGS) -c $(USER_DIR)/optix_stuff_test.cpp

primitive_test.o : $(USER_DIR)/primitive_test.cpp \
	$(USER_DIR)/primitive.h $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/primitive_test.cpp

wgeometry_test.o : $(USER_DIR)/wgeometry_test.cpp \
	$(USER_DIR)/wgeometry.h $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/wgeometry_test.cpp

whistory_test.o : $(USER_DIR)/whistory_test.cpp \
	$(USER_DIR)/whistory.h $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPTIX_FLAGS) $(CUDPP_FLAGS) $(CUDA_FLAGS) $(PYTHON_FLAGS) -c -g $(USER_DIR)/whistory_test.cpp

gtest-all.o: $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c -g \
            $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c -g \
            $(GTEST_DIR)/src/gtest_main.cc

gtest.a : gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

gtest_main.a : gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^
#/google test

libwarp.so: $(ptx_objects) $(COBJS)
	$(NVCC) --shared $(NVCC_FLAGS) $(OPTIX_FLAGS) $(CUDPP_FLAGS) $(PYTHON_FLAGS) $(PNG_FLAGS) $(CURAND_LIBS) $(OPTIX_LIBS) $(CUDPP_LIBS) $(PYTHON_LIBS) $(PNG_LIBS) $(COBJS)  -o libwarp.so

gpu: libwarp.so
	$(NVCC) $(NVCC_FLAGS) $(OPTIX_FLAGS) $(CUDPP_FLAGS) $(PYTHON_FLAGS)  -L/home/krowland/warp/warp -lwarp -loptix main.cpp -o $@

debug: libwarp.so
	$(NVCC) $(NVCC_FLAGS) $(OPTIX_FLAGS) $(CUDPP_FLAGS) $(PYTHON_FLAGS)  -L/home/krowland/warp/warp -lwarp main.cpp -o $@

optixtest: libwarp.so
	$(NVCC) -m64 $(OPTIX_FLAGS) $(CUDPP_FLAGS) $(PYTHON_FLAGS)  -L/home/krowland/warp/warp -lwarp -loptix optixtest.cpp -o $@

python: libwarp.so
	swig -python -c++ warp.i;   \
        $(CXX) -fPIC -c $(PYTHON_FLAGS) $(CUDPP_FLAGS) $(CUDA_FLAGS) warp_wrap.cxx;  \
        $(CXX) -shared libwarp.so warp_wrap.o -o _warp.so $(PYTHON_FLAGS) -lpython2.7

#google test
optix_stuff_test : libwarp.so optix_stuff.o optix_stuff_test.o gtest_main.a
	$(CXX) $(PNG_FLAGS) $(PNG_LIBS) $(OPTIX_FLAGS) $(OPTIX_LIBS) $(CUDPP_FLAGS) $(PYTHON_FLAGS) $(CPPFLAGS) $(CXXFLAGS) -pthread $^ -o $@

primitive_test : libwarp.so primitive.o primitive_test.o gtest_main.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -pthread $^ -o $@

wgeometry_test : libwarp.so wgeometry.o wgeometry_test.o gtest_main.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -pthread $^ -o $@

whistory_test : libwarp.so whistory.o whistory_test.o gtest_main.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(CUDA_FLAGS) $(CURAND_LIBS) $(OPTIX_FLAGS) $(OPTIX_LIBS) $(CUDPP_FLAGS) $(CUDPP_LIBS) $(PYTHON_FLAGS) $(PYTHON_LIBS) -lcudart -pthread $^ -o $@
#/google test
