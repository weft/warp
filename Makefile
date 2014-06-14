#
CC =  gcc
CXX = g++
OPTIX = /usr/local/OptiX-3.0.1
NVCC = nvcc 
ARCH = -arch sm_30
C_FLAGS = -O3 -m64 -fPIC
NVCC_FLAGS = -m64  -use_fast_math --compiler-options '-fPIC'
CURAND_LIBS = -lcurand
OPTIX_FLAGS = -I$(OPTIX)/include -L$(OPTIX)/lib64 
OPTIX_LIBS = -loptix 
CUDA_FLAGS = -I/usr/local/cuda/include -L/usr/local/cuda/lib
CUDPP_PATH = /usr/local/cudpp-2.1/
CUDPP_FLAGS = -I/$(CUDPP_PATH)/include -L/$(CUDPP_PATH)/lib
CUDPP_LIBS = -lcudpp_hash -lcudpp
PYTHON_FLAGS = -I/usr/local/anaconda/include/python2.7 -L/usr/local/anaconda/lib/
PYTHON_LIBS = -lpython2.7
PNG_FLAGS = -L/
PNG_LIBS = -lpng15
LIBS =


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
		reaction_edges.o \
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
		libwarp.so 

clean:
	rm -f *.ptx *.o *.so gpu debug optixtest

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

reaction_edges.o:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) -c reaction_edges.cu

write_to_file.o:
	$(NVCC) $(ARCH) $(NVCC_FLAGS) -c write_to_file.cu

whistory.o:
	$(CXX) $(C_FLAGS)  $(OPTIX_FLAGS) $(CUDPP_FLAGS) $(PNG_FLAGS) $(PYTHON_FLAGS) $(CUDA_FLAGS)  -c whistory.cpp

wgeometry.o:
	$(CXX) $(C_FLAGS)  -c wgeometry.cpp

primitive.o:
	$(CXX) $(C_FLAGS)  -c primitive.cpp

optix_stuff.o: device_copies.o
	$(CXX) -m64 -fPIC  $(OPTIX_FLAGS) $(CUDA_FLAGS) $(PNG_FLAGS) -c optix_stuff.cpp

libwarp.so: $(ptx_objects) $(COBJS)
	$(NVCC) --shared $(NVCC_FLAGS) $(OPTIX_FLAGS) $(CUDPP_FLAGS) $(PYTHON_FLAGS) $(PNG_FLAGS) $(CURAND_LIBS) $(OPTIX_LIBS) $(CUDPP_LIBS) $(PYTHON_LIBS) $(PNG_LIBS) $(COBJS)  -o libwarp.so

gpu: libwarp.so
	$(NVCC) -m64 $(OPTIX_FLAGS) $(CUDPP_FLAGS) $(PYTHON_FLAGS)  -L/mnt/scratch/gpu-cpp -lwarp main.cpp -o $@

optixtest: libwarp.so
	$(NVCC) -m64 $(OPTIX_FLAGS) $(CUDPP_FLAGS) $(PYTHON_FLAGS)  -L/Users/rmb/code/gpu-cpp -lwarp -loptix optixtest.cpp -o $@

warp.py: libwarp.so
	swig warp.i
