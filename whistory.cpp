#include <vector> 
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <cmath>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <cuda.h>
#include <curand.h>
#include <cudpp.h>
#include <Python.h>
#include <png++/png.hpp>
#include "datadef.h"
#include "wprimitive.h"
#include "wgeometry.h"
#include "optix_stuff.h"
#include "warp_cuda.h"
#include "whistory.h"

// CUDA error check wrapper
inline void check_cuda(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}
#define check_cuda(ans) { check_cuda((ans), __FILE__, __LINE__); }

// instantiate single, global optix object
optix_stuff optix_obj;

// whistory class
whistory::~whistory(){
}
void whistory::print_banner(){
	using namespace std;
	cout << "\n"; 
	cout << "\e[1;32m" << "-- Weaving All the Random Particles ---" << "\e[m \n";
	cout << "\e[1;31m";
	cout << "___       ________________________ " 	<< "\n";
	cout << "__ |     / /__    |__  __ \\__  __ \\"	<< "\n";
	cout << "__ | /| / /__  /| |_  /_/ /_  /_/ /" 	<< "\n";
	cout << "__ |/ |/ / _  ___ |  _, _/_  ____/ " 	<< "\n";
	cout << "____/|__/  /_/  |_/_/ |_| /_/      " 	<< "\n";
	cout << "\e[m \n";
	cout << "\e[1;32m" << "Monte Carlo Neutron Transport at WARP speed!  Err... hopefully :)" << "\e[m \n";
	cout << "\n";
}
whistory::whistory(unsigned Nin, wgeometry problem_geom_in){
	// do problem gemetry stuff first
	problem_geom = problem_geom_in;
	// set run flag
	RUN_FLAG = 1;
	// keff stuff
	keff_sum = 0.0;
	keff2_sum = 0.0;
	keff_err = 0.0;
	// device data stuff
	N = Nin;
	// 3 should be more than enough for criticality, might not be for fixed
	Ndataset = Nin * 1;
	reduced_yields_total = 0;
	reduced_weight_total = 0;
	accel_type = "Sbvh";
	// default device to 0
	compute_device = 0;
	// not initialized
	is_initialized = 0;
	// prints
	print_flag = 1;
	dump_flag = 0;
}
void whistory::init(){

	//
	// WELCOME!
	//
	print_banner();

	//
	// EARLY DEVICE SETUP
	//
	// device must be set BEFORE an CUDA creation (streams included)
	check_cuda(cudaSetDevice(compute_device));
	//clear device
	check_cuda(cudaDeviceReset());
	// create streams
	for (int i = 0; i < 4; ++i){
		check_cuda(cudaStreamCreate(&stream[i]));
	}

	//
	// OPTIX 
	//
	optix_obj.N=Ndataset;
	optix_obj.stack_size_multiplier=1;
	optix_obj.init(problem_geom,compute_device,accel_type);
	if(print_flag >= 1){
		optix_obj.print();
	}
	check_cuda(cudaPeekAtLastError());

	//
	// CUDA 
	//
	// set device, need to set again after optix sometimes... 
	check_cuda(cudaSetDevice(compute_device));
	if(print_flag >= 1){
		std::cout << "\e[1;32m" << "Dataset size is "<< N << "\e[m \n";
	}
	// block/thread structure
	NUM_THREADS		= 256;
	blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;
	if(print_flag >= 1){
		std::cout << "\e[1;32m" << "Compute device set to "<< compute_device << "\e[m \n";
	}
	if(print_flag >= 2){
		device_report();
	}

	// set print buffer size
	check_cuda(cudaDeviceSetLimit(cudaLimitPrintfFifoSize, (size_t) 10*1048576 ));

	// fission points array
	fiss_img  =  new long unsigned [300*300]; for(int i=0;i<300*300;i++){fiss_img[i]=0;}
	// init counters to 0
	total_bytes_scatter = 0;
	total_bytes_energy  = 0;

	//copy any info needed from the geometry read
	memcpy(outer_cell_dims,optix_obj.outer_cell_dims,6*sizeof(float));
	outer_cell 		= optix_obj.get_outer_cell();
	outer_cell_type = optix_obj.get_outer_cell_type();
	isotopes		= problem_geom.isotopes;
	n_isotopes 		= problem_geom.n_isotopes;
	n_materials		= problem_geom.n_materials;

	// map edge array and reduced values
	n_edges = 11;
	check_cuda(cudaHostAlloc( &edges          , n_edges*sizeof(unsigned), cudaHostAllocMapped));
	check_cuda(cudaHostAlloc( &reduced_yields ,       1*sizeof(unsigned), cudaHostAllocMapped));
	check_cuda(cudaHostAlloc( &reduced_weight ,       1*sizeof(float)   , cudaHostAllocMapped));
	check_cuda(cudaHostGetDevicePointer(&d_edges          ,  edges          , 0));
	check_cuda(cudaHostGetDevicePointer(&d_reduced_yields ,  reduced_yields , 0));
	check_cuda(cudaHostGetDevicePointer(&d_reduced_weight ,  reduced_weight , 0));

	// set anything else
	filename = "warp";

	//
	// init/copy the rest not done previously or in OptiX
	//
	init_host();
	init_device();
	init_RNG();
	init_CUDPP();
	init_cross_sections();

	// set inititalized flag, print done if flagged
	is_initialized = 1;
	if(print_flag >= 2){
		printf("\e[1;31mDone with init\e[m\n");
	}
}
void whistory::init_host(){

	// init tally arrays
	dh_tally	=	new tally_data[n_tallies];
	h_tally		=	new tally_data[n_tallies];

	// allocate
	h_particles.space	= new spatial_data 	[Ndataset];
	h_particles.E		= new float 		[Ndataset];
	h_particles.Q		= new float 		[Ndataset];
	h_particles.rn_bank	= new unsigned 		[Ndataset];
	h_particles.index	= new unsigned 		[Ndataset];
	h_particles.cellnum	= new unsigned 		[Ndataset];
	h_particles.matnum	= new unsigned 		[Ndataset];
	h_particles.rxn		= new unsigned 		[Ndataset];
	h_particles.isonum	= new unsigned 		[Ndataset];
	h_particles.yield	= new unsigned 		[Ndataset];
	h_particles.weight	= new float  		[Ndataset];
	remap				= new unsigned 		[Ndataset];
	zeros				= new unsigned 		[Ndataset];
	ones				= new unsigned 		[Ndataset];
	fones				= new float  		[Ndataset];

	//  init
	for(int k=0;k<Ndataset;k++){
		h_particles.space[k].x				= 0.0;
		h_particles.space[k].y				= 0.0;
		h_particles.space[k].z				= 0.0;
		h_particles.space[k].xhat			= 0.0;
		h_particles.space[k].yhat			= 0.0;
		h_particles.space[k].zhat			= 0.0;
		h_particles.space[k].surf_dist		= 10000.0;
		h_particles.space[k].macro_t		= 0.0;
		h_particles.space[k].enforce_BC		= 0;
		h_particles.space[k].norm[0]		= 0;
		h_particles.space[k].norm[1]		= 0;
		h_particles.space[k].norm[2]		= 0;
		h_particles.E[k]					= 2.5;
		h_particles.Q[k]					= 0.0;
		h_particles.cellnum[k]				= 0;
		h_particles.matnum[k]				= 0;
		h_particles.rxn[k]					= 0;
		h_particles.isonum[k]				= 0;
		h_particles.yield[k]				= 0;
		h_particles.weight[k]				= 1;
		remap[k]							= k;
		zeros[k]							= 0;
		ones[k]								= 1;
		fones[k]							= 1.0;
	}

	// set host tally size, bounds, allocate, and zero out all tallies
	for(int i=0;i<n_tallies;i++){
		h_tally[i].cell			=	0;
		h_tally[i].length		=	1024;
		h_tally[i].E_min		=	1e-12;
		h_tally[i].E_max		=	20;
		h_tally[i].score		=	new float			[h_tally[i].length];
		h_tally[i].square		=	new float			[h_tally[i].length];
		h_tally[i].count		=	new unsigned		[h_tally[i].length];
		h_tally[i].score_total	=	new double			[h_tally[i].length];
		h_tally[i].square_total	=	new double			[h_tally[i].length];
		h_tally[i].count_total	=	new long unsigned	[h_tally[i].length];
		for(int k=0;k<h_tally[k].length;k++){
			h_tally[i].score[       k]	=	0.0;
			h_tally[i].square[      k]	=	0.0;
			h_tally[i].count[       k]	=	0.0;
			h_tally[i].score_total[ k]	=	0.0;
			h_tally[i].square_total[k]	=	0.0;
			h_tally[i].count_total[ k]	=	0.0;
		}
	}

}
void whistory::init_device(){

	// copy pointers initialized by optix
	dh_particles.space		= (spatial_data*)	optix_obj.positions_ptr;
	dh_particles.cellnum	= (unsigned*)		optix_obj.cellnum_ptr;
	dh_particles.matnum		= (unsigned*)		optix_obj.matnum_ptr;
	dh_particles.rxn		= (unsigned*)		optix_obj.rxn_ptr;
	d_remap					= (unsigned*)		optix_obj.remap_ptr;

	// init others only used on CUDA side
	cudaMalloc( &d_tally				,n_tallies*sizeof(tally_data)			);
	cudaMalloc( &d_particles			,		 1*sizeof(particle_data)		);
	cudaMalloc( &d_xsdata				,        1*sizeof(cross_section_data) 	);
	cudaMalloc( &dh_particles.E			, Ndataset*sizeof(float)				);
	cudaMalloc( &dh_particles.Q			, Ndataset*sizeof(float)				);
	cudaMalloc( &dh_particles.rn_bank	, Ndataset*sizeof(float)				);
	cudaMalloc( &dh_particles.isonum	, Ndataset*sizeof(unsigned)				);
	cudaMalloc( &dh_particles.yield		, Ndataset*sizeof(unsigned)				);
	cudaMalloc( &dh_particles.weight	, Ndataset*sizeof(float)				);
	cudaMalloc( &dh_particles.index		, Ndataset*sizeof(unsigned)				);
	cudaMalloc( &d_valid_result			, Ndataset*sizeof(unsigned)				);
	cudaMalloc( &d_valid_N				,        1*sizeof(unsigned)				);
	cudaMalloc( &d_fissile_points		, Ndataset*sizeof(spatial_data)			);
	cudaMalloc( &d_fissile_energy       , Ndataset*sizeof(float)				);
	cudaMalloc( &d_scanned 				, Ndataset*sizeof(unsigned)				);
	cudaMalloc( &d_num_completed 		,        1*sizeof(unsigned)				);
	cudaMalloc( &d_num_active 			,        1*sizeof(unsigned)				);
	cudaMalloc( &d_zeros				, Ndataset*sizeof(unsigned)				);

	// copy values from initialized host arrays
	cudaMemcpy( dh_particles.space		, h_particles.space			, Ndataset*sizeof(spatial_data)	, cudaMemcpyHostToDevice );
	cudaMemcpy( dh_particles.cellnum	, h_particles.cellnum		, Ndataset*sizeof(unsigned)		, cudaMemcpyHostToDevice );
	cudaMemcpy( dh_particles.matnum		, h_particles.matnum		, Ndataset*sizeof(unsigned)		, cudaMemcpyHostToDevice );
	cudaMemcpy( dh_particles.rxn		, h_particles.rxn			, Ndataset*sizeof(unsigned)		, cudaMemcpyHostToDevice );
	cudaMemcpy( dh_particles.E			, h_particles.E				, Ndataset*sizeof(float)		, cudaMemcpyHostToDevice );
	cudaMemcpy( dh_particles.Q			, h_particles.Q				, Ndataset*sizeof(float)		, cudaMemcpyHostToDevice );
	cudaMemcpy( dh_particles.rn_bank	, h_particles.rn_bank		, Ndataset*sizeof(float)		, cudaMemcpyHostToDevice );
	cudaMemcpy( dh_particles.isonum		, h_particles.isonum		, Ndataset*sizeof(unsigned)		, cudaMemcpyHostToDevice );
	cudaMemcpy( dh_particles.yield		, h_particles.yield			, Ndataset*sizeof(unsigned)		, cudaMemcpyHostToDevice );
	cudaMemcpy( dh_particles.weight		, h_particles.weight		, Ndataset*sizeof(float)		, cudaMemcpyHostToDevice );
	cudaMemcpy( dh_particles.index		, h_particles.index			, Ndataset*sizeof(unsigned)		, cudaMemcpyHostToDevice );
	cudaMemcpy( d_valid_result			, zeros						, Ndataset*sizeof(unsigned)		, cudaMemcpyHostToDevice );
	cudaMemcpy( d_remap					, remap						, Ndataset*sizeof(unsigned)		, cudaMemcpyHostToDevice );
	cudaMemcpy( d_zeros					, zeros						, Ndataset*sizeof(unsigned)		, cudaMemcpyHostToDevice );

	// init tally containers
	for( int i=0 ; i<n_tallies ; i++ ){
		// allocate
		cudaMalloc( &dh_tally[i].score 			, h_tally[i].length*sizeof(float*)        );
		cudaMalloc( &dh_tally[i].square 		, h_tally[i].length*sizeof(float*)        );
		cudaMalloc( &dh_tally[i].count 			, h_tally[i].length*sizeof(unsigned*)     );
		cudaMalloc( &dh_tally[i].score_total 	, h_tally[i].length*sizeof(double*)       );
		cudaMalloc( &dh_tally[i].square_total 	, h_tally[i].length*sizeof(double*)       );
		cudaMalloc( &dh_tally[i].count_total 	, h_tally[i].length*sizeof(long unsigned*));
		// copy initialized values
		cudaMemcpy(  dh_tally[i].score	    	, h_tally[i].score	       , h_tally[i].length*sizeof(float)			,  cudaMemcpyHostToDevice );
		cudaMemcpy(  dh_tally[i].square	    	, h_tally[i].square	       , h_tally[i].length*sizeof(float)			,  cudaMemcpyHostToDevice );
		cudaMemcpy(  dh_tally[i].count	    	, h_tally[i].count	       , h_tally[i].length*sizeof(unsigned)			,  cudaMemcpyHostToDevice );
		cudaMemcpy(  dh_tally[i].score_total  	, h_tally[i].score_total   , h_tally[i].length*sizeof(double)			,  cudaMemcpyHostToDevice );
		cudaMemcpy(  dh_tally[i].square_total 	, h_tally[i].square_total  , h_tally[i].length*sizeof(double)			,  cudaMemcpyHostToDevice );
		cudaMemcpy(  dh_tally[i].count_total  	, h_tally[i].count_total   , h_tally[i].length*sizeof(long unsigned)	,  cudaMemcpyHostToDevice );
	}

	// copy host structures (containing the device pointers) to the device structure
	cudaMemcpy(  d_particles, 	&dh_particles	,		  1*sizeof(particle_data),		cudaMemcpyHostToDevice);
	cudaMemcpy(  d_tally,		dh_tally		, n_tallies*sizeof(tally_data),			cudaMemcpyHostToDevice);

	// check errors
	check_cuda(cudaPeekAtLastError());

}
void whistory::init_RNG(){
	unsigned seed = time( NULL );
	if(print_flag >= 2){
		std::cout << "\e[1;32m" << "Initializing random number bank on device using MTGP32 with seed of " << seed << "..." << "\e[m \n";
	}
	curandCreateGenerator( &rand_gen , CURAND_RNG_PSEUDO_MTGP32 );  //mersenne twister type
	curandSetPseudoRandomGeneratorSeed( rand_gen , 123456789ULL );
	curandGenerate( rand_gen , dh_particles.rn_bank , Ndataset  );
	cudaMemcpy(h_particles.rn_bank , dh_particles.rn_bank , Ndataset*sizeof(unsigned) , cudaMemcpyDeviceToHost); // copy bank back to keep seeds
	check_cuda(cudaPeekAtLastError());
}
void whistory::update_RNG(){

	curandGenerate( rand_gen , dh_particles.rn_bank , Ndataset );
	check_cuda(cudaPeekAtLastError());

}
void whistory::init_CUDPP(){
	
	if(print_flag >= 2){
		std::cout << "\e[1;32m" << "Initializing CUDPP..." << "\e[m \n";
	}
	// global objects
	res = cudppCreate(&theCudpp);
	if (res != CUDPP_SUCCESS){fprintf(stderr, "Error initializing CUDPP Library.\n");}
	
	if(print_flag >= 2){
		std::cout << "  configuring compact..." << "\n";
	}
	// sort stuff
	compact_config.op = CUDPP_ADD;
	compact_config.datatype = CUDPP_INT;
	compact_config.algorithm = CUDPP_COMPACT;
	compact_config.options = CUDPP_OPTION_FORWARD;
	res = cudppPlan(theCudpp, &compactplan, compact_config, Ndataset, 1, 0);
	if (CUDPP_SUCCESS != res){printf("Error creating CUDPPPlan for compact\n");exit(-1);}

	if(print_flag >= 2){
		std::cout << "  configuring reduction..." << "\n";
	}
	// int reduction stuff
	redu_int_config.op = CUDPP_ADD;
	redu_int_config.datatype = CUDPP_INT;
	redu_int_config.algorithm = CUDPP_REDUCE;
	redu_int_config.options = 0;
	res = cudppPlan(theCudpp, &reduplan_int, redu_int_config, Ndataset, 1, 0);
	if (CUDPP_SUCCESS != res){printf("Error creating CUDPPPlan for reduction\n");exit(-1);}
	
	// float reduction stuff
	redu_float_config.op = CUDPP_ADD;
	redu_float_config.datatype = CUDPP_FLOAT;
	redu_float_config.algorithm = CUDPP_REDUCE;
	redu_float_config.options = 0;
	res = cudppPlan(theCudpp, &reduplan_float, redu_float_config, Ndataset, 1, 0);
	if (CUDPP_SUCCESS != res){printf("Error creating CUDPPPlan for reduction\n");exit(-1);}

	if(print_flag >= 2){
		std::cout << "  configuring scan..." << "\n";
	}
	// int reduction stuff
	scan_int_config.op = CUDPP_ADD;
	scan_int_config.datatype = CUDPP_INT;
	scan_int_config.algorithm = CUDPP_SCAN;
	scan_int_config.options = CUDPP_OPTION_EXCLUSIVE;
	res = cudppPlan(theCudpp, &scanplan_int, scan_int_config, Ndataset, 1, 0);
	if (CUDPP_SUCCESS != res){printf("Error creating CUDPPPlan for scan\n");exit(-1);}

	if(print_flag >= 2){
		std::cout << "  configuring radix sort..." << "\n";
	}
	// int reduction stuff
	radix_config.algorithm = CUDPP_SORT_RADIX;
	radix_config.datatype = CUDPP_UINT;
	radix_config.options = CUDPP_OPTION_KEY_VALUE_PAIRS;
	res = cudppPlan(theCudpp, &radixplan, radix_config, Ndataset, 1, 0);
	if (CUDPP_SUCCESS != res){printf("Error creating CUDPPPlan for radix sort\n");exit(-1);}

	check_cuda(cudaPeekAtLastError());

}
unsigned whistory::reduce_yield(){

	unsigned reduced_yields;

	res = cudppReduce(reduplan_int, d_reduced_yields, h_particles.yield, N);
	if (res != CUDPP_SUCCESS){fprintf(stderr, "Error in reducing yield values\n");exit(-1);}
	cudaMemcpy(&reduced_yields, d_reduced_yields, 1*sizeof(unsigned), cudaMemcpyDeviceToHost);

	return reduced_yields;

}
float whistory::reduce_weight(){

	float reduced_weight;

	res = cudppReduce(reduplan_float, d_reduced_weight, h_particles.weight, N);
	if (res != CUDPP_SUCCESS){fprintf(stderr, "Error in reducing weight values\n");exit(-1);}
	cudaMemcpy(&reduced_weight, d_reduced_weight, 1*sizeof(float), cudaMemcpyDeviceToHost);

	return reduced_weight;

}
void whistory::accumulate_keff(unsigned converged, unsigned iteration, double* keff, float* keff_cycle){

	float this_count, this_square, this_mean, keff_err2;

	long unsigned reduced_yields = reduce_yield();
	double        reduced_weight  = reduce_weight();

	*keff_cycle = reduced_yields / reduced_weight;

	if(converged){
		reduced_yields_total += reduced_yields;
		reduced_weight_total += reduced_weight;
		*keff       = reduced_yields_total / reduced_weight_total;
		keff_sum   += *keff_cycle; 
		keff2_sum  += (*keff_cycle)*(*keff_cycle);
		this_count  = iteration + 1; 
		this_square = keff2_sum ;
		this_mean   = keff_sum  ;
		//printf("this_mean %10.8E this_square %10.8E iteration %d\n",this_mean,this_square, iteration+1);
		keff_err2   = (1.0/((this_count - 1.0))) * ( (this_count*this_square)/(this_mean*this_mean) -  1.0 ) ;
		if(keff_err2>0.0){ 	
			keff_err = sqrtf(keff_err2);}
		else{					
			keff_err = 0.0;}
	}

	//printf("reduced_total %lu  reduced %u keff %10.8E keff_cycle %10.8E iteration %u\n",reduced_yields_total,reduced_yields,*keff,*keff_cycle,iteration+1);

}
void whistory::accumulate_tally(){

//	//copy data to host
//	cudaMemcpy(tally_score,  d_tally_score,  n_tally*sizeof(float),    cudaMemcpyDeviceToHost);
//	cudaMemcpy(tally_square, d_tally_square, n_tally*sizeof(float),    cudaMemcpyDeviceToHost);
//	cudaMemcpy(tally_count,  d_tally_count,  n_tally*sizeof(unsigned), cudaMemcpyDeviceToHost);
//
//	//zero out vectors
//	cudaMemcpy(d_tally_score,  d_zeros, n_tally*sizeof(float),    cudaMemcpyDeviceToDevice);
//	cudaMemcpy(d_tally_square, d_zeros, n_tally*sizeof(float),    cudaMemcpyDeviceToDevice);
//	cudaMemcpy(d_tally_count,  d_zeros, n_tally*sizeof(unsigned), cudaMemcpyDeviceToDevice);
//
//	//perform sums on 64bit host side values
//	for(unsigned k=0 ; k<n_tally ; k++){
//		tally_score_total[k] 	+=  tally_score[k];
//		tally_square_total[k]	+=  tally_square[k];
//		tally_count_total[k] 	+=  tally_count[k];
//		//printf("score %10.8E %10.8E\n",tally_score[k],tally_score_total[k]);
//	}

	
}
void whistory::copy_python_buffer(float** device_pointer,float** host_pointer,std::string function_name){
	// version for two float pointers

	// misc variables for maths
	unsigned bytes,rows,columns,ndim;

	// python variables
	PyObject *pName, *pModule, *pDict, *pFunc;
	PyObject *pArgs, *pValue, *pString, *pBuffObj, *pObjList;
	PyObject *call_result;
	PyObject *call_string,*arg_string;
	PyObject *pClass;
	Py_buffer pBuff;

	// get the MT array buffer from Python
	call_string = PyString_FromString(function_name.c_str());
	call_result = PyObject_CallMethodObjArgs(xsdat_instance, call_string, NULL);
	Py_DECREF(call_string);
	if (PyObject_CheckBuffer(call_result)){
		PyObject_GetBuffer(call_result, &pBuff,PyBUF_ND);
	}
	else{
		PyErr_Print();
	    fprintf(stderr, "Returned object does not support buffer interface\n");
	    return;
	}

	// get array size info
	get_Py_buffer_dims(&rows,&columns,&bytes,&pBuff);

	// allocate xs_data pointer arrays
	*host_pointer = new float [rows*columns];
	// check to make sure bytes *= elements
	assert(bytes==rows*columns*sizeof(float));
	// copy python buffer contents to pointer
	memcpy( *host_pointer,   pBuff.buf , bytes );
	// allocate device memory now that we know the size!
	cudaMalloc(device_pointer,bytes);
	// copy the xs data to device
	cudaMemcpy(*device_pointer, *host_pointer, bytes, cudaMemcpyHostToDevice);
	// release python variable to free memory
	Py_DECREF(call_result);

}
void whistory::copy_python_buffer(unsigned** device_pointer,unsigned** host_pointer,std::string function_name){
	// version for two unsigned pointers

	// misc variables for maths
	unsigned bytes,rows,columns;

	// python variables
	PyObject *pName, *pModule, *pDict, *pFunc;
	PyObject *pArgs, *pValue, *pString, *pBuffObj, *pObjList;
	PyObject *call_result;
	PyObject *call_string,*arg_string;
	PyObject *pClass;
	Py_buffer pBuff;

	// get the MT array buffer from Python
	call_string = PyString_FromString(function_name.c_str());
	call_result = PyObject_CallMethodObjArgs(xsdat_instance, call_string, NULL);
	Py_DECREF(call_string);
	if (PyObject_CheckBuffer(call_result)){
		PyObject_GetBuffer(call_result, &pBuff,PyBUF_ND);
	}
	else{
		PyErr_Print();
	    fprintf(stderr, "Returned object does not support buffer interface\n");
	    return;
	}

	// get array size info
	get_Py_buffer_dims(&rows,&columns,&bytes,&pBuff);

	// allocate xs_data pointer arrays
	*host_pointer = new unsigned [rows*columns];
	// check to make sure bytes *= elements
	assert(bytes==rows*columns*sizeof(unsigned));
	// copy python buffer contents to pointer
	memcpy( *host_pointer,   pBuff.buf , bytes );
	// allocate device memory now that we know the size!
	cudaMalloc(device_pointer,bytes);
	// copy the xs data to device
	cudaMemcpy(*device_pointer, *host_pointer, bytes, cudaMemcpyHostToDevice);
	// release python variable to free memory
	Py_DECREF(call_result);

}
void whistory::copy_python_buffer(float** host_pointer,std::string function_name){
	// version for one (host) float pointer

	// misc variables for maths
	unsigned bytes,rows,columns;

	// python variables
	PyObject *pName, *pModule, *pDict, *pFunc;
	PyObject *pArgs, *pValue, *pString, *pBuffObj, *pObjList;
	PyObject *call_result;
	PyObject *call_string,*arg_string;
	PyObject *pClass;
	Py_buffer pBuff;

	// get the MT array buffer from Python
	call_string = PyString_FromString(function_name.c_str());
	call_result = PyObject_CallMethodObjArgs(xsdat_instance, call_string, NULL);
	Py_DECREF(call_string);
	if (PyObject_CheckBuffer(call_result)){
		PyObject_GetBuffer(call_result, &pBuff,PyBUF_ND);
	}
	else{
		PyErr_Print();
	    fprintf(stderr, "Returned object does not support buffer interface\n");
	    return;
	}

	// get array size info
	get_Py_buffer_dims(&rows,&columns,&bytes,&pBuff);

	// allocate xs_data pointer arrays
	*host_pointer = new float [rows*columns];
	// check to make sure bytes *= elements
	assert(bytes==rows*columns*sizeof(float));
	// copy python buffer contents to pointer
	memcpy( *host_pointer,   pBuff.buf , bytes );
	// release python variable to free memory
	Py_DECREF(call_result);

}
void whistory::copy_python_buffer(unsigned** host_pointer,std::string function_name){
	// version for one (host) unsigned pointer

	// misc variables for maths
	unsigned bytes,rows,columns;

	// python variables
	PyObject *pName, *pModule, *pDict, *pFunc;
	PyObject *pArgs, *pValue, *pString, *pBuffObj, *pObjList;
	PyObject *call_result;
	PyObject *call_string,*arg_string;
	PyObject *pClass;
	Py_buffer pBuff;

	// get the MT array buffer from Python
	call_string = PyString_FromString(function_name.c_str());
	call_result = PyObject_CallMethodObjArgs(xsdat_instance, call_string, NULL);
	Py_DECREF(call_string);
	if (PyObject_CheckBuffer(call_result)){
		PyObject_GetBuffer(call_result, &pBuff,PyBUF_ND);
	}
	else{
		PyErr_Print();
	    fprintf(stderr, "Returned object does not support buffer interface\n");
	    return;
	}

	// get array size info
	get_Py_buffer_dims(&rows,&columns,&bytes,&pBuff);

	// allocate xs_data pointer arrays
	*host_pointer = new unsigned [rows*columns];
	// check to make sure bytes *= elements
	assert(bytes==rows*columns*sizeof(unsigned));
	// copy python buffer contents to pointer
	memcpy( *host_pointer,   pBuff.buf , bytes );
	// release python variable to free memory
	Py_DECREF(call_result);

}
int whistory::init_python(){

	// misc variables
	int do_final;

	// python variables
	PyObject *pName, *pModule, *pDict, *pFunc;
	PyObject *pArgs, *pValue, *pString, *pBuffObj, *pObjList;
	PyObject *call_result;
	PyObject *call_string,*arg_string;
	PyObject *pClass;
	Py_buffer pBuff;

	if (Py_IsInitialized()){
		if(print_flag >= 2){
			printf("Python interpreter already initialized\n");
		}
		do_final=0;
	}
	else{
		Py_Initialize();
		do_final=1;
	}

	pName = PyString_FromString("unionize");
	pModule = PyImport_Import(pName);
	Py_DECREF(pName);

	if ( pModule==NULL ){
		PyErr_Print();
		fprintf(stderr, "Failed to import \"%s\"\n", "unionize");
		return 0;	
	}

	pName = PyString_FromString("cross_section_data");
	xsdat_instance = PyObject_CallMethodObjArgs(pModule,pName,NULL);
	PyErr_Print();
	Py_DECREF(pName);


	if (xsdat_instance != NULL) {

		// init the libraries wanted
		for (int i=0; i<n_isotopes; i++){
			call_string = PyString_FromString("_add_isotope");
			arg_string  = PyString_FromString(isotopes[i].c_str());
			call_result = PyObject_CallMethodObjArgs(xsdat_instance, call_string, arg_string, NULL);
			PyErr_Print();
			Py_DECREF(arg_string);
			Py_DECREF(call_string);
			Py_DECREF(call_result);
		}
	
		// read the tables
		call_string = PyString_FromString("_read_tables");
		arg_string  = PyString_FromString(problem_geom.datapath.c_str());
		call_result = PyObject_CallMethodObjArgs(xsdat_instance, call_string, arg_string, NULL);
		PyErr_Print();
		Py_DECREF(arg_string);
		Py_DECREF(call_string);
		Py_DECREF(call_result);

		// unionize the main energy grid across all isotopes
		call_string = PyString_FromString("_unionize");
		call_result = PyObject_CallMethodObjArgs(xsdat_instance, call_string, NULL);
		PyErr_Print();
		Py_DECREF(call_string);
		Py_DECREF(call_result);

		// make the total MT reaction list from all isotopes
		call_string = PyString_FromString("_insert_reactions");
		call_result = PyObject_CallMethodObjArgs(xsdat_instance, call_string, NULL);
		PyErr_Print();
		Py_DECREF(call_string);
		Py_DECREF(call_result);

		// allocate the unionized array
		call_string = PyString_FromString("_allocate_arrays");
		call_result = PyObject_CallMethodObjArgs(xsdat_instance, call_string, NULL);
		PyErr_Print();
		Py_DECREF(call_string);
		Py_DECREF(call_result);

		// insert and interpolate the cross sections
		call_string = PyString_FromString("_interpolate");
		call_result = PyObject_CallMethodObjArgs(xsdat_instance, call_string, NULL);
		PyErr_Print();
		Py_DECREF(call_string);
		Py_DECREF(call_result);

	}
	else {
		PyErr_Print();
		fprintf(stderr, "Failed to instanciate \"%s\"\n", "unionize.cross_section_data");
	}

	return do_final;

}
void whistory::get_Py_buffer_dims(unsigned* rows, unsigned* columns, unsigned* bytes, Py_buffer* pBuff){

	// get size
	bytes[0]	=	pBuff[0].len;

	// get dimensions
	if (pBuff[0].ndim == 0){
		rows[0]		=	1;
		columns[0]	=	1;
	}
	else if (pBuff[0].ndim == 1){
		rows[0]		=	pBuff[0].shape[0];
		columns[0]	=	1;
	}
	else if (pBuff[0].ndim == 2){
		rows[0]		=	pBuff[0].shape[0];
		columns[0]	=	pBuff[0].shape[1];
	}
	else {
		fprintf(stderr, "Python buffer returned array with dimension > 2!!!\n");
	}

}
void whistory::copy_scatter_data(){
	// get scattering ditributions from PyNE and set up data heirarchy on host and device 

	if(print_flag >= 2){
		std::cout << "\e[1;32m" << "Loading/copying scattering data and linking..." << "\e[m \n";
	}

	// python variables for arguments
	PyObject	*row_obj, *col_obj, *call_string, *obj_list;
	PyObject	*lower_erg_obj, *lower_len_obj, *lower_law_obj, *lower_intt_obj, *lower_var_obj, *lower_pdf_obj, *lower_cdf_obj; 
	PyObject	*upper_erg_obj, *upper_len_obj, *upper_law_obj, *upper_intt_obj, *upper_var_obj, *upper_pdf_obj, *upper_cdf_obj, *next_dex_obj; 
	Py_buffer	lower_var_buff, lower_pdf_buff, lower_cdf_buff, upper_var_buff, upper_pdf_buff, upper_cdf_buff;
	
	// local temp variables
	int				row, col;
	float			nu_t, nu_p;
	int				this_isotope;
	unsigned		array_index, next_dex, lower_rows;
	unsigned		lower_var_buff_bytes,	lower_pdf_buff_bytes,	lower_cdf_buff_bytes;
	unsigned		lower_var_buff_rows,	lower_pdf_buff_rows,	lower_cdf_buff_rows;
	unsigned		lower_var_buff_columns,	lower_pdf_buff_columns,	lower_cdf_buff_columns;
	unsigned		upper_var_buff_bytes,	upper_pdf_buff_bytes,	upper_cdf_buff_bytes;
	unsigned		upper_var_buff_rows,	upper_pdf_buff_rows,	upper_cdf_buff_rows;
	unsigned		upper_var_buff_columns,	upper_pdf_buff_columns,	upper_cdf_buff_columns;
	dist_data		 h_lower_dist,  h_upper_dist;
	dist_data		dh_lower_dist, dh_upper_dist;
	dist_data		*d_lower_dist, *d_upper_dist;
	dist_container	*dh_dist_scatter;

	// compute some sizes
	unsigned total_rows = h_xsdata.energy_grid_len;
	unsigned total_cols = h_xsdata.total_reaction_channels + h_xsdata.n_isotopes;

	// allocate host arrays
	h_xsdata.dist_scatter	= new dist_container [total_rows * total_cols];
	dh_dist_scatter			= new dist_container [total_rows * total_cols];

	// init to null pointers
	for (row=0 ; row<h_xsdata.energy_grid_len; row++){  //start after the total xs vectors
			for (col=0 ; col<h_xsdata.total_reaction_channels ; col++){
				array_index = row*total_cols + col;
				 h_xsdata.dist_scatter[array_index].upper = 0x0;
				 h_xsdata.dist_scatter[array_index].lower = 0x0;
				       dh_dist_scatter[array_index].upper = 0x0;
				       dh_dist_scatter[array_index].lower = 0x0;
		}
	}

	// scan through the scattering array, copying the scattering distributions and replicating pointers to 
	// them for each main energy grid points within the distribution's grid points
	//
	//  --- The current method copies the scattering distributions twice.  
	//  --- This can be eliminated via copying pointers instead of reiniting ditributions, 
	//  --- but this has been left for now for simplicity.
	//
	for( col=h_xsdata.n_isotopes ; col<total_cols ; col++ ){ // going down column to stay within to same isotope, not effcient for caching, but OK here since only done at init
		
		// rest row index
		row = 0;

		// print some nice info if flagged
		if(print_flag >= 2){
			this_isotope = -1;
			for (int i=0;i<h_xsdata.n_isotopes+1;i++){
				if(col>=(h_xsdata.rxn_numbers_total[i]+h_xsdata.n_isotopes) & col<(h_xsdata.rxn_numbers_total[i+1]+h_xsdata.n_isotopes)){
					this_isotope = i;
					break;
				}
			}
			printf("   Getting scattering data for isotope %s - reaction %u\n",isotopes[this_isotope].c_str(),h_xsdata.rxn_numbers[col]);
		}

		while(row<total_rows){

			// compute array index
			array_index = row*total_cols + col;

			// call python object, which returns the two buffers 
			// and the main e grid index that it needs to be replicated to
			row_obj     = PyInt_FromLong     (row);
			col_obj     = PyInt_FromLong     (col);
			call_string = PyString_FromString("_get_scattering_data");
			obj_list    = PyObject_CallMethodObjArgs(xsdat_instance, call_string, row_obj, col_obj, NULL);
			PyErr_Print();

			// get objects in the returned list
			lower_erg_obj	= PyList_GetItem(obj_list,0);
			lower_len_obj	= PyList_GetItem(obj_list,1);
			lower_law_obj	= PyList_GetItem(obj_list,2);
			lower_intt_obj	= PyList_GetItem(obj_list,3);
			lower_var_obj	= PyList_GetItem(obj_list,4);
			lower_pdf_obj	= PyList_GetItem(obj_list,5);
			lower_cdf_obj	= PyList_GetItem(obj_list,6);
			upper_erg_obj	= PyList_GetItem(obj_list,7);
			upper_len_obj	= PyList_GetItem(obj_list,8);
			upper_law_obj	= PyList_GetItem(obj_list,9);
			upper_intt_obj	= PyList_GetItem(obj_list,10);
			upper_var_obj	= PyList_GetItem(obj_list,11);
			upper_pdf_obj	= PyList_GetItem(obj_list,12);
			upper_cdf_obj	= PyList_GetItem(obj_list,13);
			next_dex_obj	= PyList_GetItem(obj_list,14);
			PyErr_Print();

			// copy single values to temp objects
			h_lower_dist.erg	=	PyFloat_AsDouble(	lower_erg_obj);
			h_lower_dist.len	=	PyInt_AsLong(		lower_len_obj);
			h_lower_dist.law	=	PyInt_AsLong(		lower_law_obj);
			h_lower_dist.intt	=	PyInt_AsLong(		lower_intt_obj);
			h_upper_dist.erg	=	PyFloat_AsDouble(	upper_erg_obj);
			h_upper_dist.len	=	PyInt_AsLong(		upper_len_obj);
			h_upper_dist.law	=	PyInt_AsLong(		upper_law_obj);
			h_upper_dist.intt	=	PyInt_AsLong(		upper_intt_obj);
			next_dex			=	PyInt_AsLong(		next_dex_obj);
			PyErr_Print();


			// decide what to put in the array according to length reported
			if(h_lower_dist.len==0){   

				// below threshold, do nothing except go to where the threshold is
				row = next_dex;

			}
			else if (h_lower_dist.len==-1){  

				// this is a fission reaction and this is nu_t, nu_p
				// python routine linearly interpolates them, write current values
				nu_t = PyFloat_AsDouble(	lower_erg_obj);
				nu_p = PyFloat_AsDouble(	upper_erg_obj);
				memcpy(& h_xsdata.dist_scatter[array_index].upper , &nu_t, 1*sizeof(float) );
				memcpy(& h_xsdata.dist_scatter[array_index].lower , &nu_p, 1*sizeof(float) );
				memcpy(       &dh_dist_scatter[array_index].upper , &nu_t, 1*sizeof(float) );
				memcpy(       &dh_dist_scatter[array_index].lower , &nu_p, 1*sizeof(float) );

				// go to next dex so python can do interpolation
				row++;

			}
			else{
				// normal scattering distributions
				// get new pointers for temp arrays according to length reported
				h_lower_dist.var	=	new 	float[h_lower_dist.len];
				h_lower_dist.pdf	=	new 	float[h_lower_dist.len];
				h_lower_dist.cdf	=	new 	float[h_lower_dist.len];
				h_upper_dist.var	=	new 	float[h_upper_dist.len];
				h_upper_dist.pdf	=	new 	float[h_upper_dist.len];
				h_upper_dist.cdf	=	new 	float[h_upper_dist.len];

				// get buffers from python
				if (PyObject_CheckBuffer(lower_var_obj) &
					PyObject_CheckBuffer(lower_pdf_obj) &
					PyObject_CheckBuffer(lower_cdf_obj) &
					PyObject_CheckBuffer(upper_var_obj) &
					PyObject_CheckBuffer(upper_pdf_obj) &
					PyObject_CheckBuffer(upper_cdf_obj) )
					{
					PyObject_GetBuffer( lower_var_obj, &lower_var_buff, PyBUF_ND);
					PyObject_GetBuffer( lower_pdf_obj, &lower_pdf_buff, PyBUF_ND);
					PyObject_GetBuffer( lower_cdf_obj, &lower_cdf_buff, PyBUF_ND);
					PyObject_GetBuffer( upper_var_obj, &upper_var_buff, PyBUF_ND);
					PyObject_GetBuffer( upper_pdf_obj, &upper_pdf_buff, PyBUF_ND);
					PyObject_GetBuffer( upper_cdf_obj, &upper_cdf_buff, PyBUF_ND);
					PyErr_Print();
				}
				else{
					PyErr_Print();
				    fprintf(stderr, "Returned object does not support buffer interface\n");
				    return;
				}
	
				// get array size info 
				get_Py_buffer_dims(&lower_var_buff_rows, &lower_var_buff_columns, &lower_var_buff_bytes, &lower_var_buff);
				get_Py_buffer_dims(&lower_pdf_buff_rows, &lower_pdf_buff_columns, &lower_pdf_buff_bytes, &lower_pdf_buff);
				get_Py_buffer_dims(&lower_cdf_buff_rows, &lower_cdf_buff_columns, &lower_cdf_buff_bytes, &lower_cdf_buff);
				get_Py_buffer_dims(&upper_var_buff_rows, &upper_var_buff_columns, &upper_var_buff_bytes, &upper_var_buff);
				get_Py_buffer_dims(&upper_pdf_buff_rows, &upper_pdf_buff_columns, &upper_pdf_buff_bytes, &upper_pdf_buff);
				get_Py_buffer_dims(&upper_cdf_buff_rows, &upper_cdf_buff_columns, &upper_cdf_buff_bytes, &upper_cdf_buff);

				//make shapes and lengths are all the same
				assert(lower_pdf_buff_rows		==	lower_var_buff_rows		); // shape is in elements
				assert(lower_pdf_buff_columns	==	lower_var_buff_columns	);
				assert(lower_cdf_buff_rows		==	lower_var_buff_rows		);
				assert(lower_cdf_buff_columns	==	lower_var_buff_columns	);
				assert(upper_pdf_buff_rows		==	upper_var_buff_rows		);
				assert(upper_pdf_buff_columns	==	upper_var_buff_columns	);
				assert(upper_cdf_buff_rows		==	upper_var_buff_rows		);
				assert(upper_cdf_buff_columns	==	upper_var_buff_columns	);
				assert(lower_pdf_buff_bytes==lower_var_buff_bytes); // len is in BYTES
				assert(lower_cdf_buff_bytes==lower_var_buff_bytes);
				assert(upper_pdf_buff_bytes==upper_var_buff_bytes);
				assert(upper_cdf_buff_bytes==upper_var_buff_bytes);

				// make sure len corresponds to float32!
				//printf("buffer bytes %u elements %u float bytes %lu caluclated bytes %lu\n",lower_var_buff_bytes,(lower_var_buff_rows*lower_var_buff_columns),sizeof(float),(lower_var_buff_rows*lower_var_buff_columns)*sizeof(float));
				assert(lower_var_buff_bytes==(lower_var_buff_rows*lower_var_buff_columns)*sizeof(float));
				assert(lower_pdf_buff_bytes==(lower_pdf_buff_rows*lower_pdf_buff_columns)*sizeof(float));
				assert(lower_cdf_buff_bytes==(lower_cdf_buff_rows*lower_cdf_buff_columns)*sizeof(float));
				assert(upper_var_buff_bytes==(upper_var_buff_rows*upper_var_buff_columns)*sizeof(float));
				assert(upper_pdf_buff_bytes==(upper_pdf_buff_rows*upper_pdf_buff_columns)*sizeof(float));
				assert(upper_cdf_buff_bytes==(upper_cdf_buff_rows*upper_cdf_buff_columns)*sizeof(float));

				// allocate device distribution arrays
				cudaMalloc(&dh_lower_dist.var, lower_var_buff_bytes);   // len has been verified
				cudaMalloc(&dh_lower_dist.pdf, lower_var_buff_bytes);
				cudaMalloc(&dh_lower_dist.cdf, lower_var_buff_bytes);
				cudaMalloc(&dh_upper_dist.var, upper_var_buff_bytes);
				cudaMalloc(&dh_upper_dist.pdf, upper_var_buff_bytes);
				cudaMalloc(&dh_upper_dist.cdf, upper_var_buff_bytes);
				check_cuda(cudaPeekAtLastError());

				// copy data from python buffer to host pointer in array
				memcpy(h_lower_dist.var, lower_var_buff.buf, lower_var_buff_bytes);  
				memcpy(h_lower_dist.pdf, lower_pdf_buff.buf, lower_var_buff_bytes);
				memcpy(h_lower_dist.cdf, lower_cdf_buff.buf, lower_var_buff_bytes);
				memcpy(h_upper_dist.var, upper_var_buff.buf, upper_var_buff_bytes);
				memcpy(h_upper_dist.pdf, upper_pdf_buff.buf, upper_var_buff_bytes);
				memcpy(h_upper_dist.cdf, upper_cdf_buff.buf, upper_var_buff_bytes);

				// copy data from host arrays to device arrays
				cudaMemcpy(dh_lower_dist.var, h_lower_dist.var, lower_var_buff_bytes, cudaMemcpyHostToDevice);  
				cudaMemcpy(dh_lower_dist.pdf, h_lower_dist.pdf, lower_var_buff_bytes, cudaMemcpyHostToDevice);
				cudaMemcpy(dh_lower_dist.cdf, h_lower_dist.cdf, lower_var_buff_bytes, cudaMemcpyHostToDevice);
				cudaMemcpy(dh_upper_dist.var, h_upper_dist.var, upper_var_buff_bytes, cudaMemcpyHostToDevice);
				cudaMemcpy(dh_upper_dist.pdf, h_upper_dist.pdf, upper_var_buff_bytes, cudaMemcpyHostToDevice);
				cudaMemcpy(dh_upper_dist.cdf, h_upper_dist.cdf, upper_var_buff_bytes, cudaMemcpyHostToDevice);
				check_cuda(cudaPeekAtLastError());

				// copy vals in structure
				dh_lower_dist.erg	=	h_lower_dist.erg;
				dh_lower_dist.len	=	h_lower_dist.len;
				dh_lower_dist.law	=	h_lower_dist.law;
				dh_lower_dist.intt	=	h_lower_dist.intt;
				dh_upper_dist.erg	=	h_upper_dist.erg;
				dh_upper_dist.len	=	h_upper_dist.len;
				dh_upper_dist.law	=	h_upper_dist.law;
				dh_upper_dist.intt	=	h_upper_dist.intt;

				// allocate container
				cudaMalloc(&d_lower_dist, 1*sizeof(dist_data));
				cudaMalloc(&d_upper_dist, 1*sizeof(dist_data));
				check_cuda(cudaPeekAtLastError());

				// copy device pointers to device container
				cudaMemcpy(d_lower_dist, &dh_lower_dist, 1*sizeof(dist_data), cudaMemcpyHostToDevice);
				cudaMemcpy(d_upper_dist, &dh_upper_dist, 1*sizeof(dist_data), cudaMemcpyHostToDevice);
				check_cuda(cudaPeekAtLastError());

				// replicate pointers until next index
				for (int i=row ; i<next_dex; i++){
					array_index = i*total_cols + col;
					h_xsdata.dist_scatter[array_index].upper = &h_upper_dist;
					h_xsdata.dist_scatter[array_index].lower = &h_lower_dist;
					      dh_dist_scatter[array_index].upper =  d_upper_dist;
					      dh_dist_scatter[array_index].lower =  d_lower_dist;
				}

				// go to where the next index starts
				row = next_dex;

				//printf("next dex %u\n",next_dex);

			}

			PyErr_Print();
		}
	}

	// copy host array containing device pointers to device array
	cudaMalloc(&dh_xsdata.dist_scatter,                total_rows*total_cols*sizeof(dist_container));
	cudaMemcpy( dh_xsdata.dist_scatter,dh_dist_scatter,total_rows*total_cols*sizeof(dist_container),cudaMemcpyHostToDevice);
	check_cuda(cudaPeekAtLastError());

	// free host array containing device pointers, not needed anymore
	delete dh_dist_scatter;


//					// get flattened matrix for law 61
//					if (PyObject_CheckBuffer(mu_vector_obj)){
//						PyObject_GetBuffer(      mu_vector_obj,       &muBuff, PyBUF_ND);
//					}
//					else{
//						PyErr_Print();
//					    fprintf(stderr, "Returned object does not support buffer interface\n");
//					    return;
//					}
//	
//					//  get the buffers
//					muRows     =  muBuff.shape[0];
//					muColumns  =  muBuff.shape[1];
//					muBytes    =  muBuff.len;
//					//assert( muRows==1 | muColumns==1);  // make sure 1d
//	
//					this_pointer = new float [muBytes/sizeof(float)];
//					cudaMalloc(&cuda_pointer,muBytes);
//					total_bytes_scatter +=   muBytes;  // add to total count
//					xs_data_scatter     [k*MT_columns + j] = this_pointer;
//					xs_data_scatter_host[k*MT_columns + j] = cuda_pointer;
//					PyErr_Print();
//	
//					// copy to cuda pointer and host pointer
//					memcpy(this_pointer,muBuff.buf,muBytes);
//					cudaMemcpy(cuda_pointer,this_pointer,muBytes,cudaMemcpyHostToDevice);
//	
//				}	
//
//			}


}
void whistory::copy_energy_data(){
	// get scattering ditributions from PyNE and set up data heirarchy on host and device 

	if(print_flag >= 2){
		std::cout << "\e[1;32m" << "Loading/copying energy data and linking..." << "\e[m \n";
	}

	// python variables for arguments
	PyObject	*row_obj, *col_obj, *call_string, *obj_list;
	PyObject	*lower_erg_obj, *lower_len_obj, *lower_law_obj, *lower_intt_obj, *lower_var_obj, *lower_pdf_obj, *lower_cdf_obj; 
	PyObject	*upper_erg_obj, *upper_len_obj, *upper_law_obj, *upper_intt_obj, *upper_var_obj, *upper_pdf_obj, *upper_cdf_obj, *next_dex_obj; 
	Py_buffer	lower_var_buff, lower_pdf_buff, lower_cdf_buff, upper_var_buff, upper_pdf_buff, upper_cdf_buff;
	
	// local temp variables
	int				row, col;
	float			nu_t, nu_p;
	int				this_isotope;
	unsigned		array_index, next_dex, lower_rows;
	unsigned		lower_var_buff_bytes,	lower_pdf_buff_bytes,	lower_cdf_buff_bytes;
	unsigned		lower_var_buff_rows,	lower_pdf_buff_rows,	lower_cdf_buff_rows;
	unsigned		lower_var_buff_columns,	lower_pdf_buff_columns,	lower_cdf_buff_columns;
	unsigned		upper_var_buff_bytes,	upper_pdf_buff_bytes,	upper_cdf_buff_bytes;
	unsigned		upper_var_buff_rows,	upper_pdf_buff_rows,	upper_cdf_buff_rows;
	unsigned		upper_var_buff_columns,	upper_pdf_buff_columns,	upper_cdf_buff_columns;
	dist_data		 h_lower_dist,  h_upper_dist;
	dist_data		dh_lower_dist, dh_upper_dist;
	dist_data		*d_lower_dist, *d_upper_dist;
	dist_container	*dh_dist_energy;

	// compute some sizes
	unsigned total_rows = h_xsdata.energy_grid_len;
	unsigned total_cols = h_xsdata.total_reaction_channels + h_xsdata.n_isotopes;

	// allocate host arrays
	h_xsdata.dist_energy	= new dist_container [total_rows * total_cols];
	dh_dist_energy			= new dist_container [total_rows * total_cols];

	// init to null pointers
	for (row=0 ; row<h_xsdata.energy_grid_len; row++){  //start after the total xs vectors
			for (col=0 ; col<h_xsdata.total_reaction_channels ; col++){
				array_index = row*total_cols + col;
				 h_xsdata.dist_energy[array_index].upper = 0x0;
				 h_xsdata.dist_energy[array_index].lower = 0x0;
				       dh_dist_energy[array_index].upper = 0x0;
				       dh_dist_energy[array_index].lower = 0x0;
		}
	}

	// scan through the scattering array, copying the scattering distributions and replicating pointers to 
	// them for each main energy grid points within the distribution's grid points
	//
	//  --- The current method copies the scattering distributions twice.  
	//  --- This can be eliminated via copying pointers instead of reiniting ditributions, 
	//  --- but this has been left for now for simplicity.
	//
	for( col=h_xsdata.n_isotopes ; col<total_cols ; col++ ){ // going down column to stay within to same isotope, not effcient for caching, but OK here since only done at init
		
		// rest row index
		row = 0;

		// print some nice info if flagged
		if(print_flag >= 2){
			this_isotope = -1;
			for (int i=0;i<h_xsdata.n_isotopes+1;i++){
				if(col>=(h_xsdata.rxn_numbers_total[i]+h_xsdata.n_isotopes) & col<(h_xsdata.rxn_numbers_total[i+1]+h_xsdata.n_isotopes)){
					this_isotope = i;
					break;
				}
			}
			printf("   Getting energy data for isotope %s - reaction %u\n",isotopes[this_isotope].c_str(),h_xsdata.rxn_numbers[col]);
		}

		while(row<total_rows){

			// compute array index
			array_index = row*total_cols + col;

			// call python object, which returns the two buffers 
			// and the main e grid index that it needs to be replicated to
			row_obj     = PyInt_FromLong     (row);
			col_obj     = PyInt_FromLong     (col);
			call_string = PyString_FromString("_get_energy_data");
			obj_list    = PyObject_CallMethodObjArgs(xsdat_instance, call_string, row_obj, col_obj, NULL);
			PyErr_Print();

			// get objects in the returned list
			lower_erg_obj	= PyList_GetItem(obj_list,0);
			lower_len_obj	= PyList_GetItem(obj_list,1);
			lower_law_obj	= PyList_GetItem(obj_list,2);
			lower_intt_obj	= PyList_GetItem(obj_list,3);
			lower_var_obj	= PyList_GetItem(obj_list,4);
			lower_pdf_obj	= PyList_GetItem(obj_list,5);
			lower_cdf_obj	= PyList_GetItem(obj_list,6);
			upper_erg_obj	= PyList_GetItem(obj_list,7);
			upper_len_obj	= PyList_GetItem(obj_list,8);
			upper_law_obj	= PyList_GetItem(obj_list,9);
			upper_intt_obj	= PyList_GetItem(obj_list,10);
			upper_var_obj	= PyList_GetItem(obj_list,11);
			upper_pdf_obj	= PyList_GetItem(obj_list,12);
			upper_cdf_obj	= PyList_GetItem(obj_list,13);
			next_dex_obj	= PyList_GetItem(obj_list,14);
			PyErr_Print();

			// copy single values to temp objects
			h_lower_dist.erg	=	PyFloat_AsDouble(	lower_erg_obj);
			h_lower_dist.len	=	PyInt_AsLong(		lower_len_obj);
			h_lower_dist.law	=	PyInt_AsLong(		lower_law_obj);
			h_lower_dist.intt	=	PyInt_AsLong(		lower_intt_obj);
			h_upper_dist.erg	=	PyFloat_AsDouble(	upper_erg_obj);
			h_upper_dist.len	=	PyInt_AsLong(		upper_len_obj);
			h_upper_dist.law	=	PyInt_AsLong(		upper_law_obj);
			h_upper_dist.intt	=	PyInt_AsLong(		upper_intt_obj);
			next_dex			=	PyInt_AsLong(		next_dex_obj);
			PyErr_Print();


			// decide what to put in the array according to length reported
			if(h_lower_dist.len==0){   

				// below threshold, do nothing except go to where the threshold is
				row = next_dex;

			}
			else if (h_lower_dist.len==-1){  

				// this is a fission reaction and this is nu_t, nu_p
				// python routine linearly interpolates them, write current values
				nu_t = PyFloat_AsDouble(	lower_erg_obj);
				nu_p = PyFloat_AsDouble(	upper_erg_obj);
				memcpy(& h_xsdata.dist_energy[array_index].upper , &nu_t, 1*sizeof(float) );
				memcpy(& h_xsdata.dist_energy[array_index].lower , &nu_p, 1*sizeof(float) );
				memcpy(       &dh_dist_energy[array_index].upper , &nu_t, 1*sizeof(float) );
				memcpy(       &dh_dist_energy[array_index].lower , &nu_p, 1*sizeof(float) );

				// go to next dex so python can do interpolation
				row++;

			}
			else{
				// normal scattering distributions
				// get new pointers for temp arrays according to length reported
				h_lower_dist.var	=	new 	float[h_lower_dist.len];
				h_lower_dist.pdf	=	new 	float[h_lower_dist.len];
				h_lower_dist.cdf	=	new 	float[h_lower_dist.len];
				h_upper_dist.var	=	new 	float[h_upper_dist.len];
				h_upper_dist.pdf	=	new 	float[h_upper_dist.len];
				h_upper_dist.cdf	=	new 	float[h_upper_dist.len];

				// get buffers from python
				if (PyObject_CheckBuffer(lower_var_obj) &
					PyObject_CheckBuffer(lower_pdf_obj) &
					PyObject_CheckBuffer(lower_cdf_obj) &
					PyObject_CheckBuffer(upper_var_obj) &
					PyObject_CheckBuffer(upper_pdf_obj) &
					PyObject_CheckBuffer(upper_cdf_obj) )
					{
					PyObject_GetBuffer( lower_var_obj, &lower_var_buff, PyBUF_ND);
					PyObject_GetBuffer( lower_pdf_obj, &lower_pdf_buff, PyBUF_ND);
					PyObject_GetBuffer( lower_cdf_obj, &lower_cdf_buff, PyBUF_ND);
					PyObject_GetBuffer( upper_var_obj, &upper_var_buff, PyBUF_ND);
					PyObject_GetBuffer( upper_pdf_obj, &upper_pdf_buff, PyBUF_ND);
					PyObject_GetBuffer( upper_cdf_obj, &upper_cdf_buff, PyBUF_ND);
					PyErr_Print();
				}
				else{
					PyErr_Print();
				    fprintf(stderr, "Returned object does not support buffer interface\n");
				    return;
				}
	
				// get array size info 
				get_Py_buffer_dims(&lower_var_buff_rows, &lower_var_buff_columns, &lower_var_buff_bytes, &lower_var_buff);
				get_Py_buffer_dims(&lower_pdf_buff_rows, &lower_pdf_buff_columns, &lower_pdf_buff_bytes, &lower_pdf_buff);
				get_Py_buffer_dims(&lower_cdf_buff_rows, &lower_cdf_buff_columns, &lower_cdf_buff_bytes, &lower_cdf_buff);
				get_Py_buffer_dims(&upper_var_buff_rows, &upper_var_buff_columns, &upper_var_buff_bytes, &upper_var_buff);
				get_Py_buffer_dims(&upper_pdf_buff_rows, &upper_pdf_buff_columns, &upper_pdf_buff_bytes, &upper_pdf_buff);
				get_Py_buffer_dims(&upper_cdf_buff_rows, &upper_cdf_buff_columns, &upper_cdf_buff_bytes, &upper_cdf_buff);

				//make shapes and lengths are all the same
				assert(lower_pdf_buff_rows		==	lower_var_buff_rows		); // shape is in elements
				assert(lower_pdf_buff_columns	==	lower_var_buff_columns	);
				assert(lower_cdf_buff_rows		==	lower_var_buff_rows		);
				assert(lower_cdf_buff_columns	==	lower_var_buff_columns	);
				assert(upper_pdf_buff_rows		==	upper_var_buff_rows		);
				assert(upper_pdf_buff_columns	==	upper_var_buff_columns	);
				assert(upper_cdf_buff_rows		==	upper_var_buff_rows		);
				assert(upper_cdf_buff_columns	==	upper_var_buff_columns	);
				assert(lower_pdf_buff_bytes==lower_var_buff_bytes); // len is in BYTES
				assert(lower_cdf_buff_bytes==lower_var_buff_bytes);
				assert(upper_pdf_buff_bytes==upper_var_buff_bytes);
				assert(upper_cdf_buff_bytes==upper_var_buff_bytes);

				// make sure len corresponds to float32!
				//printf("buffer bytes %u elements %u float bytes %lu caluclated bytes %lu\n",lower_var_buff_bytes,(lower_var_buff_rows*lower_var_buff_columns),sizeof(float),(lower_var_buff_rows*lower_var_buff_columns)*sizeof(float));
				assert(lower_var_buff_bytes==(lower_var_buff_rows*lower_var_buff_columns)*sizeof(float));
				assert(lower_pdf_buff_bytes==(lower_pdf_buff_rows*lower_pdf_buff_columns)*sizeof(float));
				assert(lower_cdf_buff_bytes==(lower_cdf_buff_rows*lower_cdf_buff_columns)*sizeof(float));
				assert(upper_var_buff_bytes==(upper_var_buff_rows*upper_var_buff_columns)*sizeof(float));
				assert(upper_pdf_buff_bytes==(upper_pdf_buff_rows*upper_pdf_buff_columns)*sizeof(float));
				assert(upper_cdf_buff_bytes==(upper_cdf_buff_rows*upper_cdf_buff_columns)*sizeof(float));

				// allocate device distribution arrays
				cudaMalloc(&dh_lower_dist.var, lower_var_buff_bytes);   // len has been verified
				cudaMalloc(&dh_lower_dist.pdf, lower_var_buff_bytes);
				cudaMalloc(&dh_lower_dist.cdf, lower_var_buff_bytes);
				cudaMalloc(&dh_upper_dist.var, upper_var_buff_bytes);
				cudaMalloc(&dh_upper_dist.pdf, upper_var_buff_bytes);
				cudaMalloc(&dh_upper_dist.cdf, upper_var_buff_bytes);
				check_cuda(cudaPeekAtLastError());

				// copy data from python buffer to host pointer in array
				memcpy(h_lower_dist.var, lower_var_buff.buf, lower_var_buff_bytes);  
				memcpy(h_lower_dist.pdf, lower_pdf_buff.buf, lower_var_buff_bytes);
				memcpy(h_lower_dist.cdf, lower_cdf_buff.buf, lower_var_buff_bytes);
				memcpy(h_upper_dist.var, upper_var_buff.buf, upper_var_buff_bytes);
				memcpy(h_upper_dist.pdf, upper_pdf_buff.buf, upper_var_buff_bytes);
				memcpy(h_upper_dist.cdf, upper_cdf_buff.buf, upper_var_buff_bytes);

				// copy data from host arrays to device arrays
				cudaMemcpy(dh_lower_dist.var, h_lower_dist.var, lower_var_buff_bytes, cudaMemcpyHostToDevice);  
				cudaMemcpy(dh_lower_dist.pdf, h_lower_dist.pdf, lower_var_buff_bytes, cudaMemcpyHostToDevice);
				cudaMemcpy(dh_lower_dist.cdf, h_lower_dist.cdf, lower_var_buff_bytes, cudaMemcpyHostToDevice);
				cudaMemcpy(dh_upper_dist.var, h_upper_dist.var, upper_var_buff_bytes, cudaMemcpyHostToDevice);
				cudaMemcpy(dh_upper_dist.pdf, h_upper_dist.pdf, upper_var_buff_bytes, cudaMemcpyHostToDevice);
				cudaMemcpy(dh_upper_dist.cdf, h_upper_dist.cdf, upper_var_buff_bytes, cudaMemcpyHostToDevice);
				check_cuda(cudaPeekAtLastError());

				// copy vals in structure
				dh_lower_dist.erg	=	h_lower_dist.erg;
				dh_lower_dist.len	=	h_lower_dist.len;
				dh_lower_dist.law	=	h_lower_dist.law;
				dh_lower_dist.intt	=	h_lower_dist.intt;
				dh_upper_dist.erg	=	h_upper_dist.erg;
				dh_upper_dist.len	=	h_upper_dist.len;
				dh_upper_dist.law	=	h_upper_dist.law;
				dh_upper_dist.intt	=	h_upper_dist.intt;

				// allocate container
				cudaMalloc(&d_lower_dist, 1*sizeof(dist_data));
				cudaMalloc(&d_upper_dist, 1*sizeof(dist_data));
				check_cuda(cudaPeekAtLastError());

				// copy device pointers to device container
				cudaMemcpy(d_lower_dist, &dh_lower_dist, 1*sizeof(dist_data), cudaMemcpyHostToDevice);
				cudaMemcpy(d_upper_dist, &dh_upper_dist, 1*sizeof(dist_data), cudaMemcpyHostToDevice);
				check_cuda(cudaPeekAtLastError());

				// replicate pointers until next index
				for (int i=row ; i<next_dex; i++){
					array_index = i*total_cols + col;
					h_xsdata.dist_energy[array_index].upper = &h_upper_dist;
					h_xsdata.dist_energy[array_index].lower = &h_lower_dist;
					      dh_dist_energy[array_index].upper =  d_upper_dist;
					      dh_dist_energy[array_index].lower =  d_lower_dist;
				}

				// go to where the next index starts
				row = next_dex;

				//printf("next dex %u\n",next_dex);

			}

			PyErr_Print();
		}
	}

	// copy host array containing device pointers to device array
	cudaMalloc(&dh_xsdata.dist_energy,               total_rows*total_cols*sizeof(dist_container));
	cudaMemcpy( dh_xsdata.dist_energy,dh_dist_energy,total_rows*total_cols*sizeof(dist_container),cudaMemcpyHostToDevice);
	check_cuda(cudaPeekAtLastError());

	// free host array containing device pointers, not needed anymore
	delete dh_dist_energy;

}
void whistory::init_cross_sections(){
	
	if(print_flag >= 2){
		std::cout << "\e[1;32m" << "Loading cross sections and unionizing..." << "\e[m \n";
	}

	// Python variables
	int i, do_final;
	unsigned* length_numbers = new unsigned [3];

	// Make python object ready for queries - init python, load cross sections, etc.
	do_final = init_python();

	// copy various cross section arrays from python
	copy_python_buffer(&dh_xsdata.xs,         			&h_xsdata.xs,         			"_get_MT_array_pointer");
	copy_python_buffer(&dh_xsdata.energy_grid,			&h_xsdata.energy_grid,			"_get_main_Egrid_pointer");
	copy_python_buffer(&dh_xsdata.rxn_numbers,			&h_xsdata.rxn_numbers,			"_get_MT_numbers_pointer");
	copy_python_buffer(&dh_xsdata.rxn_numbers_total,	&h_xsdata.rxn_numbers_total,	"_get_MT_numbers_total_pointer");
	copy_python_buffer(&dh_xsdata.awr,					&h_xsdata.awr,					"_get_awr_pointer");
	copy_python_buffer(&dh_xsdata.temp,					&h_xsdata.temp,					"_get_temp_pointer");
	copy_python_buffer(&dh_xsdata.Q,					&h_xsdata.Q,					"_get_Q_pointer");
	copy_python_buffer(									&length_numbers,				"_get_length_numbers_pointer");
	dh_xsdata.n_isotopes				= h_xsdata.n_isotopes				= length_numbers[0];				
	dh_xsdata.energy_grid_len			= h_xsdata.energy_grid_len			= length_numbers[1];	
	dh_xsdata.total_reaction_channels	= h_xsdata.total_reaction_channels	= length_numbers[2];

	// copy scattering data
	copy_scatter_data();
	check_cuda(cudaPeekAtLastError());

	// copy energy data
	copy_energy_data();
	check_cuda(cudaPeekAtLastError());

	// intialization complete, copy host structure (containing device pointers) to device structure
	cudaMemcpy( d_xsdata,	&dh_xsdata,	1*sizeof(cross_section_data),	cudaMemcpyHostToDevice);
	check_cuda(cudaPeekAtLastError());

	// launch a test kernel...
	test_function( NUM_THREADS, 1, d_xsdata, d_particles, d_tally);
	check_cuda(cudaPeekAtLastError());

	// finalize python if initialized by warp
	if(do_final){
		Py_Finalize();
	}

	if(print_flag >= 2){//

	}

	std::cout << "Done." << "\e[m \n";

	std::cout << "\e[1;32m" << "Making material table..." << "\e[m \n";

	//pass awr pointer to geometry object, make the number density table, copy pointers back
	problem_geom.awr_list = h_xsdata.awr;
	problem_geom.make_material_table();
	problem_geom.get_material_table(&n_materials,&n_isotopes,&number_density_matrix);  
	assert(n_isotopes == h_xsdata.n_isotopes);

	// allocate material matrix (now that we have n_materials) and copy to device
	cudaMalloc(&d_number_density_matrix,                        n_materials*n_isotopes*sizeof(float));
	cudaMemcpy( d_number_density_matrix, number_density_matrix, n_materials*n_isotopes*sizeof(float), cudaMemcpyHostToDevice);
	check_cuda(cudaPeekAtLastError());
	
	std::cout << "Done." << "\e[m \n";

}
void whistory::print_xs_data(){  // 0=isotopes, 1=main E points, 2=total numer of reaction channels, 3=matrix E points, 4=angular cosine points, 5=outgoing energy points
//	unsigned dsum = 0;
//	printf("\e[1;32m%-6s\e[m \n","Cross section data info:");
//	std::cout << "--- Bytes ---" << "\n";
//	std::cout << "  xs_length_numbers:        " << 6*sizeof(unsigned) << "\n";			  dsum += (6*sizeof(unsigned) );
//	std::cout << "  xs_MT_numbers_total:      " << xs_length_numbers[0]*sizeof(unsigned) 	<< "\n";  dsum += (xs_length_numbers[0]*sizeof(unsigned) );
//	std::cout << "  xs_MT_numbers:            " << (xs_length_numbers[2]+xs_length_numbers[0])*sizeof(unsigned)<< "\n";
//	dsum += (xs_length_numbers[2]*sizeof(unsigned) );
//	std::cout << "  xs_data_main_E_grid:      " << xs_length_numbers[1]	*sizeof(float)	<< "\n";  dsum += (xs_length_numbers[1]*sizeof(float)	  );
//	std::cout << "  xs_data_MT:               " << MT_rows*MT_columns*sizeof(float)		<< "\n";  dsum += (MT_rows*MT_columns)*sizeof(float);
//	std::cout << "  xs_data_scatter_pointers: " << MT_rows*MT_columns*sizeof(float*)		<< "\n";  dsum += (MT_rows*MT_columns)*sizeof(float*);
//	std::cout << "  xs_data_energy_pointers:  " << MT_rows*MT_columns*sizeof(float*)		<< "\n";  dsum += (MT_rows*MT_columns)*sizeof(float*);
//	std::cout << "  scatter data:             " << total_bytes_scatter				<< "\n";  dsum += (total_bytes_scatter);
//	std::cout << "  energy data:              " << total_bytes_energy				<< "\n";  dsum += (total_bytes_energy);
//	std::cout << "  TOTAL:                    " << dsum << " bytes \n";
//	std::cout << "  TOTAL:                    " << dsum/1048576 << " MB \n";
}
void whistory::write_xs_data(std::string filename){

//	std::cout << "\e[1;32m" << "Writing xs_data to " << filename << "... ";
//
//	std::string this_name;
//	// write MT array
//	this_name = filename + ".MTarray";
//	FILE* xsfile = fopen(this_name.c_str(),"w");
//	this_name = filename + ".scatterptr";
//	FILE* scatterfile = fopen(this_name.c_str(),"w");
//	this_name = filename + ".energyptr";
//	FILE* energyfile = fopen(this_name.c_str(),"w");
//	for (int k=0;k<MT_rows;k++){
//		for(int j=0;j<MT_columns;j++){
//			fprintf(xsfile,"% 10.8E ",xs_data_MT[k*MT_columns+j]);
//			fprintf(scatterfile,"%p ",xs_data_scatter_host[k*MT_columns+j]);
//			fprintf(energyfile,"%p ",xs_data_energy_host[k*MT_columns+j]);
//		}
//		fprintf(xsfile,"\n");
//		fprintf(scatterfile,"\n");
//		fprintf(energyfile,"\n");
//	}
//	fclose(xsfile);
//	fclose(scatterfile);
//	fclose(energyfile);
//
//	// write unionized E grid
//	this_name = filename + ".Egrid";
//	xsfile = fopen(this_name.c_str(),"w");
//	for (int k=0;k<MT_rows;k++){
//		fprintf(xsfile,"%10.8E\n",xs_data_main_E_grid[k]);
//	}
//	fclose(xsfile);
//
//	// write MT number array
//	this_name = filename + ".MTnums";
//	xsfile = fopen(this_name.c_str(),"w");
//	for(int j=0;j<MT_columns;j++){
//		fprintf(xsfile,"%u\n",xs_MT_numbers[j]);
//	}
//	fclose(xsfile);
//
//	// write (hopefully) covnerged fission source
//	cudaMemcpy(space, d_fissile_points, Ndataset*sizeof(spatial_data), cudaMemcpyDeviceToHost);
//	cudaMemcpy(E,     d_fissile_energy, Ndataset*sizeof(float),        cudaMemcpyDeviceToHost);
//	this_name = filename + ".fission_source";
//	xsfile = fopen(this_name.c_str(),"w");
//	for(int j=0;j<Ndataset;j++){
//		fprintf(xsfile,"% 6.4E % 6.4E % 6.4E %6.4E\n",space[j].x,space[j].y,space[j].z,E[j]);
//	}
//	fclose(xsfile); 


}
void whistory::print_pointers(){
//	std::cout << "\e[1;32m" << "Pointer Info:" << "\e[m \n";
//	std::cout << "--- HOST ---" << "\n";
//	std::cout << "  space:               " <<   space   << "\n";
//	std::cout << "  E:                   " <<   E       << "\n";
//	std::cout << "  Q:                   " <<   Q       << "\n";
//	std::cout << "  rn_bank:             " <<   rn_bank << "\n";
//	std::cout << "  cellnum:             " <<   cellnum << "\n";
//	std::cout << "  matnum:              " <<   matnum  << "\n";
//	std::cout << "  isonum:              " <<   isonum  << "\n";
//	std::cout << "  rxn:                 " <<   rxn     << "\n";
//	std::cout << "  done:                " <<   done    << "\n";
//	std::cout << "  yield:               " <<   yield   << "\n";
//	std::cout << "  xs_length_numbers:   " << xs_length_numbers   << "\n"; 
//	std::cout << "  xs_MT_numbers_total: " << xs_MT_numbers_total << "\n";
//	std::cout << "  xs_MT_numbers:       " << xs_MT_numbers       << "\n";
//	std::cout << "  xs_data_MT:          " << xs_data_MT          << "\n";
//	std::cout << "  xs_data_main_E_grid: " << xs_data_main_E_grid << "\n";
//	std::cout << "--- DEVICE ---" << "\n";
//	std::cout << "d_space:               " << d_space   << "\n";
//	std::cout << "d_E:                   " << d_E       << "\n";
//	std::cout << "d_Q:                   " << d_Q       << "\n";
//	std::cout << "h_particles.rn_bank:             " << h_particles.rn_bank << "\n";
//	std::cout << "d_cellnum:             " << d_cellnum << "\n";
//	std::cout << "d_matnum:              " << d_matnum  << "\n";
//	std::cout << "d_isonum:              " << d_isonum  << "\n";
//	std::cout << "d_rxn:                 " << d_rxn     << "\n";
//	std::cout << "d_done:                " << d_done    << "\n";
//	std::cout << "d_yield:               " << d_yield   << "\n";
//	std::cout << "d_xs_length_numbers:   " << d_xs_length_numbers   << "\n"; 
//	std::cout << "d_xs_MT_numbers_total: " << d_xs_MT_numbers_total << "\n";
//	std::cout << "d_xs_MT_numbers:       " << d_xs_MT_numbers       << "\n";
//	std::cout << "d_xs_data_MT:          " << d_xs_data_MT          << "\n";
//	std::cout << "d_xs_data_main_E_grid: " << d_xs_data_main_E_grid << "\n";
}
void whistory::trace(unsigned type){

	//make sure device is ready
	cudaDeviceSynchronize();

	//trace
	optix_obj.trace(type);

}
void whistory::trace(unsigned type, unsigned n_active){

	//make sure device is ready
	cudaDeviceSynchronize();

	//trace
	optix_obj.trace(type,n_active);

}
void whistory::print_materials_table(){

	// print materials
	problem_geom.print_materials_table();

	//print temperatures:
	printf("  --------------   \n");
	for (int i=0; i<n_isotopes;i++){
		printf("  Isotope %3d = %10s , temp: %7.2F K\n",i,isotopes[i].c_str(),h_xsdata.temp[i]/8.617332478e-11);
	}

}
void whistory::sample_fissile_points(){

	if(print_flag >= 2){
		std::cout << "\e[1;32m" << "Sampling initial fissile starting points uniformly... " << "\e[m \n";
	}

	// iterate
	unsigned current_index = 0;
	unsigned valid_N = 0;

	while (current_index < N){
		
		// advance RN bank
		update_RNG();

		// set uniformly random positions on GPU
		set_positions_rand ( NUM_THREADS, N , outer_cell_type, dh_particles.space , dh_particles.rn_bank, outer_cell_dims);
		
		//run OptiX to get cell number, set as a hash run for fissile, writes 1/0 to matnum
		trace(3, N);

		// compact
		res = cudppCompact(compactplan, d_valid_result, (size_t*)d_valid_N , d_remap , dh_particles.matnum , N);
		if (res != CUDPP_SUCCESS){fprintf(stderr, "Error in compacting\n");exit(-1);}

		// copy back and add
		cudaMemcpy( &valid_N, d_valid_N, 1*sizeof(unsigned), cudaMemcpyDeviceToHost);
		current_index += valid_N;

		//copy in new values, keep track of index, copies positions and direction
		copy_points( NUM_THREADS, valid_N, d_valid_N, current_index, d_valid_result, d_fissile_points, dh_particles.space, d_fissile_energy, dh_particles.E); 
		
		// print how far along we are
		std::cout << (float)current_index/(float)N*100.0 <<" \% done\r";

		if((float)current_index/(float)N > 1){
			if(print_flag >= 2){
				std::cout << "  100.00 \% done     \n";
			}
		} 
	}

	float a=0.7;
	float b=1.0;

	if(print_flag >= 2){
		printf("  Sampling Watt fission spectrum with a=%6.4f and b=%6.4f...\n",a,b);
	}

	sample_fissile_energy(NUM_THREADS, N, a, b, dh_particles.rn_bank, d_fissile_energy);


	if(print_flag >= 2){
		std::cout << "  Copying to starting points...\n";
	}

	cudaMemcpy(dh_particles.space,	d_fissile_points	,	N*sizeof(spatial_data),	cudaMemcpyDeviceToDevice);
	cudaMemcpy(dh_particles.E,		d_fissile_energy	,	N*sizeof(float),		cudaMemcpyDeviceToDevice);
	cudaMemcpy(dh_particles.rxn,	zeros 				,	N*sizeof(unsigned),		cudaMemcpyHostToDevice);
	//cudaFree(d_fissile_points);

	if(print_flag >= 2){
		std::cout << "  Done.\n";
	}

	//write starting positions to file
	if(dump_flag >= 3){
		cudaMemcpy(h_particles.space,dh_particles.space,N*sizeof(spatial_data),cudaMemcpyDeviceToHost);
		FILE* positionsfile = fopen("starting_positions","w");
		for(int k=0;k<N;k++){
			fprintf(positionsfile,"% 10.8E % 10.8E % 10.8E % 10.8E % 10.8E % 10.8E\n",h_particles.space[k].x,h_particles.space[k].y,h_particles.space[k].z,h_particles.space[k].xhat,h_particles.space[k].yhat,h_particles.space[k].zhat);
		}
		fclose(positionsfile);
	}

	// advance RN bank
	update_RNG();

}
void whistory::write_to_file(spatial_data* array_in , unsigned N , std::string filename, std::string opentype){

	FILE* f  = fopen(filename.c_str(),opentype.c_str());
	spatial_data * hostdata = new spatial_data [N];
	cudaMemcpy(hostdata,array_in,N*sizeof(spatial_data),cudaMemcpyDeviceToHost);

	for(unsigned k = 0;  k<N ;k++){
		fprintf(f,"% 6.4E % 6.4E % 6.4E\n",hostdata[k].x,hostdata[k].y,hostdata[k].z);
	}

	delete hostdata;
	fclose(f);

}
void whistory::write_to_file(spatial_data* array_in , float* array_in2, unsigned N , std::string filename, std::string opentype){

	FILE* f  = fopen(filename.c_str(),opentype.c_str());
	spatial_data * hostdata = new spatial_data [N];
	float * hostdata2 = new float [N];
	cudaMemcpy(hostdata,array_in,N*sizeof(spatial_data),cudaMemcpyDeviceToHost);
	cudaMemcpy(hostdata2,array_in2,N*sizeof(float),cudaMemcpyDeviceToHost);

	for(unsigned k = 0;  k<N ;k++){
		fprintf(f,"% 6.4E % 6.4E % 6.4E % 6.4E\n",hostdata[k].x,hostdata[k].y,hostdata[k].z,hostdata2[k]);
	}

	delete hostdata;
	delete hostdata2;
	fclose(f);

}
void whistory::write_to_file(unsigned* array_in , unsigned N , std::string filename, std::string opentype){

	FILE* f  = fopen(filename.c_str(),opentype.c_str());
	unsigned * hostdata = new unsigned [N];
	cudaMemcpy(hostdata,array_in,N*sizeof(unsigned),cudaMemcpyDeviceToHost);

	for(unsigned k = 0;  k<N ;k++){
		fprintf(f,"%u\n",hostdata[k]);
	}

	delete hostdata;
	fclose(f);

}
void whistory::write_to_file(unsigned* array_in , unsigned* array_in2, unsigned N , std::string filename, std::string opentype){

	FILE* f  = fopen(filename.c_str(),opentype.c_str());
	unsigned * hostdata = new unsigned [N];
	unsigned * hostdata2 = new unsigned [N];
	cudaMemcpy(hostdata,array_in,N*sizeof(unsigned),cudaMemcpyDeviceToHost);
	cudaMemcpy(hostdata2,array_in2,N*sizeof(unsigned),cudaMemcpyDeviceToHost);

	for(unsigned k = 0;  k<N ;k++){
		fprintf(f,"%u %u\n",hostdata[k],hostdata2[k]);
	}

	delete hostdata;
	delete hostdata2;
	fclose(f);

}
void whistory::reset_cycle(float keff_cycle){
//
//	// re-base the yield so keff is 1
//	rebase_yield( NUM_THREADS, N,  keff_cycle, dh_particles.rn_bank, dh_particles.yield);
//	check_cuda(cudaPeekAtLastError());
//
//	// prefix sum scan (non-inclusive) yields to see where to write
//	res = cudppScan( scanplan_int, d_scanned,  dh_particles.yield,  Ndataset );
//	if (res != CUDPP_SUCCESS){fprintf(stderr, "Error in scanning yield values\n");exit(-1);}
//	check_cuda(cudaPeekAtLastError());
//
//	// swap the key/values to sort the rxn vector by tid to align the rest of data with the right reactions, done so 18 doesn't have to be assumed.  Sorting is necessary for prefix sum to work (order-dependent)
//	res = cudppRadixSort(radixplan, d_remap, dh_particles.rxn, N);  
//	if (res != CUDPP_SUCCESS){fprintf(stderr, "Error in sorting reactions\n");exit(-1);}
//	check_cuda(cudaPeekAtLastError());
//
//	//pop them in!  should be the right size now due to keff rebasing  
//	pop_source( NUM_THREADS, N,  d_remap, d_isonum, d_scanned, d_yield, d_index, d_rxn, d_space, d_E , dh_particles.rn_bank , d_xs_data_energy, d_xs_data_scatter, d_fissile_points, d_fissile_energy, d_awr_list, d_weight);
//	check_cuda(cudaPeekAtLastError());
//
// 	// rest run arrays
//	cudaMemcpy( d_space,	d_fissile_points,	N*sizeof(spatial_data),	cudaMemcpyDeviceToDevice );
//	cudaMemcpy( d_E,		d_fissile_energy,	N*sizeof(unsigned),		cudaMemcpyDeviceToDevice );
//	cudaMemcpy( d_done,		zeros,				N*sizeof(unsigned),		cudaMemcpyHostToDevice );
//	cudaMemcpy( d_cellnum,	zeros,				N*sizeof(unsigned),		cudaMemcpyHostToDevice );
//	cudaMemcpy( d_matnum,	zeros,				N*sizeof(unsigned),		cudaMemcpyHostToDevice );
//	cudaMemcpy( d_isonum,	zeros,				N*sizeof(unsigned),		cudaMemcpyHostToDevice );
//	cudaMemcpy( d_yield,	zeros,				N*sizeof(unsigned),		cudaMemcpyHostToDevice );
//	cudaMemcpy( d_rxn,		zeros,				N*sizeof(unsigned),		cudaMemcpyHostToDevice );
//	cudaMemcpy( d_remap,	remap,				N*sizeof(unsigned),		cudaMemcpyHostToDevice );
//	cudaMemcpy( d_index,	zeros,				N*sizeof(unsigned),		cudaMemcpyHostToDevice );
//	cudaMemcpy( d_weight,	fones,				N*sizeof(float),		cudaMemcpyHostToDevice );
//	check_cuda(cudaPeekAtLastError());
//
//	// update RNG seeds
//	update_RNG();
//	check_cuda(cudaPeekAtLastError());
//
//	// sync, these H2D and D2D copies aren't strictly synchronous
//	cudaDeviceSynchronize();
//	check_cuda(cudaPeekAtLastError());

}
void whistory::reset_fixed(){

//	// rest read-in run arrays (ie not ones that are written to in-between)
//	cudaMemcpy( d_space,		space,		Ndataset*sizeof(spatial_data),	cudaMemcpyHostToDevice );
//	cudaMemcpy( d_done,			done,		Ndataset*sizeof(unsigned),		cudaMemcpyHostToDevice );
//	cudaMemcpy( d_cellnum,		cellnum,	Ndataset*sizeof(unsigned),		cudaMemcpyHostToDevice );
//	cudaMemcpy( d_matnum,		matnum,		Ndataset*sizeof(unsigned),		cudaMemcpyHostToDevice );
//	cudaMemcpy( d_isonum,		isonum,		Ndataset*sizeof(unsigned),		cudaMemcpyHostToDevice );
//	cudaMemcpy( d_yield,		zeros,		Ndataset*sizeof(unsigned),		cudaMemcpyHostToDevice );
//	cudaMemcpy( d_rxn,			zeros,		Ndataset*sizeof(unsigned),		cudaMemcpyHostToDevice );
//	cudaMemcpy( d_active,		remap,		Ndataset*sizeof(unsigned),		cudaMemcpyHostToDevice );
//
//	//set position, direction, energy
//	sample_fixed_source( NUM_THREADS, N, d_active, dh_particles.rn_bank, d_E, d_space);
//
//	// update RNG seeds
//	update_RNG();
//
//	// sync, these H2D and D2D copies aren't strictly synchronous
//	cudaDeviceSynchronize();

}
void whistory::run(){

	std::string runtype = "UNDEFINED";
	if     (RUN_FLAG==0){runtype="fixed";}
	else if(RUN_FLAG==1){runtype="criticality";}

	double keff = 0.0;
	float keff_cycle = 0.0;
	float it = 0.0;
	unsigned current_fission_index = 0;
	unsigned current_fission_index_temp = 0;
	unsigned Nrun = N;
	int difference = 0;
	int iteration = 0;
	int iteration_total=0;
	unsigned converged = 0;
	unsigned active_size1, active_gap, escatter_N, escatter_start, iscatter_N, iscatter_start, cscatter_N, cscatter_start, fission_N, fission_start;
	std::string fiss_name, stats_name;
	FILE* statsfile;
	float runtime = get_time();

	if(RUN_FLAG==0){
		reset_fixed();
		converged=1;
		n_skip=0;
	}
	else if(RUN_FLAG==1){
		sample_fissile_points();
	}

	// init run vars for cycle
	if(RUN_FLAG==0){
		Nrun=Ndataset;
	}
	else if(RUN_FLAG==1){
		Nrun=N;
	}
	keff_cycle = 0;

	if(print_flag>=1){
		printf("\e[1;32m--- Running in\e[m \e[1;31m%s\e[m \e[1;32msource mode --- \e[m \n",runtype.c_str());
		printf("\e[1;32m--- Skipping %u cycles, Running %u ACTIVE CYCLES, %u histories each--- \e[m \n",n_skip,n_cycles,N);
	}

	// make sure fissile_points file is cleared
	//if(dump_flag>=3){          // level 3 includes fissile points
	//	fiss_name=filename;
	//	fiss_name.append(".fission_points.png");
	//	FILE* ffile = fopen(fiss_name.c_str(),"w");
	//	fclose(ffile);
	//}

	// open file for run stats
	if(dump_flag>=2){         // level 2 includes stats files
		stats_name=filename;
		stats_name.append(".stats");
		statsfile = fopen(stats_name.c_str(),"w");
	}

//	while(iteration<n_cycles){
//
//		//write source positions to file if converged
//		if(converged){
//			if(dump_flag>=3){
//				//write_to_file(d_space,d_E,N,fiss_name,"a+");
//				bin_fission_points(dh_particles.space,N);
//			}
//		}
//
//		// reset edges and active number
//		Nrun=N;
//		edges[0]  = 0; 
//		edges[1]  = 0; 
//		edges[2]  = 0;
//		edges[3]  = 0;
//		edges[4]  = 0;
//		edges[5]  = 0;
//		edges[6]  = 0; 
//		edges[7]  = 0;
//		edges[8]  = 0;
//		edges[9]  = 0;
//		edges[10] = 0;
//
//		// check any previous errors
//		check_cuda(cudaPeekAtLastError());
//
//		// process batch
//		while(Nrun>0){
//
//			//record stats
//			if(dump_flag >= 2){
//				fprintf(statsfile,"%u %10.8E\n",Nrun,get_time());
//			}
//			if(print_flag>=3){printf("TOP OF CYCLE, Nrun = %u\n",Nrun);
//			}
//	
//			// find what material we are in and nearest surface distance
//			trace(2, Nrun);
//			check_cuda(cudaPeekAtLastError());
//
//			//find the main E grid index
//			find_E_grid_index( NUM_THREADS, Nrun, xs_length_numbers[1],  d_remap, d_xs_data_main_E_grid, d_E, d_index, d_rxn);
//			check_cuda(cudaPeekAtLastError());
//
//			// run macroscopic kernel to find interaction length, macro_t, and reaction isotope, move to interactino length, set resample flag, 
//			macroscopic( NUM_THREADS, Nrun,  n_isotopes, n_materials, MT_columns, outer_cell, d_remap, d_space, d_isonum, d_cellnum, d_index, d_matnum, d_rxn, d_xs_data_main_E_grid, dh_particles.rn_bank, d_E, d_xs_data_MT , d_number_density_matrix, d_done);
//			check_cuda(cudaPeekAtLastError());
//
//			// run tally kernel to compute spectra
//			if(converged){
//				tally_spec( NUM_THREADS, Nrun, n_tally, tally_cell, d_remap, d_space, d_E, d_tally_score, d_tally_square, d_tally_count, d_done, d_cellnum, d_rxn, d_weight);
//			}
//			check_cuda(cudaPeekAtLastError());
//
//			// run microscopic kernel to find reaction type
//			microscopic( NUM_THREADS, Nrun, n_isotopes, MT_columns, d_remap, d_isonum, d_index, d_xs_data_main_E_grid, dh_particles.rn_bank, d_E, d_xs_data_MT , d_xs_MT_numbers_total, d_xs_MT_numbers, d_xs_data_Q, d_rxn, d_Q, d_done);
//			check_cuda(cudaPeekAtLastError());
//
//			// remap threads to still active data
//			remap_active(&Nrun, &escatter_N, &escatter_start, &iscatter_N, &iscatter_start, &cscatter_N, &cscatter_start, &fission_N, &fission_start);
//			check_cuda(cudaPeekAtLastError());
//
//			// concurrent calls to do escatter/iscatter/cscatter/fission
//			cudaThreadSynchronize();
//			cudaDeviceSynchronize();
//			escatter( stream[0], NUM_THREADS,   escatter_N, escatter_start, d_remap, d_isonum, d_index, dh_particles.rn_bank, d_E, d_space, d_rxn, d_awr_list, d_temp_list, d_done, d_xs_data_scatter);
//			iscatter( stream[1], NUM_THREADS,   iscatter_N, iscatter_start, d_remap, d_isonum, d_index, dh_particles.rn_bank, d_E, d_space, d_rxn, d_awr_list, d_Q, d_done, d_xs_data_scatter, d_xs_data_energy);
//			cscatter( stream[2], NUM_THREADS,1, cscatter_N, cscatter_start, d_remap, d_isonum, d_index, dh_particles.rn_bank, d_E, d_space, d_rxn, d_awr_list, d_Q, d_done, d_xs_data_scatter, d_xs_data_energy); // 1 is for transport run mode, as opposed to 'pop' mode
//			fission ( stream[3], NUM_THREADS,   fission_N,  fission_start , d_remap, d_isonum, d_index, dh_particles.rn_bank, d_E, d_space, d_rxn, d_awr_list, d_yield, d_weight, d_xs_data_scatter, d_xs_data_energy);  
//			cudaDeviceSynchronize();
//			check_cuda(cudaPeekAtLastError());
//
//			if(RUN_FLAG==0){  //fixed source
//				// pop secondaries back in
//				keff_cycle += reduce_yield();
//				prep_secondaries();
//				pop_secondaries( NUM_THREADS, Ndataset, d_completed, d_scanned, d_yield, d_done, d_index, d_rxn, d_space, d_E , dh_particles.rn_bank , d_xs_data_energy);
//				//if(reduce_yield()!=0.0){printf("pop_secondaries did not reset all yields!\n");}
//			}
//			check_cuda(cudaPeekAtLastError());
//
//			//std::cout << "press enter to continue...\n";
//			//std::cin.ignore();
//
//		}
//
//		//reduce yield and reset cycle
//		if(RUN_FLAG==0){
//			keff_cycle = 1.0 - 1.0/(keff_cycle+1.0);   // based on: Ntotal = Nsource / (1-k) 
//			reset_fixed();
//			Nrun=Ndataset;
//		}
//		else if (RUN_FLAG==1){	
//			accumulate_keff(converged, iteration, &keff, &keff_cycle);
//			check_cuda(cudaPeekAtLastError());
//			accumulate_tally();
//			check_cuda(cudaPeekAtLastError());
//			reset_cycle(keff_cycle);
//			check_cuda(cudaPeekAtLastError());
//			Nrun=N;
//		}
//
//		// update active iteration
//		if (converged){
//			iteration++;
//		}
//
//		// print whatever's clever
//		if(print_flag >= 1){
//			if(converged){
//				     if(RUN_FLAG==0){std::cout << "Cumulative keff/sc-mult = " << keff << " / " << 1.0/(1.0-keff) << ", ACTIVE cycle " << iteration << ", cycle keff/sc-mult = " << keff_cycle << " / " << 1.0/(1.0-keff_cycle) << "\n";}
//				else if(RUN_FLAG==1){printf("Cumulative keff =  %8.6E +- %6.4E , ACTIVE cycle %4u, cycle keff = %8.6E\n",keff,keff_err,iteration,keff_cycle);}
//			}
//			else{
//				printf("Converging fission source... skipped cycle %4u\n",iteration_total+1);
//			}
//		}
//
//		if(dump_flag>=2){
//			fprintf(statsfile,"---- iteration %u done ----\n",iteration);
//		}
//
//		//std::cout << "press enter to continue...\n";
//		//std::cin.ignore();
//		//exit(0);
//
//		// set convergence flag
//		if( iteration_total == n_skip-1){ 
//			converged=1;
//			reduced_yields_total = 0;
//		}
//
//		// advance iteration number, reset cycle keff
//		iteration_total++;
//		keff_cycle = 0.0;
//
//	}
//
//
//	// print total transport running time
//	runtime = get_time() - runtime;
//	if(print_flag >= 1){
//		if(runtime>60.0){
//			std::cout << "RUNTIME = " << runtime/60.0 << " minutes.\n";
//		}
//		else{
//			std::cout << "RUNTIME = " << runtime << " seconds.\n";
//		}
//	}
//
//	if(dump_flag>=1){
//		write_results(runtime,keff,"w");
//	}
//	if(dump_flag>=2){
//		fclose(statsfile);
//	}
//	if(dump_flag>=3){
//		write_fission_points();
//	}

}
void whistory::write_results(float runtime, float keff, std::string opentype){

	std::string resname = filename;
	resname.append(".results");
	FILE* f  = fopen(resname.c_str(),opentype.c_str());
	fprintf(f,"runtime = % 10.8E m    k_eff = % 10.8E   rel. err. = % 10.8E     \ncycles total %u skipped %u scored %u \n %u source neutrons per cycle\n",runtime/60.0,keff,keff_err,n_skip+n_cycles,n_skip,n_cycles,N);
	fclose(f);

}
void whistory::write_tally(unsigned tallynum){

//	//tallynum is unused at this point
//	float tally_err 	= 0;
//	float tally_err_sq 	= 0;
//	float this_square 	= 0;
//	float this_mean 	= 0;
//	float this_count 	= 0;
//	float Emin 			= 1e-11;
//	float Emax 			= 20.0;
//	float edge 			= 0.0;
//
//	// write tally values
//	std::string tallyname = filename+".tally";
//	FILE* tfile = fopen(tallyname.c_str(),"w");
//	for (int k=0;k<n_tally;k++){
//		this_count  	= (float) 	(N*n_cycles);//tally_count_total[k];
//		this_mean 		= 			tally_score_total[k];
//		this_square 	= 			tally_square_total[k];
//		tally_err_sq 	= 			(1.0/((this_count - 1.0))) * ( (this_count*this_square)/(this_mean*this_mean) -  1.0 ) ;
//		if(tally_err_sq>0.0){ 	
//			tally_err = sqrtf(tally_err_sq);}
//		else{					
//			tally_err = 0.0;}
//		fprintf(tfile,"%10.8E %10.8E %lu\n", this_mean/reduced_weight_total, tally_err, tally_count_total[k]);
//	}
//	fclose(tfile);
//
//	//write spacing, since there is n_tally+1 values, including it with the tally bins in a flat text messes with reading in as a matrix
//	std::string binsname = filename+".tallybins";
//	tfile = fopen(binsname.c_str(),"w");
//	for (int k=0;k<=n_tally;k++){
//		edge = Emin*powf((Emax/Emin), ((float)k)/ ((float)n_tally) );
//		fprintf(tfile,"%10.8E\n",edge);
//	}
//	fclose(tfile);
}
void whistory::prep_secondaries(){

//	// scan the yields to determine where the individual threads write into the done data
//	res = cudppScan( scanplan_int, d_scanned,  d_yield,  Ndataset );
//	if (res != CUDPP_SUCCESS){fprintf(stderr, "Error in scanning yield values\n");exit(-1);}
//
//	// compact the done data to know where to write
//	res = cudppCompact( compactplan, d_completed , (size_t*) d_num_completed , d_remap , d_done , Ndataset);
//	if (res != CUDPP_SUCCESS){fprintf(stderr, "Error in compacting done values\n");exit(-1);}

	//unsigned * tmp1 = new unsigned [Ndataset];
	//unsigned * tmp2 = new unsigned [Ndataset];
	//unsigned * tmp3 = new unsigned [Ndataset];
	//unsigned * tmp4 = new unsigned [Ndataset];
	//cudaMemcpy(tmp1,d_scanned,Ndataset*sizeof(unsigned),cudaMemcpyDeviceToHost);
	//cudaMemcpy(tmp2,d_completed,Ndataset*sizeof(unsigned),cudaMemcpyDeviceToHost);
	//cudaMemcpy(tmp3,d_done,Ndataset*sizeof(unsigned),cudaMemcpyDeviceToHost);
	//cudaMemcpy(tmp4,d_yield,Ndataset*sizeof(unsigned),cudaMemcpyDeviceToHost);
	//for(int k=0; k<Ndataset ; k++ ){printf("(tid,done,scanned,completed,yield) %u %u %u %u %u\n",k,tmp3[k],tmp1[k],tmp2[k],tmp4[k]);}

}
void whistory::remap_active(unsigned* num_active, unsigned* escatter_N, unsigned* escatter_start, unsigned* iscatter_N, unsigned* iscatter_start, unsigned* cscatter_N, unsigned* cscatter_start, unsigned* fission_N, unsigned* fission_start){

//	unsigned resamp_N = 0;
//	unsigned resamp_start = 0;
//
//	// sort key/value of rxn/tid
//	if(*num_active>1){
//		res = cudppRadixSort(radixplan, d_rxn, d_remap, *num_active );  //everything in 900s doesn't need to be sorted anymore
//		if (res != CUDPP_SUCCESS){fprintf(stderr, "Error in sorting reactions\n");exit(-1);}
//	}
//
//	// launch edge detection kernel, writes mapped d_edges array
//	reaction_edges(NUM_THREADS, *num_active, d_edges, d_rxn);
//
//	// calculate lengths and starting indicies for the blocks, 0 indicates not found
//	if (edges[0]){
//		*escatter_N 	= (edges[1]-edges[0])+1;
//		*escatter_start	= edges[0]-1;
//	}
//	else{
//		*escatter_N 	= 0;
//		*escatter_start	= 0;
//	}
//	if (edges[2]){
//		*iscatter_N 	= (edges[3]-edges[2])+1;
//		*iscatter_start	= edges[2]-1;
//	}
//	else{
//		*iscatter_N 	= 0;
//		*iscatter_start	= 0;
//	}
//	if (edges[4]){
//		*cscatter_N 	= (edges[5]-edges[4])+1;
//		*cscatter_start	= edges[4]-1;
//	}
//	else{
//		*cscatter_N 	= 0;
//		*cscatter_start	= 0;
//	}
//	if (edges[6]){
//		   resamp_N 	= (edges[7]-edges[6])+1;
//		   resamp_start	= edges[6]-1;
//	}
//	else{
//		   resamp_N 	= 0;
//		   resamp_start	= 0;
//	}
//	if (edges[8]){
//		 *fission_N 	= (edges[9]-edges[8])+1;
//		 *fission_start	= edges[8]-1;
//	}
//	else{
//		 *fission_N 	= 0;
//		 *fission_start	= 0;
//	}
//
//	//calculate total active, [10] is the starting index of >=811, so if==1, means that we are done!
//	*num_active=*escatter_N + *iscatter_N + *cscatter_N + resamp_N + *fission_N;   // fission includes multiplicity, have to reprocess
//
//	// debug
//	if(print_flag>=4){
//		//print
//		printf("num_active %u , edges[9] %u\n",*num_active,edges[9]);
//		printf("nactive = %u, edges %u %u %u %u %u %u %u %u %u %u %u \n",*num_active,edges[0],edges[1],edges[2],edges[3],edges[4],edges[5],edges[6],edges[7],edges[8],edges[9],edges[10]);
//		printf("escatter s %u n %u, iscatter s %u n %u, cscatter s %u n %u, resamp s %u n %u, fission s %u n %u \n\n",*escatter_start,*escatter_N,*iscatter_start,*iscatter_N,*cscatter_start,*cscatter_N,resamp_start,resamp_N, *fission_start, *fission_N);
//		if(dump_flag>=2){//dump
//			write_to_file(d_remap, d_rxn, N,"remap","w");
//		}
//	}
//
//	// ensure order
//	if(*iscatter_N>0){ assert(*iscatter_start >= *escatter_start);}
//	if(*cscatter_N>0){ assert(*cscatter_start >= *iscatter_start);}
//	if(resamp_N>0){    assert(   resamp_start >= *cscatter_start);}
//	if(*fission_N>0){  assert( *fission_start >=    resamp_start);}
//
//	// assert totals
//
//
//	// rezero edge vector (mapped, so is reflected on GPU)
//	edges[0]  = 0; 
//	edges[1]  = 0; 
//	edges[2]  = 0;
//	edges[3]  = 0;
//	edges[4]  = 0;
//	edges[5]  = 0;
//	edges[6]  = 0;
//	edges[7]  = 0;
//	edges[8]  = 0;
//	edges[9]  = 0;
//	edges[10] = 0;

}
void whistory::set_run_type(unsigned type_in){

	RUN_FLAG = type_in;

}
void whistory::set_run_type(std::string type_in){

	if(type_in.compare("fixed")==0){
		RUN_FLAG = 0;
	}
	else if(type_in.compare("criticality")==0){
		// check of there are fissile materials
		if(problem_geom.check_fissile()){
			//set flag to criticality
			RUN_FLAG = 1;
		}
		else{
			RUN_FLAG = 0;
			std::cout << "\e[1;31m" << "!!! No materials marked as fissile, criticality source mode rejected, will run in fixed source mode !!!" << "\e[m \n";
		}
	}
	else{
		std::cout << "Run type \"" << type_in << "\" not recognized." << "\n";
	}

}
void whistory::set_run_param(unsigned n_cycles_in, unsigned n_skip_in){

	n_skip = n_skip_in;
	n_cycles = n_cycles_in;

}
void whistory::set_device(unsigned dev_in){

	if (is_initialized){
		printf("Cannot change device after initialization.  Device already set to %u\n",compute_device);
		return;
	}

	//get number to make sure this is a valid device
	int 			n_devices;
	cudaGetDeviceCount(&n_devices);

	// set obj
	if(dev_in < n_devices){
		compute_device = dev_in;
	}
	else{
		std::cout << "!!!! Device " << dev_in << " does not exist.  Max devices is " << n_devices <<"\n";
	}

}
void whistory::device_report(){

	// vars
	cudaDeviceProp 	device_prop;
	int 		n_devices;
	int 		enabled_dev;
	float 		compute_cap;
	std::string 	con_string;
	std::string 	enabled_string;

	// get the number of devices
	cudaGetDeviceCount(&n_devices);
	cudaGetDevice(&enabled_dev);

	// loop over and print
	std::cout << "\e[1;32m" << "--- Compute Devices Present ---" << "\e[m \n";
	std::cout << "  --------------------------------------------------------------------------------------------------------------------\n";
	std::cout << "  | Device  | Model                       |  SMs  | Global Mem | SM Freq | Mem Freq | Compute Cap. | Concurrent Kern |" << "\n";
	std::cout << "  --------------------------------------------------------------------------------------------------------------------\n";
	for(unsigned k=0;k<n_devices;k++){
		cudaGetDeviceProperties(&device_prop,k);
		compute_cap = (float)device_prop.major + (float)device_prop.minor/10.0;
		con_string = "no ";
		enabled_string=' ';
		if(device_prop.concurrentKernels){con_string="yes";}
		if(k==enabled_dev){enabled_string='*';}
		printf(  "%1s | %2d      | %25s   |  %3d  | %6.4f  | %6.1f  | %6.1f   | %2.1f          | %4s            |\n", enabled_string.c_str(),k, device_prop.name, device_prop.multiProcessorCount, (float)device_prop.totalGlobalMem/(1024*1024), (float)device_prop.clockRate/1e3, (float)device_prop.memoryClockRate/1e3, compute_cap, con_string.c_str());
		std::cout << "  --------------------------------------------------------------------------------------------------------------------\n";
	}
		
}
void whistory::set_filename(std::string filename_in){

	filename = filename_in;

}
void whistory::set_acceration(std::string accel_in){

}
float whistory::get_time(){

	return ((float)clock())/((float)CLOCKS_PER_SEC);

}
void whistory::set_tally_cell(unsigned cell){

	// deprecated as of jan 12, 2016
	//tally_cell = cell;

}
void whistory::write_histories(unsigned iteration){

//	//allocate
//	unsigned*  done2;	
//	unsigned*  cellnum2;
//	unsigned*  matnum2;	
//	unsigned*  isonum2;	
//	unsigned*  yield2;
//	unsigned*  rxn2;
//	unsigned*  dex2;
//	spatial_data* space2;
//	float* E2;
//	done2 		= new unsigned [N];
//	cellnum2	= new unsigned [N];
//	matnum2		= new unsigned [N];
//	isonum2		= new unsigned [N];
//	yield2 		= new unsigned [N];
//	rxn2		= new unsigned [N];	
//	dex2 		= new unsigned [N];
//	space2 		= new spatial_data [N];
//	E2 			= new float [N];
//
//	// print actions
//	std::string histfile_name = "history_file";
//	char numstr[5];
//	sprintf(numstr,"%u",iteration);
//	histfile_name += numstr;
//	printf("Writing to \"%s\"\n",histfile_name.c_str());
//
//	// open file, appending to a new file.
//	FILE* history_file = fopen(histfile_name.c_str(),"w");
//
//	// write iteration delimiter
//	//fprintf(history_file,"==== ITERATION %u ====\n",iteration);
//
//	// copy gemetrical (positions / directions)
//	cudaMemcpy(space2,d_space,N*sizeof(spatial_data),cudaMemcpyDeviceToHost);
//
//	// copy energies
//	cudaMemcpy(E2,d_E,N*sizeof(float),cudaMemcpyDeviceToHost);
//
//	// copy cell numbeer
//	cudaMemcpy(cellnum2,d_cellnum,N*sizeof(unsigned),cudaMemcpyDeviceToHost);
//
//	// copy material numbers
//	cudaMemcpy(matnum2,d_matnum,N*sizeof(unsigned),cudaMemcpyDeviceToHost);
//
//	// copy done flags
//	cudaMemcpy(done2,d_done,N*sizeof(unsigned),cudaMemcpyDeviceToHost);
//
//	// copy isotope numbers
//	cudaMemcpy(isonum2,d_isonum,N*sizeof(unsigned),cudaMemcpyDeviceToHost);
//
//	// copy yields
//	cudaMemcpy(yield2,d_yield,N*sizeof(unsigned),cudaMemcpyDeviceToHost);
//
//	// copy reaction numbers
//	cudaMemcpy(rxn2,d_rxn,N*sizeof(unsigned),cudaMemcpyDeviceToHost);
//
//	// copy index numbers
//	cudaMemcpy(dex2,d_index,N*sizeof(unsigned),cudaMemcpyDeviceToHost);
//
//	// sync device before write and return
//	cudaDeviceSynchronize();
//
//	// write history data to file
//	fprintf(history_file,"a=zeros(%u,7)\n",N);
//	for(unsigned k=0;k<N;k++){
//		//fprintf(history_file,"tid %u (x,y,z) %8.6E %8.6E %8.6E (x,y,z)-hat %8.6E %8.6E %8.6E surf_dist %8.6E macro_t %8.6E enforce_BC %u E %8.6E cellnum %u matnum %u isonum %u rxn %u done %u yield %u\n",k,space2[k].x,space2[k].y,space2[k].z,space2[k].xhat,space2[k].yhat,space2[k].zhat,space2[k].surf_dist,space2[k].macro_t,space2[k].enforce_BC,E2[k],cellnum2[k],matnum2[k],isonum2[k],rxn2[k],done2[k],yield2[k]);
//		fprintf(history_file,"a(%u,1:7)=[%8.6E,%8.6E,%8.6E,%8.6E,%u,%u,%u];\n",k+1,space2[k].x,space2[k].y,space2[k].z,E2[k],rxn2[k],yield2[k],dex2[k]);
//	}
// 	fclose(history_file);
//
// 	//deallocate so can be alloaed again next time
// 	delete done2 	;
//	delete cellnum2	;
//	delete matnum2	;
//	delete isonum2	;
//	delete yield2 	;
//	delete rxn2		;
//	delete space2 	;
//	delete E2 		;
//	delete dex2		;

}
void whistory::set_print_level(unsigned level){
	print_flag = level;
}
void whistory::set_dump_level(unsigned level){
	dump_flag = level;
}
void whistory::bin_fission_points( spatial_data * d_space, unsigned N){

	// init
	unsigned res_x = 300;
	unsigned res_y = 300;
	unsigned dex_x, dex_y;
	float width_x = 1.41421356*( outer_cell_dims[3] - outer_cell_dims[0]);
	float width_y = 1.41421356*( outer_cell_dims[4] - outer_cell_dims[1]);
	float x_min   = 1.41421356 * outer_cell_dims[0];
	float y_min   = 1.41421356 * outer_cell_dims[1];

	// copy from device
	spatial_data * positions_local = new spatial_data[N];
	cudaMemcpy(positions_local,dh_particles.space,N*sizeof(spatial_data),cudaMemcpyDeviceToHost);

	// bin to grid
	for(unsigned i=0; i<N; i++){
		dex_x = unsigned(  ( positions_local[i].x - x_min )*res_x/width_x  );
		dex_y = unsigned(  ( positions_local[i].y - y_min )*res_y/width_y  );
		fiss_img[ dex_y*res_x + dex_x ] += 1; 
	}

	// free
	delete positions_local;

}
void whistory::write_fission_points(){

	unsigned res_x = 300;
	unsigned res_y = 300;
	size_t x, y;
	long unsigned minnum = 0; 
	long unsigned maxnum = 0;
	float * colormap = new float[3];
	png::image< png::rgb_pixel > image(res_x, res_y);

	// find maximum
	for ( y = 0; y < res_y; ++y)
	{
	    for ( x = 0; x < res_x; ++x)
	    {
	    	if(fiss_img[y*res_x+x] >= maxnum){
	    		maxnum = fiss_img[y*res_x+x];
	    	}
	    }
	}

	// make image via colormap
	for ( y = 0; y < res_y; ++y)
	{
	    for ( x = 0; x < res_x; ++x)
	    {
	    	hot2(colormap,fiss_img[y*res_x+x],minnum,maxnum);
	        image[image.get_height()-1-y][x] = png::rgb_pixel(colormap[0],colormap[1],colormap[2]);
	    }
	}

	// write
	std::string fiss_name=filename;
	fiss_name.append(".fission_points.png");
	image.write(fiss_name.c_str());
	printf("  Fission point image written to %s.\n",fiss_name.c_str());


	// clear
	delete fiss_img;
	delete colormap;

}
void whistory::plot_geom(std::string type_in){

	if(is_initialized!=1){
		printf("  ! geometry plotting must be done AFTER init, skipping.\n");
		return;
	}

	printf("\e[1;32mPlotting Geometry... \e[m \n");

	srand (time(NULL));

	// type logic
	char this_filename[50];
	unsigned* type_array;
	unsigned width, height, width_in, height_in, index, N_plot;
	float dx, dy, dz;
	float aspect, mu, theta;
	float resolution = 1024;
	float pi = 3.14159;
	spatial_data * positions_local = new spatial_data[N];
	unsigned * image_local = new unsigned[N];
	unsigned minnum, maxnum;
	if (type_in.compare("cell")==0){
		printf("  color set to \e[1;31mCELL\e[m\n");
		type_array = dh_particles.cellnum;
		minnum = problem_geom.get_minimum_cell();
		maxnum = problem_geom.get_maximum_cell();
	}
	else if (type_in.compare("material")==0){
		printf("  color set to \e[1;31mMATERIAL\e[m\n");
		type_array = dh_particles.matnum;
		minnum = problem_geom.get_minimum_material();
		maxnum = problem_geom.get_maximum_material();
	}
	else{
		printf("  ! plot type '%s' unknown, skipping.\n",type_in.c_str());
		return;
	}

	// get outer cell dims
	float xmin = outer_cell_dims[0] * 1.41421356;
	float ymin = outer_cell_dims[1] * 1.41421356;
	float zmin = outer_cell_dims[2] * 1.41421356;
	float xmax = outer_cell_dims[3] * 1.41421356;
	float ymax = outer_cell_dims[4] * 1.41421356;
	float zmax = outer_cell_dims[5] * 1.41421356;
	
	//
	// xy
	//
	aspect = (xmax-xmin)/(ymax-ymin);
	printf("  xy aspect ratio of outer cell: %6.4f\n",aspect);
	width_in  = resolution*aspect;
	height_in = resolution;
	if (width_in*height_in > N){
		width  = sqrtf(N*aspect); 
		height = sqrtf(N/aspect);
		printf("  !resolution reduced by dataset size ->");
		printf("  (height,width) (%u,%u) px\n",height,width);
	}else{
		width = width_in;
		height = height_in;
		printf(" (height,width) (%u,%u) px\n",height,width);
	}
	N_plot = width*height;
	dx = (xmax-xmin)/width;
	dy = (ymax-ymin)/height;
	float r1, r2;
	for(int j=0;j<height;j++){
		for(int k=0;k<width;k++){
			r1 = (rand()/RAND_MAX);
			r2 = (rand()/RAND_MAX);
			mu = 2.0*r1-1.0;
			theta = 2.0*pi*r2;
			index = j * width + k;
			positions_local[index].x = xmin + dx/2 + k*dx;
			positions_local[index].y = ymin + dy/2 + j*dy;
			positions_local[index].z = 0.0;
			positions_local[index].xhat = sqrtf(1.0-mu*mu) * cosf( theta ); 
			positions_local[index].yhat = sqrtf(1.0-mu*mu) * sinf( theta ); 
			positions_local[index].zhat =       mu; 
			positions_local[index].surf_dist = 9999999999.9; 
		}
	}

	// copy starting positions data to pointer
	cudaMemcpy(dh_particles.space,positions_local,N_plot*sizeof(spatial_data),cudaMemcpyHostToDevice);
	
	// trace with whereami?
	trace(4, N_plot);
	
	//copy to local buffer
	cudaMemcpy(image_local,type_array,N_plot*sizeof(unsigned),cudaMemcpyDeviceToHost);

	// make image
	png::image< png::rgb_pixel > image(width, height);
	float * colormap = new float[3];
	for (size_t y = 0; y < image.get_height(); ++y)
	{
	    for (size_t x = 0; x < image.get_width(); ++x)
	    {
	    	make_color(colormap,image_local[y*width+x],minnum,maxnum);
	        image[image.get_height()-1-y][x] = png::rgb_pixel(colormap[0],colormap[1],colormap[2]);
	    }
	}

	sprintf(this_filename,"%s-xy.png",filename.c_str());
	image.write(this_filename);
	printf("  Written to %s.\n",this_filename);

	//
	// xz
	//
	aspect = (xmax-xmin)/(zmax-zmin);
	printf("  xz aspect ratio of outer cell: %6.4f\n",aspect);
	width_in  = resolution*aspect;
	height_in = resolution;
	if (width_in*height_in > N){
		width  = sqrtf(N*aspect); 
		height = sqrtf(N/aspect);
		printf("  !resolution reduced by dataset size ->");
		printf("  (height,width) (%u,%u) px\n",height,width);
	}else{
		width = width_in;
		height = height_in;
		printf(" (height,width) (%u,%u) px\n",height,width);
	}
	N_plot = width*height;
	dx = (xmax-xmin)/width;
	dz = (zmax-zmin)/height;
	for(int j=0;j<height;j++){
		for(int k=0;k<width;k++){
			r1 = (rand()/RAND_MAX);
			r2 = (rand()/RAND_MAX);
			mu = 2.0*r1-1.0;
			theta = 2.0*pi*r2;
			index = j * width + k;
			positions_local[index].x = xmin + dx/2 + k*dx;
			positions_local[index].y = 0.0;
			positions_local[index].z = zmin + dz/2 + j*dz;
			positions_local[index].xhat = sqrtf(1.0-mu*mu) * cosf( theta ); 
			positions_local[index].yhat = sqrtf(1.0-mu*mu) * sinf( theta ); 
			positions_local[index].zhat =       mu; 
			positions_local[index].surf_dist = 9999999999.9; 
		}
	}

	// copy starting positions data to pointer
	cudaMemcpy(dh_particles.space,positions_local,N_plot*sizeof(spatial_data),cudaMemcpyHostToDevice);
	
	// trace with whereami?
	trace(4, N_plot);
	
	//copy to local buffer
	cudaMemcpy(image_local,type_array,N_plot*sizeof(unsigned),cudaMemcpyDeviceToHost);

	// make image
	image = png::image< png::rgb_pixel > (width, height);
	for (size_t y = 0; y < image.get_height(); ++y)
	{
	    for (size_t x = 0; x < image.get_width(); ++x)
	    {
	    	make_color(colormap,image_local[y*width+x],minnum,maxnum);
	        image[image.get_height()-1-y][x] = png::rgb_pixel(colormap[0],colormap[1],colormap[2]);
	    }
	}

	sprintf(this_filename,"%s-xz.png",filename.c_str());
	image.write(this_filename);
	printf("  Written to %s.\n",this_filename);

	//
	// yz
	//
	aspect = (ymax-ymin)/(zmax-zmin);
	printf("  yz aspect ratio of outer cell: %6.4f\n",aspect);
	width_in  = resolution*aspect;
	height_in = resolution;
	if (width_in*height_in > N){
		width  = sqrtf(N*aspect); 
		height = sqrtf(N/aspect);
		printf("  !resolution reduced by dataset size ->");
		printf("  (height,width) (%u,%u) px\n",height,width);
	}else{
		width = width_in;
		height = height_in;
		printf(" (height,width) (%u,%u) px\n",height,width);
	}
	N_plot = width*height;
	dy = (ymax-ymin)/width;
	dz = (zmax-zmin)/height;
	for(int j=0;j<height;j++){
		for(int k=0;k<width;k++){
			r1 = (rand()/RAND_MAX);
			r2 = (rand()/RAND_MAX);
			mu = 2.0*r1-1.0;
			theta = 2.0*pi*r2;
			index = j * width + k;
			positions_local[index].x = 0.0;
			positions_local[index].y = ymin + dy/2 + k*dy;
			positions_local[index].z = zmin + dz/2 + j*dz;
			positions_local[index].xhat = sqrtf(1.0-mu*mu) * cosf( theta ); 
			positions_local[index].yhat = sqrtf(1.0-mu*mu) * sinf( theta ); 
			positions_local[index].zhat =       mu; 
			positions_local[index].surf_dist = 9999999999.9; 
		}
	}

	// copy starting positions data to pointer
	cudaMemcpy(dh_particles.space,positions_local,N_plot*sizeof(spatial_data),cudaMemcpyHostToDevice);
	
	// trace with whereami?
	trace(4, N_plot);
	
	//copy to local buffer
	cudaMemcpy(image_local,type_array,N_plot*sizeof(unsigned),cudaMemcpyDeviceToHost);

	// make image
	image = png::image< png::rgb_pixel > (width, height);
	for (size_t y = 0; y < image.get_height(); ++y)
	{
	    for (size_t x = 0; x < image.get_width(); ++x)
	    {
	    	make_color(colormap,image_local[y*width+x],minnum,maxnum);
	        image[image.get_height()-1-y][x] = png::rgb_pixel(colormap[0],colormap[1],colormap[2]);
	    }
	}

	sprintf(this_filename,"%s-yz.png",filename.c_str());
	image.write(this_filename);
	printf("  Written to %s.\n",this_filename);


	// clear
	delete image_local;
	delete colormap;
	delete positions_local;

}
void whistory::make_color(float* color, unsigned x, unsigned min, unsigned max){
	// red linear blue linear green sin colormap
	if (x==4294967295){ //miss value
		color[0]=0.0;
		color[1]=0.0;
		color[2]=0.0;
	}
	else{
		float normed_value = (float) (x-min+1)/(max+2-min);
		color[0] = normed_value;              //red
		color[1] = sin(normed_value*3.14159); //green
		color[2] = 1.0-normed_value;          //blue
	}
	
	//bring up to 256?
	color[0]=color[0]*256;
	color[1]=color[1]*256;
	color[2]=color[2]*256;

}
void whistory::hot2(float* color, long unsigned val, long unsigned min, long unsigned max){

	float val_norm = (float) (val-min)/(max-min);
	if(val_norm>1.0){
		printf("val_norm>1, min %lu max %lu val %lu\n",min,max,val);
	}

	if( val_norm < 0.025){
		color[0]=0.3*val_norm/0.025;
		color[1]=0.0;
		color[2]=0.0;
	}
	else if(val_norm < 0.3 & val_norm >= 0.025){
		color[0]= (0.7/0.275)*(val_norm-0.025) + 0.3;
		color[1]=0.0;
		color[2]=0.0;
	}
	else if(val_norm < 0.7 & val_norm >= 0.3){
		color[0]=1.0;
		color[1]=1.0/(0.4)*(val_norm-0.3) + 0.0;
		color[2]=0.0;
	}
	else {
		color[0]=1.0;
		color[1]=1.0;
		color[2]=1.0/(0.3)*(val_norm-0.7) + 0.0;
	}

	if(color[0]>1){
		printf("red=%6.4E!  val_norm = %6.4E\n",color[0],val_norm);
	}
	if(color[1]>1){
		printf("green>%6.4E  val_norm = %6.4E!\n",color[1],val_norm);
	}
	if(color[2]>1){
		printf("blue>%6.4E  val_norm = %6.4E!\n",color[2],val_norm);
	}

	float total_norm = sqrtf(color[0]*color[0]+color[1]*color[1]+color[2]*color[2]);

	//bring up to 256?
	color[0] = color[0]*256/total_norm;
	color[1] = color[1]*256/total_norm;
	color[2] = color[2]*256/total_norm;

}
void whistory::nonzero(float* color, unsigned val, unsigned min, unsigned max){

	// white if not zero
	if( val > 0 ){
		//printf("nonzero\n");
		color[0]=0.57735026918963;
		color[1]=0.57735026918963;
		color[2]=0.57735026918963;
	}
	else{
		//printf("ZERO\n");
		color[0]=0.0;
		color[1]=0.0;
		color[2]=0.0;
	}

	//bring up to 256?
	color[0] = color[0]*256;
	color[1] = color[1]*256;
	color[2] = color[2]*256;

}