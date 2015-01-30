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
#include <cudpp_hash.h>
#include <Python.h>
#include <png++/png.hpp>
#include "datadef.h"
#include "primitive.h"
#include "wgeometry.h"
#include "optix_stuff.h"
#include "warp_cuda.h"
#include "whistory.h"

optix_stuff optix_obj;
//wgeometry geom_obj;

whistory::whistory(unsigned Nin, wgeometry problem_geom_in){
	//clear device
	cudaDeviceReset();
	// do problem gemetry stuff first
	problem_geom = problem_geom_in;
	// set tally vector length
	n_tally = 1024;
	RUN_FLAG = 1;
	// keff stuff
	keff_sum = 0.0;
	keff2_sum = 0.0;
	keff_err = 0.0;
	// device data stuff
	N = Nin;
	Ndataset = Nin * 5;
	n_qnodes = 0;
	compute_device = 0;	
	reduced_yields_total = 0;
	accel_type = "Sbvh";
	for (int i = 0; i < 4; ++i){
		cudaStreamCreate(&stream[i]);
	}
}
void whistory::init(){
	// set device first
	cudaSetDevice(compute_device);
	// init optix stuff second
	optix_obj.N=Ndataset;
	optix_obj.stack_size_multiplier=1;
	optix_obj.set_image_type("cell");
	optix_obj.init(problem_geom,compute_device,accel_type);
	optix_obj.print();
	// CUDA stuff
	std::cout << "\e[1;32m" << "Dataset size is "<< N << "\e[m \n";
	NUM_THREADS = 256;
	RNUM_PER_THREAD = 1;
	blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;
	std::cout << "\e[1;32m" << "Compute device set to "<< compute_device << ". Available devices shown below:" <<"\e[m \n";
	device_report();
	cudaDeviceSetLimit(cudaLimitPrintfFifoSize, (size_t) 10*1048576 );
	//device data
	d_space 	= (source_point*) optix_obj.positions_ptr;
	d_cellnum 	= (unsigned*)     optix_obj.cellnum_ptr;
	d_matnum 	= (unsigned*)     optix_obj.matnum_ptr;
	d_rxn 		= (unsigned*)     optix_obj.rxn_ptr;
	d_done 		= (unsigned*)     optix_obj.done_ptr;
	d_remap 	= (unsigned*)     optix_obj.remap_ptr;
	cudaMalloc( &d_xs_length_numbers	, 6*sizeof(unsigned) );		 
	cudaMalloc( &d_E 			, Ndataset*sizeof(float) );
	cudaMalloc( &d_Q 			, Ndataset*sizeof(float) );
	cudaMalloc( &d_rn_bank  		, Ndataset*RNUM_PER_THREAD*sizeof(float) );
	cudaMalloc( &d_isonum   		, Ndataset*sizeof(unsigned) );
	cudaMalloc( &d_yield			, Ndataset*sizeof(unsigned) );
	cudaMalloc( &d_index			, Ndataset*sizeof(unsigned) );
	cudaMalloc( &d_tally_score  		, n_tally*sizeof(float));
	cudaMalloc( &d_tally_square  		, n_tally*sizeof(float));
	cudaMalloc( &d_tally_count  		, n_tally*sizeof(unsigned));
	cudaMalloc( &d_reduced_yields 		, 1*sizeof(unsigned));
	cudaMalloc( &d_reduced_done 		, 1*sizeof(unsigned));
	cudaMalloc( &d_valid_result		, Ndataset*sizeof(unsigned));
	cudaMalloc( &d_valid_N			, 1*sizeof(unsigned));
	//cudaMalloc( &d_remap			, Ndataset*sizeof(unsigned));
	cudaMalloc( &d_fissile_points		, Ndataset*sizeof(source_point));
	cudaMalloc( &d_fissile_energy       	, Ndataset*sizeof(float));
	cudaMalloc( &d_mask 			, Ndataset*sizeof(unsigned));
	cudaMalloc( &d_completed 		, Ndataset*sizeof(unsigned));	
	cudaMalloc( &d_scanned 			, Ndataset*sizeof(unsigned));
	cudaMalloc( &d_num_completed 		, 1*sizeof(unsigned));
	cudaMalloc( &d_active 			, Ndataset*sizeof(unsigned));
	cudaMalloc( &d_rxn_remap		, Ndataset*sizeof(unsigned));
	cudaMalloc( &d_num_active 		, 1*sizeof(unsigned));
	cudaMalloc( &d_zeros			, Ndataset*sizeof(unsigned));
	// host data stuff
	//xs_length_numbers 	= new unsigned [6];
	space 			= new source_point 	[Ndataset];
	E 				= new float 		[Ndataset];
	Q 				= new float 		[Ndataset];
	rn_bank  		= new unsigned 		[Ndataset*RNUM_PER_THREAD];
	tally_score 		= new float 		[n_tally];
	tally_square 		= new float 		[n_tally];
	tally_count 		= new unsigned 		[n_tally];
	tally_score_total	= new double 		[n_tally];
	tally_square_total	= new double 		[n_tally];
	tally_count_total	= new long unsigned	[n_tally];
	index     		= new unsigned 		[Ndataset];
	cellnum 		= new unsigned 		[Ndataset];
	matnum 			= new unsigned 		[Ndataset];
	rxn 			= new unsigned 		[Ndataset];
	done 			= new unsigned 		[Ndataset];
	isonum   		= new unsigned 		[Ndataset];
	yield	   		= new unsigned 		[Ndataset];
	remap 			= new unsigned 		[Ndataset];
	zeros 			= new unsigned 		[Ndataset];
	ones 			= new unsigned 		[Ndataset];
	// init counters to 0
	total_bytes_scatter = 0;
	total_bytes_energy  = 0;
	//copy any info needed
	memcpy(outer_cell_dims,optix_obj.outer_cell_dims,6*sizeof(float));
	outer_cell = optix_obj.get_outer_cell();
	outer_cell_type = optix_obj.get_outer_cell_type();
	xs_isotope_string = problem_geom.isotope_list;
	//  map edge array
	n_edges = 11;
	cudaHostAlloc(&edges,n_edges*sizeof(unsigned),cudaHostAllocMapped);
	cudaHostGetDevicePointer(&d_edges,edges,0);
	// init host values
	filename = "warp";
	init_host();
	init_RNG();
	init_CUDPP();
	load_cross_sections();
	//create_quad_tree();
	copy_data_to_device();
	printf("Done with init\n");
}
whistory::~whistory(){
	//cudaFree( d_xs_length_numbers 	);
	//cudaFree( d_xs_MT_numbers_total );
	//cudaFree( d_xs_MT_numbers 		);
	//cudaFree( d_xs_data_MT 			);
	//cudaFree( d_xs_data_main_E_grid );
	//cudaFree( d_xs_data_scatter     );
	//cudaFree( d_xs_data_energy      );
	//cudaFree( d_tally_score 		);
    //cudaFree( d_tally_count 		);
    //cudaFree( d_xs_data_Q    		);
	//cudaFree( d_index   );
	//cudaFree( d_E         );
	//cudaFree( d_Q         );
	//cudaFree( d_rn_bank   );
	//cudaFree( d_isonum    );
	//cudaFree( d_yield     );
	//cudaFree( d_awr_list);
/*	delete xs_length_numbers; 
	delete xs_MT_numbers_total;
	delete xs_MT_numbers;
	delete xs_data_MT;
	delete xs_data_main_E_grid;
	delete xs_data_Q;
	delete space;
	delete index;
	delete awr_list;
	delete E;
	delete Q;
	delete rn_bank;
	delete cellnum;
	delete matnum;
	delete rxn;
	delete done;
	delete isonum;
	delete yield; 
	delete tally_count;
	delete tally_score;
	delete zeros;
	delete ones;
	// for loops to deallocate everything in the pointer arrays
	for (int j=0 ; j < MT_columns ; j++){  //start after the total xs and total abs vectors
		for (int k=0 ; k < MT_rows ; k++){
			// scatter
			//std::cout << "j,k = " << j << ", " << k << " colums,rows = " << MT_columns << ", " << MT_rows << "\n";
			float * this_pointer =   xs_data_scatter     [k*MT_columns + j];
			float * cuda_pointer =   xs_data_scatter_host[k*MT_columns + j];
			if(this_pointer!=NULL){
				while(xs_data_scatter[(k+1)*MT_columns + j ]==this_pointer){
					k++; //push k to the end of the copies so don't try to free it twice
				}
				//std::cout << "j,k " << j << ", " << k << " - " ;
				//std::cout << "freeing " << this_pointer << " " << cuda_pointer << "\n";
				delete this_pointer;
				//cudaFree(cuda_pointer);
			}
		}
	}
	//delete pointer arrays themselves
	delete xs_data_scatter;
	delete xs_data_energy;
	delete xs_data_scatter_host;
	delete xs_data_energy_host;
*/}
void whistory::init_host(){

	for(int k=0;k<N;k++){
		remap[k]		= k;
		zeros[k]		= 0;
		ones[k]			= 1;
		space[k].x 		= 0.0;
		space[k].y 		= 0.0;
		space[k].z 		= 0.0;
		space[k].xhat 		= 0.0;
		space[k].yhat 		= 0.0;
		space[k].zhat 		= 0.0;
		space[k].surf_dist 	= 10000.0;
		space[k].macro_t 	= 0.0;
		space[k].enforce_BC     = 0;
		E[k]			= 2.5;
		Q[k]			= 0.0;
		cellnum[k]		= 0;
		matnum[k]		= 0;
		rxn[k]			= 0;
		done[k]			= 0;
		isonum[k]		= 0;
		yield[k]		= 0;
	}
	for(int k=N;k<Ndataset;k++){
		remap[k]		= k;
		zeros[k]		= 0;
		ones[k]			= 1;
		space[k].x 		= 0.0;
		space[k].y 		= 0.0;
		space[k].z 		= 0.0;
		space[k].xhat 		= 0.0;
		space[k].yhat 		= 0.0;
		space[k].zhat 		= 0.0;
		space[k].surf_dist 	= 10000.0;
		space[k].macro_t 	= 0.0;
		space[k].enforce_BC = 0;
		E[k]			= 0.0;
		Q[k]			= 0.0;
		cellnum[k]		= 0;
		matnum[k]		= 0;
		rxn[k]			= 0;
		done[k]			= 1;
		isonum[k]		= 0;
		yield[k]		= 0;
	}
	for(int k=0;k<n_tally;k++){
		tally_score_total[k]=0.0;
		tally_square_total[k]=0.0;
		tally_count_total[k]=0;
	}

}
void whistory::init_RNG(){
	unsigned seed = time( NULL );
	std::cout << "\e[1;32m" << "Initializing random number bank on device using MTGP32 with seed of " << seed << "..." << "\e[m \n";
	curandCreateGenerator( &rand_gen , CURAND_RNG_PSEUDO_MTGP32 );  //mersenne twister type
	curandSetPseudoRandomGeneratorSeed( rand_gen , 123456789ULL );
	curandGenerate( rand_gen , d_rn_bank , Ndataset * RNUM_PER_THREAD );
	cudaMemcpy(rn_bank , d_rn_bank , Ndataset * RNUM_PER_THREAD *sizeof(unsigned) , cudaMemcpyDeviceToHost); // copy bank back to keep seeds
}
void whistory::update_RNG(){

	curandGenerate( rand_gen , d_rn_bank , Ndataset * RNUM_PER_THREAD );

}
void whistory::init_CUDPP(){
	
	std::cout << "\e[1;32m" << "Initializing CUDPP..." << "\e[m \n";
	// global objects
	res = cudppCreate(&theCudpp);
	if (res != CUDPP_SUCCESS){fprintf(stderr, "Error initializing CUDPP Library.\n");}
	
	std::cout << "  configuring compact..." << "\n";
	// sort stuff
	compact_config.op = CUDPP_ADD;
	compact_config.datatype = CUDPP_INT;
	compact_config.algorithm = CUDPP_COMPACT;
	compact_config.options = CUDPP_OPTION_FORWARD;
	res = cudppPlan(theCudpp, &compactplan, compact_config, Ndataset, 1, 0);
	if (CUDPP_SUCCESS != res){printf("Error creating CUDPPPlan for compact\n");exit(-1);}

	std::cout << "  configuring reduction..." << "\n";
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

	std::cout << "  configuring scan..." << "\n";
	// int reduction stuff
	scan_int_config.op = CUDPP_ADD;
	scan_int_config.datatype = CUDPP_INT;
	scan_int_config.algorithm = CUDPP_SCAN;
	scan_int_config.options = CUDPP_OPTION_EXCLUSIVE;
	res = cudppPlan(theCudpp, &scanplan_int, scan_int_config, Ndataset, 1, 0);
	if (CUDPP_SUCCESS != res){printf("Error creating CUDPPPlan for scan\n");exit(-1);}

	std::cout << "  configuring radix sort..." << "\n";
	// int reduction stuff
	radix_config.algorithm = CUDPP_SORT_RADIX;
	radix_config.datatype = CUDPP_UINT;
	radix_config.options = CUDPP_OPTION_KEY_VALUE_PAIRS;
	res = cudppPlan(theCudpp, &radixplan, radix_config, Ndataset, 1, 0);
	if (CUDPP_SUCCESS != res){printf("Error creating CUDPPPlan for radix sort\n");exit(-1);}

}
unsigned whistory::reduce_yield(){

	unsigned reduced_yields;

	//printf("reducing yield for %u elements\n",N);
	//cudaMemcpy(yield,d_yield,N*sizeof(unsigned),cudaMemcpyDeviceToHost);
	//for(unsigned k=0;k<N;k++){printf("yield[%u]=%u\n",k,yield[k]);}

	res = cudppReduce(reduplan_int, d_reduced_yields, d_yield, N);
	if (res != CUDPP_SUCCESS){fprintf(stderr, "Error in reducing yield values\n");exit(-1);}
	cudaMemcpy(&reduced_yields, d_reduced_yields, 1*sizeof(unsigned), cudaMemcpyDeviceToHost);

	return reduced_yields;

}
void whistory::accumulate_keff(unsigned converged, unsigned iteration, double* keff, float* keff_cycle){

	float this_count, this_square, this_mean, keff_err2;

	long unsigned reduced_yields = reduce_yield();

	*keff_cycle = reduced_yields / ( (float) N );

	if(converged){
		reduced_yields_total += reduced_yields;
		*keff       = reduced_yields_total / ( (double) (iteration + 1) * N ) ;
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

	//copy data to host
	cudaMemcpy(tally_score,  d_tally_score,  n_tally*sizeof(float),    cudaMemcpyDeviceToHost);
	cudaMemcpy(tally_square, d_tally_square, n_tally*sizeof(float),    cudaMemcpyDeviceToHost);
	cudaMemcpy(tally_count,  d_tally_count,  n_tally*sizeof(unsigned), cudaMemcpyDeviceToHost);

	//zero out vectors
	cudaMemcpy(d_tally_score,  d_zeros, n_tally*sizeof(float),    cudaMemcpyDeviceToDevice);
	cudaMemcpy(d_tally_square, d_zeros, n_tally*sizeof(float),    cudaMemcpyDeviceToDevice);
	cudaMemcpy(d_tally_count,  d_zeros, n_tally*sizeof(unsigned), cudaMemcpyDeviceToDevice);

	//perform sums on 64bit host side values
	for(unsigned k=0 ; k<n_tally ; k++){
		tally_score_total[k] 	+=  tally_score[k];
		tally_square_total[k]	+=  tally_square[k];
		tally_count_total[k] 	+=  tally_count[k];
		//printf("score %10.8E %10.8E\n",tally_score[k],tally_score_total[k]);
	}

	
}
unsigned whistory::reduce_done(){

	unsigned reduced_done = 0;

	printf("reducing done\n");
	res = cudppReduce(reduplan_int, d_reduced_done, d_done, Ndataset);
	if (res != CUDPP_SUCCESS){fprintf(stderr, "Error in reducing done values\n");exit(-1);}
	cudaMemcpy(&reduced_done, d_reduced_done, 1*sizeof(unsigned), cudaMemcpyDeviceToHost);

	return reduced_done;

}
void whistory::copy_data_to_device(){

	float * this_pointer;
	float * cuda_pointer;
	int vlen, next_vlen;
	float * temp = new float [128];
	for(int g=0;g<128;g++){temp[g]=123456789;}

	std::cout << "\e[1;32m" << "Copying data to device (number?)..." << "\e[m \n";

	// copy history data
	std::cout << "  History data... ";
	cudaMemcpy( d_space,		space,		Ndataset*sizeof(source_point),	cudaMemcpyHostToDevice );
	cudaMemcpy( d_E,			E,			Ndataset*sizeof(float),		cudaMemcpyHostToDevice );
	cudaMemcpy( d_Q,    		Q,			Ndataset*sizeof(float),		cudaMemcpyHostToDevice );
	cudaMemcpy( d_done,			done,		Ndataset*sizeof(unsigned),	cudaMemcpyHostToDevice );
	cudaMemcpy( d_cellnum,		cellnum,	Ndataset*sizeof(unsigned),	cudaMemcpyHostToDevice );
	cudaMemcpy( d_matnum,		matnum,		Ndataset*sizeof(unsigned),	cudaMemcpyHostToDevice );
	cudaMemcpy( d_isonum,		isonum,		Ndataset*sizeof(unsigned),	cudaMemcpyHostToDevice );
	cudaMemcpy( d_yield,		yield,		Ndataset*sizeof(unsigned),	cudaMemcpyHostToDevice );
	cudaMemcpy( d_rxn,			rxn,		Ndataset*sizeof(unsigned),	cudaMemcpyHostToDevice );
    cudaMemcpy( d_remap, 		remap,    	Ndataset*sizeof(unsigned),	cudaMemcpyHostToDevice );
    cudaMemcpy( d_active,		remap,		Ndataset*sizeof(unsigned),	cudaMemcpyHostToDevice );
    cudaMemcpy( d_zeros,		zeros,		Ndataset*sizeof(unsigned),	cudaMemcpyHostToDevice );
    std::cout << "Done.\n";
    std::cout << "  Unionized cross sections... ";
    // copy xs_data,  0=isotopes, 1=main E points, 2=total numer of reaction channels
    cudaMemcpy( d_xs_length_numbers,     xs_length_numbers,     3*sizeof(unsigned),						  cudaMemcpyHostToDevice );
    cudaMemcpy( d_xs_MT_numbers_total,   xs_MT_numbers_total,   xs_length_numbers[0]*sizeof(unsigned),			  cudaMemcpyHostToDevice );
    cudaMemcpy( d_xs_MT_numbers,	     xs_MT_numbers,	    (xs_length_numbers[2]+xs_length_numbers[0])*sizeof(unsigned), cudaMemcpyHostToDevice );
    cudaMemcpy( d_xs_data_MT,	     xs_data_MT,	    MT_rows*MT_columns *sizeof(float), 				  cudaMemcpyHostToDevice );
	cudaMemcpy( d_xs_data_main_E_grid,   xs_data_main_E_grid,   xs_length_numbers[1]*sizeof(float),				  cudaMemcpyHostToDevice );
	cudaMemcpy( d_awr_list, 	     awr_list,   	    xs_length_numbers[0]*sizeof(float),				  cudaMemcpyHostToDevice );
	cudaMemcpy( d_material_list,         material_list,         n_materials*sizeof(unsigned), 				  cudaMemcpyHostToDevice );
	cudaMemcpy( d_isotope_list,          isotope_list,          xs_length_numbers[0]*sizeof(unsigned), 			  cudaMemcpyHostToDevice );
	cudaMemcpy( d_number_density_matrix, number_density_matrix, n_materials*xs_length_numbers[0]*sizeof(float),    		  cudaMemcpyHostToDevice );
	cudaMemcpy( d_xs_data_Q,	     xs_data_Q, 	    (xs_length_numbers[2]+xs_length_numbers[0])*sizeof(float),    cudaMemcpyHostToDevice );
	std::cout << "Done.\n";
	// copy device pointer array to device array
	std::cout << "  Linked pointer arrays... ";
	cudaMemcpy( d_xs_data_scatter, 	xs_data_scatter_host,	MT_rows*MT_columns*sizeof(float*), cudaMemcpyHostToDevice); 	
	cudaMemcpy( d_xs_data_energy, 	xs_data_energy_host,	MT_rows*MT_columns*sizeof(float*), cudaMemcpyHostToDevice); 	
	std::cout << "Done.\n";
	// copy scattering data to device array pointers
	//std::cout << "  Scattering data... ";
	//for (int j=0 ; j < MT_columns ; j++){  //start after the total xs and total abs vectors
	//	for (int k=0 ; k < MT_rows ; k++){
	//		// scatter
	//		this_pointer =   xs_data_scatter     [k*MT_columns + j];
	//		cuda_pointer =   xs_data_scatter_host[k*MT_columns + j];
	//		printf("cpu=%p gpu=%p\n",this_pointer,cuda_pointer);
	//		if (this_pointer != NULL & k<MT_rows-1) {
	//			memcpy(&vlen,     &this_pointer[2],1*sizeof(int));
	//			memcpy(&next_vlen,&this_pointer[3],1*sizeof(int));
	//			while(xs_data_scatter[(k+1)*MT_columns + j ]==this_pointer){
	//				k++; //push k to the end of the copies so don't try to free it twice
	//			}
	//			cudaMemcpy(cuda_pointer,this_pointer,(2*vlen+2*next_vlen+4)*sizeof(float),cudaMemcpyHostToDevice);
	//		}
	//	}
	//}
	//std::cout << " Done.\n";
	// zero out tally arrays
	std::cout << "  Zeroing tally arrays... ";
	cudaMemcpy( d_tally_score, 	d_zeros,	n_tally*sizeof(float),    cudaMemcpyDeviceToDevice); 	
	cudaMemcpy( d_tally_square, d_zeros,	n_tally*sizeof(float),    cudaMemcpyDeviceToDevice); 	
	cudaMemcpy( d_tally_count,	d_zeros,	n_tally*sizeof(unsigned), cudaMemcpyDeviceToDevice); 	
	std::cout << "Done.\n";


}
void whistory::load_cross_sections(){
	
	std::cout << "\e[1;32m" << "Loading cross sections and unionizing..." << "\e[m \n";

	// set the string, make ints list
	std::istringstream ss(xs_isotope_string);
	std::string token;
	unsigned utoken;
	unsigned bytes,rows,columns;

	while(std::getline(ss, token, ',')) {
		utoken = std::atoi(token.c_str());
    		xs_isotope_ints.push_back(utoken);
	}

	// get data from python
	/* 	 need to do
	xs = unionize.cross_section_data()
	xs._init_from_string(this_string)
	xs._read_tables()
	xs._unionize()
	xs._insert_reactions()
	xs._allocate_arrays()
	xs._interpolate() */

	// variables
	PyObject *pName, *pModule, *pDict, *pFunc;
    PyObject *pArgs, *pValue, *pString, *pBuffObj, *pObjList;
    PyObject *call_result;
    PyObject *call_string,*arg_string;
    PyObject *xsdat_instance;
    PyObject *pClass;
    Py_buffer pBuff;
    int i, do_final;

    if (Py_IsInitialized()){
    	printf("Python interpreter already initialized\n");
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
        	return;	
    	}

    	pName = PyString_FromString("cross_section_data");
    	xsdat_instance = PyObject_CallMethodObjArgs(pModule,pName,NULL);
    	PyErr_Print();
    	Py_DECREF(pName);


    	if (xsdat_instance != NULL) {

		// init the libraries wanted
		char tope_string_c[256];
		call_string = PyString_FromString("_init_from_string");
		arg_string  = PyString_FromString(xs_isotope_string.c_str());
		call_result = PyObject_CallMethodObjArgs(xsdat_instance, call_string, arg_string, NULL);
		PyErr_Print();
		Py_DECREF(arg_string);
		Py_DECREF(call_string);
		Py_DECREF(call_result);
	
		// read the tables
		call_string = PyString_FromString("_read_tables");
		call_result = PyObject_CallMethodObjArgs(xsdat_instance, call_string, NULL);
		PyErr_Print();
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
        	return;
    	}

    	// get the MT array buffer
    	call_string = PyString_FromString("_get_MT_array_pointer");
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

    //
    // get and copy the unionized MT array
    //
	MT_rows    = pBuff.shape[0];
	MT_columns = pBuff.shape[1];
	bytes   = pBuff.len;
	//std::cout << "unionized MT array " << MT_rows << " " << MT_columns << " " << bytes << "\n";
    // allocate xs_data pointer arrays
    xs_data_MT       = new float  [MT_rows*MT_columns];
    // check to make sure bytes *= elements
    assert(bytes==MT_rows*MT_columns*4);
    // copy python buffer contents to pointer
    memcpy( xs_data_MT,   pBuff.buf , bytes );
    // cudaallocate device memory now that we know the size!
    cudaMalloc(&d_xs_data_MT,bytes);
    // release python variable to free memory
	Py_DECREF(call_result);

    // get the unionized main energy grid buffer
    call_string = PyString_FromString("_get_main_Egrid_pointer");
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

    //
    // get and copy unionized main energy grid
    //
	rows    = pBuff.shape[0];
	columns = pBuff.shape[1];
	bytes   = pBuff.len;
	//std::cout << "main e grid " << rows << " " << columns << " " << bytes << "\n";
    	// allocate xs_data pointer arrays
    xs_data_main_E_grid  = new float  [rows];
    // check to make sure bytes *= elements
    assert(bytes==rows*4);
    // copy python buffer contents to pointer
    memcpy( xs_data_main_E_grid,   pBuff.buf , bytes );
    // cudaallocate device memory now that we know the size!
    cudaMalloc(&d_xs_data_main_E_grid,bytes);
    // release python variable to free memory
    Py_DECREF(call_result);

    // mt number vector
    call_string = PyString_FromString("_get_MT_numbers_pointer");
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

    //
    // get and copy mt number vector
    //
	rows    = pBuff.shape[0];
	columns = pBuff.shape[1];
	bytes   = pBuff.len;
	//std::cout << "mt nums " << rows << " " << columns << " " << bytes << "\n";
    // allocate xs_data pointer arrays
    xs_MT_numbers      = new unsigned  [rows];
    // check to make sure bytes *= elements
    assert(bytes==rows*4);
    // copy python buffer contents to pointer
    memcpy( xs_MT_numbers,   pBuff.buf , bytes );
    // cudaallocate device memory now that we know the size!
    cudaMalloc(&d_xs_MT_numbers,bytes);
    // release python variable to free memory
    Py_DECREF(call_result);

    // mt number total vector
    call_string = PyString_FromString("_get_MT_numbers_total_pointer");
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

    //
    // get and copy unionized totals vector
    //
	rows    = pBuff.shape[0];
	columns = pBuff.shape[1];
	bytes   = pBuff.len;
	//std::cout << "totals " << rows << " " << columns << " " << bytes << "\n";
    // allocate xs_data pointer arrays
    xs_MT_numbers_total      = new unsigned  [rows];
    // check to make sure bytes *= elements
    assert(bytes==rows*4);
    // copy python buffer contents to pointer
    memcpy( xs_MT_numbers_total,   pBuff.buf , bytes );
    // cudaallocate device memory now that we know the size!
    cudaMalloc(&d_xs_MT_numbers_total,bytes);
    // release python variable to free memory
    Py_DECREF(call_result);

    // lengths vector
    call_string = PyString_FromString("_get_length_numbers_pointer");
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

	//
	// get and copy lengths vector
	//
	rows    = pBuff.shape[0];
	columns = pBuff.shape[1];
	bytes   = pBuff.len;
	//std::cout << "lengths " << rows << " " << columns << " " << bytes << "\n";
	// allocate xs_data pointer arrays
	xs_length_numbers     = new unsigned  [rows];
	// check to make sure bytes *= elements
	assert(bytes==rows*4);
	// copy python buffer contents to pointer
	memcpy( xs_length_numbers,   pBuff.buf , bytes );
	// cudaallocate device memory now that we know the size!
	cudaMalloc(&d_xs_length_numbers,bytes);
	// release python variable to free memory
	Py_DECREF(call_result);

	// AWR vector
	call_string = PyString_FromString("_get_awr_pointer");
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

	//
	// get and copy AWR vector
	//
	rows    = pBuff.shape[0];
	columns = pBuff.shape[1];
	bytes   = pBuff.len;
	//std::cout << "lengths " << rows << " " << columns << " " << bytes << "\n";
	// allocate xs_data pointer arrays
	awr_list     = new float  [rows];
	// check to make sure bytes *= elements
	assert(bytes==rows*4);
	// copy python buffer contents to pointer
	memcpy( awr_list,   pBuff.buf , bytes );
	// cudaallocate device memory now that we know the size!
	cudaMalloc(&d_awr_list,bytes);
	// release python variable to free memory
	Py_DECREF(call_result);

	// Q vector
	call_string = PyString_FromString("_get_Q_pointer");
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

	//
	// get and copy Q vector
	//
	rows    = pBuff.shape[0];
	columns = pBuff.shape[1];
	bytes   = pBuff.len;
	//std::cout << "lengths " << rows << " " << columns << " " << bytes << "\n";
	// allocate xs_data pointer arrays
	xs_data_Q     = new float  [rows];
	// check to make sure bytes *= elements
	assert(bytes==rows*4);
	// copy python buffer contents to pointer
	memcpy( xs_data_Q,   pBuff.buf , bytes );
	// cudaallocate device memory now that we know the size!
	cudaMalloc(&d_xs_data_Q,bytes);
	// release python variable to free memory
	Py_DECREF(call_result);

	std::cout << "\e[1;32m" << "Loading/copying scattering data and linking..." << "\e[m \n";

	////////////////////////////////////
	// do scattering stuff
	////////////////////////////////////

    float * temp = new float [128];
	for(int g=0;g<128;g++){temp[g]=123456789;}
	unsigned vlen;

    //ALLOCATE THE ARRAYS.
    xs_data_scatter      = new float* [MT_rows*MT_columns];
    xs_data_energy       = new float* [MT_rows*MT_columns];
    xs_data_scatter_host = new float* [MT_rows*MT_columns];
    xs_data_energy_host  = new float* [MT_rows*MT_columns];
    cudaMalloc(&d_xs_data_scatter,MT_rows*MT_columns*sizeof(float*));
    cudaMalloc(&d_xs_data_energy, MT_rows*MT_columns*sizeof(float*));
    // nu container
    float * nuBuff 	     = new float [MT_rows];
    // python variables for arguments
    PyObject 	*E_obj, *row_obj, *col_obj;
    PyObject 	*cdf_vector_obj, *mu_vector_obj , *vector_length_obj, *nextDex_obj, *nextE_obj, *lastE_obj; 
    PyObject 	*next_cdf_vector_obj, *next_mu_vector_obj , *next_vector_length_obj;
    PyObject 	*obj_list;
    Py_buffer 	muBuff, cdfBuff, next_muBuff, next_cdfBuff;
    float 		*this_pointer,*cuda_pointer;
    float  		nextE, lastE;
    float       this_energy;
    float 		nu_test;
    unsigned	this_MT, this_tope, vector_length_L;
    int 		vector_length,next_vector_length;
    int 		minusone = -1;
    unsigned 	muRows,  muColumns,  muBytes;
    unsigned 	cdfRows, cdfColumns, cdfBytes;
    unsigned 	next_muRows,  next_muColumns,  next_muBytes;
    unsigned 	next_cdfRows, next_cdfColumns, next_cdfBytes;
    unsigned 	nextDex;

    //set total cross sections to NULL
    for (int j=0 ; j<1*xs_length_numbers[0] ; j++){  //start after the total xs vectors
    		for (int k=0 ; k<MT_rows ; k++){
    			xs_data_scatter     [k*MT_columns + j] = 0;//NULL;
			xs_data_scatter_host[k*MT_columns + j] = 0;//NULL;
		}
	}

    	// do the rest of the MT numbers
    	for (int j=1*xs_length_numbers[0] ; j<MT_columns ; j++){  //start after the total xs vectors
    		for (int k=0 ; k<MT_rows ; k++){
    			// call cross_section_data instance to get buffer
    			row_obj     = PyInt_FromLong     (k);
    			col_obj     = PyInt_FromLong     (j);
    			call_string = PyString_FromString("_get_scattering_data");
			obj_list    = PyObject_CallMethodObjArgs(xsdat_instance, call_string, row_obj, col_obj, NULL);
			PyErr_Print();

			// get objects in the returned list  [nextDex,next_E,vlen,nextvlen,mu,cdf,nextmu,nextcdf]
			nextDex_obj  			= PyList_GetItem(obj_list,0);
			lastE_obj  				= PyList_GetItem(obj_list,1);
			nextE_obj 				= PyList_GetItem(obj_list,2);
			vector_length_obj 		= PyList_GetItem(obj_list,3);
			next_vector_length_obj 	= PyList_GetItem(obj_list,4);
			mu_vector_obj 			= PyList_GetItem(obj_list,5);
			cdf_vector_obj 			= PyList_GetItem(obj_list,6);
			next_mu_vector_obj 		= PyList_GetItem(obj_list,7);
			next_cdf_vector_obj 	= PyList_GetItem(obj_list,8);
			PyErr_Print();

			// expand list to c variables
			nextDex 	  		= PyInt_AsLong 	  (nextDex_obj);
			lastE 		  		= PyFloat_AsDouble(lastE_obj);
			nextE 		  		= PyFloat_AsDouble(nextE_obj);
			vector_length 		= PyInt_AsLong    (vector_length_obj);
			next_vector_length 	= PyInt_AsLong    (next_vector_length_obj);
			PyErr_Print();

			// set this pointer
			if(vector_length==0){   // EMPTY
				//printf("(%u,%u) empty\n",k,j);
				xs_data_scatter     [k*MT_columns + j] = NULL;
				xs_data_scatter_host[k*MT_columns + j] = NULL;
				PyErr_Print();
			}
			else if (vector_length==-1){  // NU!!!!!!!
				//printf("(%u,%u) nu\n",k,j);
				// this nu, not scatter, copy the entire column.
				// get data buffer from numpy array
				if (PyObject_CheckBuffer(mu_vector_obj) & PyObject_CheckBuffer(cdf_vector_obj)){
					PyObject_GetBuffer( mu_vector_obj,  &muBuff, PyBUF_ND);
					PyErr_Print();
				}
				else{
					PyErr_Print();
				    	fprintf(stderr, "Returned object does not support buffer interface\n");
				    	return;
				}
			
				// shape info
				muRows     =  muBuff.shape[0];
				muColumns  =  muBuff.shape[1];
				muBytes    =  muBuff.len;
			
				//make sure every is ok
				assert(muRows==MT_rows);
			
				// copy to regular array
				memcpy(nuBuff, muBuff.buf , MT_rows*sizeof(float));
			
				//write the returned values into the array, byte-copy into 64bit pointer
				for(k;k<MT_rows;k++){
					//std::cout << "copying nu value "<< nuBuff[k] <<" at energy "<< xs_data_main_E_grid[k]<< " MeV\n";
					memcpy( &xs_data_scatter_host[ k*MT_columns + j ] , &nuBuff[k] , 1*sizeof(float) );
					memcpy( &nu_test , &xs_data_scatter_host[ k*MT_columns + j ], 1*sizeof(float) );
					//std::cout << "copying nu value "<< nu_test <<" at energy "<< xs_data_main_E_grid[k]<< " MeV\n";
				}
			
				PyErr_Print();
			
				break;
			
			}
			else{
				// get data buffer from numpy array
				if (PyObject_CheckBuffer(mu_vector_obj) & PyObject_CheckBuffer(cdf_vector_obj) & PyObject_CheckBuffer(next_mu_vector_obj) & PyObject_CheckBuffer(next_cdf_vector_obj) ){
					PyObject_GetBuffer(      mu_vector_obj,       &muBuff, PyBUF_ND);
					PyObject_GetBuffer(     cdf_vector_obj,      &cdfBuff, PyBUF_ND);
					PyObject_GetBuffer( next_mu_vector_obj,  &next_muBuff, PyBUF_ND);
					PyObject_GetBuffer(next_cdf_vector_obj, &next_cdfBuff, PyBUF_ND);
					PyErr_Print();
				}
				else{
					PyErr_Print();
    			    		fprintf(stderr, "Returned object does not support buffer interface\n");
    			    		return;
				}
	
				// shape info
				muRows     =  muBuff.shape[0];
				muColumns  =  muBuff.shape[1];
				muBytes    =  muBuff.len;
				cdfRows    = cdfBuff.shape[0];
				cdfColumns = cdfBuff.shape[1];
				cdfBytes   = cdfBuff.len;
				next_muRows     = next_muBuff.shape[0];
				next_muColumns  = next_muBuff.shape[1];
				next_muBytes    = next_muBuff.len;
				next_cdfRows    = next_cdfBuff.shape[0];
				next_cdfColumns = next_cdfBuff.shape[1];
				next_cdfBytes   = next_cdfBuff.len;
	
				//make sure every is ok
				assert(muRows==   cdfRows);
				assert(muColumns==cdfColumns);
				assert(muBytes==  cdfBytes);
				assert(next_muRows==   next_cdfRows);
				assert(next_muColumns==next_cdfColumns);
				assert(next_muBytes==  next_cdfBytes);
	
				//allocate pointer, write into array
				//for cuda too
				this_pointer = new float [muRows+cdfRows+next_muRows+next_cdfRows+4];
				cudaMalloc(&cuda_pointer,(muRows+cdfRows+next_muRows+next_cdfRows+4)*sizeof(float));
				total_bytes_scatter += (muRows+cdfRows+next_muRows+next_cdfRows  +4)*sizeof(float);  // add to total count
				xs_data_scatter     [k*MT_columns + j] = this_pointer;
				xs_data_scatter_host[k*MT_columns + j] = cuda_pointer;
	
				//copy data from python buffer to pointer in array
				memcpy(&this_pointer[0], 			&lastE,   	  sizeof(float));
				memcpy(&this_pointer[1], 			&nextE,   	  sizeof(float));    // nextE   to first position
				memcpy(&this_pointer[2], 			&muRows,   	  sizeof(unsigned)); // len to third position
				memcpy(&this_pointer[3], 			&next_muRows, 	  sizeof(unsigned));
				memcpy(&this_pointer[4],			muBuff.buf,  	  muRows*sizeof(float));     // mu  to len bytes after
				memcpy(&this_pointer[4+  muRows],		cdfBuff.buf, 	  cdfRows*sizeof(float));     // cdf to len bytes after that
				memcpy(&this_pointer[4+ 2*muRows],		next_muBuff.buf,  next_muRows*sizeof(float));
				memcpy(&this_pointer[4+ 2*muRows+next_muRows],	next_cdfBuff.buf, next_cdfRows*sizeof(float));
				PyErr_Print();

				cudaMemcpy(cuda_pointer,this_pointer,(muRows+cdfRows+next_muRows+next_cdfRows+4)*sizeof(float),cudaMemcpyHostToDevice);

			}

			// replicate this pointer into array until nextE
			// for cuda too
			if (k < (MT_rows-1) ){
				while(k<nextDex-1){
					xs_data_scatter     [(k+1)*MT_columns + j] = this_pointer;
					xs_data_scatter_host[(k+1)*MT_columns + j] = cuda_pointer;
					k++;
				}
			}
	
			this_pointer = NULL;
			cuda_pointer = NULL;

		}
	}


	std::cout << "\e[1;32m" << "Loading/copying energy distribution data and linking..." << "\e[m \n";


    	////////////////////////////////////
    	// do energy stuff
    	////////////////////////////////////
	PyObject 	*law_obj,*MT_obj,*tope_obj,*pdf_vector_obj,*next_pdf_vector_obj;
	Py_buffer 	pdfBuff, next_pdfBuff;
	unsigned 	law=0;

     	//set total cross sections to NULL
    	for (int j=0 ; j<1*xs_length_numbers[0] ; j++){  //start after the total xs vectors
    		for (int k=0 ; k<MT_rows ; k++){
    			xs_data_energy     [k*MT_columns + j] = NULL;
			xs_data_energy_host[k*MT_columns + j] = NULL;
		}
	}

    	// do the rest of the MT numbers
    	for (int j=1*xs_length_numbers[0] ; j<MT_columns ; j++){  //start after the total xs vectors
    		//std::cout << "  at energy column " << j+1 << " of " << MT_columns<< "\n";
    		for (int k=0 ; k<MT_rows ; k++){

    			// call cross_section_data instance to get buffer
    			row_obj     = PyInt_FromLong     (k);
    			col_obj     = PyInt_FromLong     (j);
    			call_string = PyString_FromString("_get_energy_data");
			obj_list    = PyObject_CallMethodObjArgs(xsdat_instance, call_string, row_obj, col_obj, NULL);
			PyErr_Print();

			// get objects in the returned list
			nextDex_obj  		= PyList_GetItem(obj_list,0);
			lastE_obj  		= PyList_GetItem(obj_list,1);
			nextE_obj 		= PyList_GetItem(obj_list,2);
			vector_length_obj 	= PyList_GetItem(obj_list,3);
			next_vector_length_obj 	= PyList_GetItem(obj_list,4);
			law_obj			= PyList_GetItem(obj_list,5);
			mu_vector_obj 		= PyList_GetItem(obj_list,6);
			cdf_vector_obj 		= PyList_GetItem(obj_list,7);
			pdf_vector_obj 		= PyList_GetItem(obj_list,8);
			next_mu_vector_obj 	= PyList_GetItem(obj_list,9);
			next_cdf_vector_obj 	= PyList_GetItem(obj_list,10);
			next_pdf_vector_obj 	= PyList_GetItem(obj_list,11);
			PyErr_Print();

			// expand list to c variables
			nextDex 	  	= PyInt_AsLong 	  (nextDex_obj);
			lastE 		  	= PyFloat_AsDouble(lastE_obj);
			nextE 		  	= PyFloat_AsDouble(nextE_obj);
			vector_length 		= PyInt_AsLong    (vector_length_obj);
			next_vector_length 	= PyInt_AsLong    (next_vector_length_obj);
			law 			= PyInt_AsLong    (law_obj);

			// set this pointer
			if(vector_length==0){
				xs_data_energy     [k*MT_columns + j] = NULL;
				xs_data_energy_host[k*MT_columns + j] = NULL;
				PyErr_Print();
			}
			else{
				// get data buffer from numpy array
				if (PyObject_CheckBuffer(mu_vector_obj) & PyObject_CheckBuffer(cdf_vector_obj)){
					PyObject_GetBuffer( mu_vector_obj,  &muBuff, PyBUF_ND);
					PyObject_GetBuffer(cdf_vector_obj, &cdfBuff, PyBUF_ND);
					PyObject_GetBuffer(pdf_vector_obj, &pdfBuff, PyBUF_ND);
					PyObject_GetBuffer( next_mu_vector_obj,    &next_muBuff, PyBUF_ND);
					PyObject_GetBuffer( next_cdf_vector_obj,  &next_cdfBuff, PyBUF_ND);
					PyObject_GetBuffer( next_pdf_vector_obj,  &next_pdfBuff, PyBUF_ND);
					PyErr_Print();
				}
				else{
					PyErr_Print();
    			    		fprintf(stderr, "Returned object does not support buffer interface\n");
    			    		return;
				}
	
				// shape info
				muRows     =  muBuff.shape[0];
				muColumns  =  muBuff.shape[1];
				muBytes    =  muBuff.len;
				cdfRows    = cdfBuff.shape[0];
				cdfColumns = cdfBuff.shape[1];
				cdfBytes   = cdfBuff.len;
				next_muRows     = next_muBuff.shape[0];
				next_muColumns  = next_muBuff.shape[1];
				next_muBytes    = next_muBuff.len;
				next_cdfRows    = next_cdfBuff.shape[0];
				next_cdfColumns = next_cdfBuff.shape[1];
				next_cdfBytes   = next_cdfBuff.len;

				//std::cout << "C++ vlen/next " << muRows << " " << next_muRows << "\n";
	
				//make sure every is ok
				assert(muRows==   cdfRows);
				assert(muColumns==cdfColumns);
				assert(muBytes==  cdfBytes);
				assert(next_muRows==   next_cdfRows);
				assert(next_muColumns==next_cdfColumns);
				assert(next_muBytes==  next_cdfBytes);
	
				//allocate pointer, write into array
				//for cuda too
				this_pointer = new float [muRows + 2*cdfRows + next_muRows + 2*next_cdfRows + 5];
				cudaMalloc(&cuda_pointer,(muRows + 2*cdfRows + next_muRows + 2*next_cdfRows + 5)*sizeof(float));
				total_bytes_energy +=    (muRows + 2*cdfRows + next_muRows + 2*next_cdfRows + 5)*sizeof(float);    // add to total count
				xs_data_energy     [k*MT_columns + j] = this_pointer;
				xs_data_energy_host[k*MT_columns + j] = cuda_pointer;
	
				//copy data from python buffer to pointer in array
				memcpy(&this_pointer[0], 				&lastE,   	  sizeof(float));
				memcpy(&this_pointer[1], 				&nextE,   	  sizeof(float));    // nextE   to first position
				memcpy(&this_pointer[2], 				&muRows,   	  sizeof(unsigned)); // len to third position
				memcpy(&this_pointer[3], 				&next_muRows, 	  sizeof(unsigned));
				memcpy(&this_pointer[4], 				&law, 		  sizeof(unsigned));
				memcpy(&this_pointer[5],				muBuff.buf,  	  muRows*sizeof(float));     // mu  to len bytes after
				memcpy(&this_pointer[5+   muRows],			cdfBuff.buf, 	  cdfRows*sizeof(float));
				memcpy(&this_pointer[5+ 2*muRows],			pdfBuff.buf, 	  cdfRows*sizeof(float));
				memcpy(&this_pointer[5+ 3*muRows],			next_muBuff.buf,  next_muRows*sizeof(float));
				memcpy(&this_pointer[5+ 3*muRows +   next_muRows],	next_cdfBuff.buf, next_cdfRows*sizeof(float));
				memcpy(&this_pointer[5+ 3*muRows + 2*next_muRows],	next_pdfBuff.buf, next_cdfRows*sizeof(float));

				cudaMemcpy(cuda_pointer,this_pointer,(muRows+2*cdfRows+next_muRows+2*next_cdfRows+5)*sizeof(float),cudaMemcpyHostToDevice);

				//for(unsigned n=0;n<(muRows+2*cdfRows+next_muRows+2*next_cdfRows+5);n++){
				//	printf("%6.4E\n",this_pointer[n]);
				//}


				PyErr_Print();

			}

			// replicate this pointer into array until nextE
			// for cuda too
			if (k < (MT_rows-1) ){
				while(k<nextDex-1){
					xs_data_energy     [(k+1)*MT_columns + j] = this_pointer;
					xs_data_energy_host[(k+1)*MT_columns + j] = cuda_pointer;
					k++;
				}
			}
	
			this_pointer = NULL;
			cuda_pointer = NULL;

		}
	}

	if(do_final){
    	Py_Finalize();
    }

	std::cout << "\e[1;32m" << "Making material table..." << "\e[m \n";

    	//pass awr pointer to geometry object, make the number density table, copy pointers back
    	problem_geom.awr_list = awr_list;
    	problem_geom.make_material_table();
    	problem_geom.get_material_table(&n_materials,&n_isotopes,&material_list,&isotope_list,&number_density_matrix);  

    	assert(n_isotopes == xs_length_numbers[0]);

    	//do cudamalloc for these arrays
    	cudaMalloc(&d_material_list , 			n_materials*sizeof(unsigned) );
    	cudaMalloc(&d_isotope_list , 			n_isotopes*sizeof(unsigned) );
    	cudaMalloc(&d_number_density_matrix , 	n_materials*n_isotopes*sizeof(unsigned) );

}
void whistory::print_xs_data(){  // 0=isotopes, 1=main E points, 2=total numer of reaction channels, 3=matrix E points, 4=angular cosine points, 5=outgoing energy points
	unsigned dsum = 0;
	printf("\e[1;32m%-6s\e[m \n","Cross section data info:");
	std::cout << "--- Bytes ---" << "\n";
	std::cout << "  xs_length_numbers:        " << 6*sizeof(unsigned) << "\n";			  dsum += (6*sizeof(unsigned) );
	std::cout << "  xs_MT_numbers_total:      " << xs_length_numbers[0]*sizeof(unsigned) 	<< "\n";  dsum += (xs_length_numbers[0]*sizeof(unsigned) );
	std::cout << "  xs_MT_numbers:            " << (xs_length_numbers[2]+xs_length_numbers[0])*sizeof(unsigned)<< "\n";
	dsum += (xs_length_numbers[2]*sizeof(unsigned) );
	std::cout << "  xs_data_main_E_grid:      " << xs_length_numbers[1]	*sizeof(float)	<< "\n";  dsum += (xs_length_numbers[1]*sizeof(float)	  );
	std::cout << "  xs_data_MT:               " << MT_rows*MT_columns*sizeof(float)		<< "\n";  dsum += (MT_rows*MT_columns)*sizeof(float);
	std::cout << "  xs_data_scatter_pointers: " << MT_rows*MT_columns*sizeof(float*)		<< "\n";  dsum += (MT_rows*MT_columns)*sizeof(float*);
	std::cout << "  xs_data_energy_pointers:  " << MT_rows*MT_columns*sizeof(float*)		<< "\n";  dsum += (MT_rows*MT_columns)*sizeof(float*);
	std::cout << "  scatter data:             " << total_bytes_scatter				<< "\n";  dsum += (total_bytes_scatter);
	std::cout << "  energy data:              " << total_bytes_energy				<< "\n";  dsum += (total_bytes_energy);
	std::cout << "  TOTAL:                    " << dsum << " bytes \n";
	std::cout << "  TOTAL:                    " << dsum/1048576 << " MB \n";
}
void whistory::write_xs_data(std::string filename){

	std::cout << "\e[1;32m" << "Writing xs_data to " << filename << "... ";

	std::string this_name;
	// write MT array
	this_name = filename + ".MTarray";
	FILE* xsfile = fopen(this_name.c_str(),"w");
	this_name = filename + ".scatterptr";
	FILE* scatterfile = fopen(this_name.c_str(),"w");
	this_name = filename + ".energyptr";
	FILE* energyfile = fopen(this_name.c_str(),"w");
	for (int k=0;k<MT_rows;k++){
		for(int j=0;j<MT_columns;j++){
			fprintf(xsfile,"% 10.8E ",xs_data_MT[k*MT_columns+j]);
			fprintf(scatterfile,"%p ",xs_data_scatter_host[k*MT_columns+j]);
			fprintf(energyfile,"%p ",xs_data_energy_host[k*MT_columns+j]);
		}
		fprintf(xsfile,"\n");
		fprintf(scatterfile,"\n");
		fprintf(energyfile,"\n");
	}
	fclose(xsfile);
	fclose(scatterfile);
	fclose(energyfile);

	// write unionized E grid
	this_name = filename + ".Egrid";
	xsfile = fopen(this_name.c_str(),"w");
	for (int k=0;k<MT_rows;k++){
		fprintf(xsfile,"%10.8E\n",xs_data_main_E_grid[k]);
	}
	fclose(xsfile);

	// write MT number array
	this_name = filename + ".MTnums";
	xsfile = fopen(this_name.c_str(),"w");
	for(int j=0;j<MT_columns;j++){
		fprintf(xsfile,"%u\n",xs_MT_numbers[j]);
	}
	fclose(xsfile);

	// write (hopefully) covnerged fission source
	cudaMemcpy(space, d_fissile_points, Ndataset*sizeof(source_point), cudaMemcpyDeviceToHost);
	cudaMemcpy(E,     d_fissile_energy, Ndataset*sizeof(float),        cudaMemcpyDeviceToHost);
	this_name = filename + ".fission_source";
	xsfile = fopen(this_name.c_str(),"w");
	for(int j=0;j<Ndataset;j++){
		fprintf(xsfile,"% 6.4E % 6.4E % 6.4E %6.4E\n",space[j].x,space[j].y,space[j].z,E[j]);
	}
	fclose(xsfile); 

	std::cout << "Done." << "\e[m \n";

}
void whistory::print_pointers(){
	std::cout << "\e[1;32m" << "Pointer Info:" << "\e[m \n";
	std::cout << "--- HOST ---" << "\n";
	std::cout << "  space:               " <<   space   << "\n";
	std::cout << "  E:                   " <<   E       << "\n";
	std::cout << "  Q:                   " <<   Q       << "\n";
	std::cout << "  rn_bank:             " <<   rn_bank << "\n";
	std::cout << "  cellnum:             " <<   cellnum << "\n";
	std::cout << "  matnum:              " <<   matnum  << "\n";
	std::cout << "  isonum:              " <<   isonum  << "\n";
	std::cout << "  rxn:                 " <<   rxn     << "\n";
	std::cout << "  done:                " <<   done    << "\n";
	std::cout << "  yield:               " <<   yield   << "\n";
	std::cout << "  xs_length_numbers:   " << xs_length_numbers   << "\n"; 
	std::cout << "  xs_MT_numbers_total: " << xs_MT_numbers_total << "\n";
	std::cout << "  xs_MT_numbers:       " << xs_MT_numbers       << "\n";
	std::cout << "  xs_data_MT:          " << xs_data_MT          << "\n";
	std::cout << "  xs_data_main_E_grid: " << xs_data_main_E_grid << "\n";
	std::cout << "--- DEVICE ---" << "\n";
	std::cout << "d_space:               " << d_space   << "\n";
	std::cout << "d_E:                   " << d_E       << "\n";
	std::cout << "d_Q:                   " << d_Q       << "\n";
	std::cout << "d_rn_bank:             " << d_rn_bank << "\n";
	std::cout << "d_cellnum:             " << d_cellnum << "\n";
	std::cout << "d_matnum:              " << d_matnum  << "\n";
	std::cout << "d_isonum:              " << d_isonum  << "\n";
	std::cout << "d_rxn:                 " << d_rxn     << "\n";
	std::cout << "d_done:                " << d_done    << "\n";
	std::cout << "d_yield:               " << d_yield   << "\n";
	std::cout << "d_xs_length_numbers:   " << d_xs_length_numbers   << "\n"; 
	std::cout << "d_xs_MT_numbers_total: " << d_xs_MT_numbers_total << "\n";
	std::cout << "d_xs_MT_numbers:       " << d_xs_MT_numbers       << "\n";
	std::cout << "d_xs_data_MT:          " << d_xs_data_MT          << "\n";
	std::cout << "d_xs_data_main_E_grid: " << d_xs_data_main_E_grid << "\n";
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

	problem_geom.print_materials_table();

}
void whistory::sample_fissile_points(){

	std::cout << "\e[1;32m" << "Sampling initial fissile starting points uniformly... " << "\e[m \n";

	// iterate
	unsigned current_index = 0;
	unsigned valid_N = 0;

	while (current_index < N){
		
		// advance RN bank
		update_RNG();

		// set uniformly random positions on GPU
		set_positions_rand ( NUM_THREADS, N , outer_cell_type, d_space , d_rn_bank, outer_cell_dims);
		
		//run OptiX to get cell number, set as a hash run for fissile, writes 1/0 to matnum, trace_type=4
		trace(3, N);

		// compact
		res = cudppCompact(compactplan, d_valid_result, (size_t*)d_valid_N , d_remap , d_matnum , N);
		if (res != CUDPP_SUCCESS){fprintf(stderr, "Error in compacting\n");exit(-1);}

		// copy back and add
		cudaMemcpy( &valid_N, d_valid_N, 1*sizeof(unsigned), cudaMemcpyDeviceToHost);
		current_index += valid_N;

		//copy in new values, keep track of index, copies positions and direction
		copy_points( NUM_THREADS, valid_N, d_valid_N, current_index, d_valid_result, d_fissile_points, d_space, d_fissile_energy, d_E); 
		
		// print how far along we are
		std::cout << (float)current_index/(float)N*100.0 <<" \% done\r";

		if((float)current_index/(float)N > 1){
			std::cout << "100.00 \% done     \n";
		} 
	}

	std::cout << "Copying to starting points...\n";

	cudaMemcpy(d_space,	d_fissile_points	,N*sizeof(source_point),	cudaMemcpyDeviceToDevice);
	cudaMemcpy(d_E,    	d_fissile_energy	,N*sizeof(float),		cudaMemcpyDeviceToDevice);
	cudaMemcpy(d_rxn,  	rxn 			,N*sizeof(unsigned),		cudaMemcpyHostToDevice);
	//cudaFree(d_fissile_points);

	std::cout << "Done.\n";

	//write starting positions to file
	cudaMemcpy(space,d_space,N*sizeof(source_point),cudaMemcpyDeviceToHost);
	FILE* positionsfile = fopen("starting_positions","w");
	for(int k=0;k<N;k++){
		fprintf(positionsfile,"% 10.8E % 10.8E % 10.8E % 10.8E % 10.8E % 10.8E\n",space[k].x,space[k].y,space[k].z,space[k].xhat,space[k].yhat,space[k].zhat);
	}
	fclose(positionsfile);

	// advance RN bank
	update_RNG();

}
void whistory::write_to_file(source_point* array_in , unsigned N , std::string filename, std::string opentype){

	FILE* f  = fopen(filename.c_str(),opentype.c_str());
	source_point * hostdata = new source_point [N];
	cudaMemcpy(hostdata,array_in,N*sizeof(source_point),cudaMemcpyDeviceToHost);

	for(unsigned k = 0;  k<N ;k++){
		fprintf(f,"% 6.4E % 6.4E % 6.4E\n",hostdata[k].x,hostdata[k].y,hostdata[k].z);
	}

	delete hostdata;
	fclose(f);

}
void whistory::write_to_file(source_point* array_in , float* array_in2, unsigned N , std::string filename, std::string opentype){

	FILE* f  = fopen(filename.c_str(),opentype.c_str());
	source_point * hostdata = new source_point [N];
	float * hostdata2 = new float [N];
	cudaMemcpy(hostdata,array_in,N*sizeof(source_point),cudaMemcpyDeviceToHost);
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

	// re-base the yield so keff is 1
	//printf("CUDA ERROR4, %s\n",cudaGetErrorString(cudaPeekAtLastError()));
	rebase_yield( NUM_THREADS, N,  keff_cycle, d_rn_bank, d_yield);

	//reduce to check, unecessary in real runs
	//float keff_new = reduce_yield();
	//std::cout << "real keff "<< keff_cycle <<", artificially rebased keff " << keff_new <<"\n";

	// prefix sum scan (non-inclusive) yields to see where to write
	//printf("CUDA ERROR4.1, %s\n",cudaGetErrorString(cudaPeekAtLastError()));
	res = cudppScan( scanplan_int, d_scanned,  d_yield,  Ndataset );
	if (res != CUDPP_SUCCESS){fprintf(stderr, "Error in scanning yield values\n");exit(-1);}

	// swap the key/values to sort the rxn vector by tid to align the rest of data with the right reactions, done so 18 doesn't have to be assumed.  Sorting is necessary for prefix sum to work (order-dependent)
	//printf("CUDA ERROR4.2, %s\n",cudaGetErrorString(cudaPeekAtLastError()));
	res = cudppRadixSort(radixplan, d_remap, d_rxn, N);  
	if (res != CUDPP_SUCCESS){fprintf(stderr, "Error in sorting reactions\n");exit(-1);}

	//write_histories(0);

	//pop them in!  should be the right size now due to keff rebasing  
	//printf("CUDA ERROR4.3, %s\n",cudaGetErrorString(cudaPeekAtLastError()));
	pop_source( NUM_THREADS, N, d_isonum, d_remap, d_scanned, d_remap, d_yield, d_done, d_index, d_rxn, d_space, d_E , d_rn_bank , d_xs_data_energy, d_xs_data_scatter, d_fissile_points, d_fissile_energy, d_awr_list);
	//printf("CUDA ERROR4.4, %s\n",cudaGetErrorString(cudaPeekAtLastError()));
 
 	// rest run arrays
	cudaMemcpy( d_space,	d_fissile_points,	N*sizeof(source_point),	cudaMemcpyDeviceToDevice );
	cudaMemcpy( d_E,		d_fissile_energy,	N*sizeof(unsigned),		cudaMemcpyDeviceToDevice );
	cudaMemcpy( d_done,		zeros,				N*sizeof(unsigned),		cudaMemcpyHostToDevice );
	cudaMemcpy( d_cellnum,	zeros,				N*sizeof(unsigned),		cudaMemcpyHostToDevice );
	cudaMemcpy( d_matnum,	zeros,				N*sizeof(unsigned),		cudaMemcpyHostToDevice );
	cudaMemcpy( d_isonum,	zeros,				N*sizeof(unsigned),		cudaMemcpyHostToDevice );
	cudaMemcpy( d_yield,	zeros,				N*sizeof(unsigned),		cudaMemcpyHostToDevice );
	cudaMemcpy( d_rxn,		zeros,				N*sizeof(unsigned),		cudaMemcpyHostToDevice );
	cudaMemcpy( d_remap,	remap,				N*sizeof(unsigned),		cudaMemcpyHostToDevice );
	cudaMemcpy( d_index,	zeros,				N*sizeof(unsigned),		cudaMemcpyHostToDevice );

	//printf("CUDA ERROR4.6, %s\n",cudaGetErrorString(cudaPeekAtLastError()));
	// update RNG seeds
	update_RNG();
	//printf("CUDA ERROR4.7, %s\n",cudaGetErrorString(cudaPeekAtLastError()));
	// sync, these H2D and D2D copies aren't strictly synchronous
	cudaDeviceSynchronize();
	//printf("CUDA ERROR4.8, %s\n",cudaGetErrorString(cudaPeekAtLastError()));

}
void whistory::reset_fixed(){

	// rest read-in run arrays (ie not ones that are written to in-between)
	cudaMemcpy( d_space,		space,		Ndataset*sizeof(source_point),	cudaMemcpyHostToDevice );
	cudaMemcpy( d_done,			done,		Ndataset*sizeof(unsigned),		cudaMemcpyHostToDevice );
	cudaMemcpy( d_cellnum,		cellnum,	Ndataset*sizeof(unsigned),		cudaMemcpyHostToDevice );
	cudaMemcpy( d_matnum,		matnum,		Ndataset*sizeof(unsigned),		cudaMemcpyHostToDevice );
	cudaMemcpy( d_isonum,		isonum,		Ndataset*sizeof(unsigned),		cudaMemcpyHostToDevice );
	cudaMemcpy( d_yield,		zeros,		Ndataset*sizeof(unsigned),		cudaMemcpyHostToDevice );
	cudaMemcpy( d_rxn,			zeros,		Ndataset*sizeof(unsigned),		cudaMemcpyHostToDevice );
	cudaMemcpy( d_active,		remap,		Ndataset*sizeof(unsigned),		cudaMemcpyHostToDevice );

	//set position, direction, energy
	sample_fixed_source( NUM_THREADS, N, RNUM_PER_THREAD, d_active, d_rn_bank, d_E, d_space);

	// update RNG seeds
	update_RNG();

	// sync, these H2D and D2D copies aren't strictly synchronous
	cudaDeviceSynchronize();

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

	std::cout << "\e[1;32m" << "--- Running in " << runtype << " source mode --- " << "\e[m \n";
	std::cout << "\e[1;32m" << "--- Skipping "<< n_skip << " cycles, Running "<< n_cycles << " ACTIVE CYCLES, "<< N << " histories each--- " << "\e[m \n";

	// make sure fissile_points file is cleared
	fiss_name=filename;
	fiss_name.append(".fission_points");
	FILE* ffile = fopen(fiss_name.c_str(),"w");
	fclose(ffile);

	// open file for run stats
	stats_name=filename;
	stats_name.append(".stats");
	FILE* statsfile = fopen(stats_name.c_str(),"w");

	while(iteration<n_cycles){

		//printf("CUDA ERROR0, %s\n",cudaGetErrorString(cudaPeekAtLastError()));
		//write source positions to file if converged
		if(converged){
			write_to_file(d_space,d_E,N,fiss_name,"a+");
		}

		//printf("CUDA ERROR1, %s\n",cudaGetErrorString(cudaPeekAtLastError()));

		// reset edges and active number
		Nrun=N;
		edges[0]  = 0; 
		edges[1]  = 0; 
		edges[2]  = 0;
		edges[3]  = 0;
		edges[4]  = 0;
		edges[5]  = 0;
		edges[6]  = 0; 
		edges[7]  = 0;
		edges[8]  = 0;
		edges[9]  = 0;
		edges[10] = 0;

		//printf("CUDA ERROR2, %s\n",cudaGetErrorString(cudaPeekAtLastError()));
		//print_data(stream[0], NUM_THREADS, N, d_space, d_E, d_cellnum, d_matnum, d_isonum, d_rxn, d_done, d_yield);

		while(Nrun>0){
			//printf("CUDA ERROR, %s\n",cudaGetErrorString(cudaPeekAtLastError()));

			//if(RUN_FLAG==0){
			//	Nrun=Ndataset;
			//}
			//else if(RUN_FLAG==1){
			//	Nrun=N;
			//}

			//record stats
			fprintf(statsfile,"%u %10.8E\n",Nrun,get_time());
			
			// find what material we are in and nearest surface distance
			trace(2, Nrun);

			//find the main E grid index
			//find_E_grid_index_quad( NUM_THREADS, N,  qnodes_depth,  qnodes_width, d_active, d_qnodes_root, d_E, d_index, d_done);
			find_E_grid_index( NUM_THREADS, Nrun, xs_length_numbers[1],  d_remap, d_xs_data_main_E_grid, d_E, d_index, d_rxn);

			// run macroscopic kernel to find interaction length, macro_t, and reaction isotope, move to interactino length, set resample flag, 
			macroscopic( NUM_THREADS, Nrun,  n_isotopes, n_materials, MT_columns, outer_cell, d_remap, d_space, d_isonum, d_cellnum, d_index, d_matnum, d_rxn, d_xs_data_main_E_grid, d_rn_bank, d_E, d_xs_data_MT , d_number_density_matrix, d_done);

			// run tally kernel to compute spectra
			if(converged){
				tally_spec( NUM_THREADS, Nrun, n_tally, tally_cell, d_remap, d_space, d_E, d_tally_score, d_tally_square, d_tally_count, d_done, d_cellnum, d_rxn);
			}

			// run microscopic kernel to fi`nd reaction type
			microscopic( NUM_THREADS, Nrun, n_isotopes, MT_columns, d_remap, d_isonum, d_index, d_xs_data_main_E_grid, d_rn_bank, d_E, d_xs_data_MT , d_xs_MT_numbers_total, d_xs_MT_numbers, d_xs_data_Q, d_rxn, d_Q, d_done);

			// remap threads to still active data
			remap_active(&Nrun, &escatter_N, &escatter_start, &iscatter_N, &iscatter_start, &cscatter_N, &cscatter_start, &fission_N, &fission_start);

			// concurrent calls to do escatter/iscatter/abs/fission
			cudaThreadSynchronize();cudaDeviceSynchronize();
			escatter( stream[0], NUM_THREADS,   escatter_N, escatter_start , d_remap, d_isonum, d_index, d_rn_bank, d_E, d_space, d_rxn, d_awr_list, d_done, d_xs_data_scatter);
			iscatter( stream[1], NUM_THREADS,   iscatter_N, iscatter_start , d_remap, d_isonum, d_index, d_rn_bank, d_E, d_space, d_rxn, d_awr_list, d_Q, d_done, d_xs_data_scatter, d_xs_data_energy);
			cscatter( stream[2], NUM_THREADS,1, cscatter_N, cscatter_start , d_remap, d_isonum, d_index, d_rn_bank, d_E, d_space, d_rxn, d_awr_list, d_Q, d_done, d_xs_data_scatter, d_xs_data_energy); // 1 is for transport run mode, as opposed to 'pop' mode
			fission ( stream[3], NUM_THREADS,   fission_N,  fission_start,   d_remap,  d_rxn ,  d_index, d_yield , d_rn_bank, d_done, d_xs_data_scatter);  
			cudaDeviceSynchronize();
			//printf("CUDA ERROR12, %s\n",cudaGetErrorString(cudaPeekAtLastError()));

			if(RUN_FLAG==0){  //fixed source
				// pop secondaries back in
				keff_cycle += reduce_yield();
				prep_secondaries();
				pop_secondaries( NUM_THREADS, Ndataset, RNUM_PER_THREAD, d_completed, d_scanned, d_yield, d_done, d_index, d_rxn, d_space, d_E , d_rn_bank , d_xs_data_energy);
				//if(reduce_yield()!=0.0){printf("pop_secondaries did not reset all yields!\n");}
			}

		}

		//reduce yield and reset cycle
		if(RUN_FLAG==0){
			keff_cycle = 1.0 - 1.0/(keff_cycle+1.0);   // based on: Ntotal = Nsource / (1-k) 
			reset_fixed();
			Nrun=Ndataset;
		}
		else if (RUN_FLAG==1){	
			accumulate_keff(converged, iteration, &keff, &keff_cycle);
			//printf("CUDA ERROR3, %s\n",cudaGetErrorString(cudaPeekAtLastError()));
			accumulate_tally();
			//printf("CUDA ERROR4, %s\n",cudaGetErrorString(cudaPeekAtLastError()));
			reset_cycle(keff_cycle);
			//printf("CUDA ERROR5, %s\n",cudaGetErrorString(cudaPeekAtLastError()));
			Nrun=N;
			//write_histories(1);
			//printf("CUDA ERROR6, %s\n",cudaGetErrorString(cudaPeekAtLastError()));
		}

		// update active iteration
		if (converged){
			iteration++;
		}

		// print whatever's clever
		if(converged){
			     if(RUN_FLAG==0){std::cout << "Cumulative keff/sc-mult = " << keff << " / " << 1.0/(1.0-keff) << ", ACTIVE cycle " << iteration << ", cycle keff/sc-mult = " << keff_cycle << " / " << 1.0/(1.0-keff_cycle) << "\n";}
			else if(RUN_FLAG==1){printf("Cumulative keff =  %8.6E +- %6.4E , ACTIVE cycle %4u, cycle keff = %8.6E\n",keff,keff_err,iteration,keff_cycle);}
		}
		else{
			printf("Converging fission source... skipped cycle %4u\n",iteration_total+1);
		}

		fprintf(statsfile,"---- iteration %u done ----\n",iteration);
		
		//write_histories(2);
		//printf("CUDA ERROR7, %s\n",cudaGetErrorString(cudaPeekAtLastError()));
		//std::cout << "press enter to continue...\n";
		//std::cin.ignore();

		// set convergence flag
		if( iteration_total == n_skip-1){ 
			converged=1;
			reduced_yields_total = 0;
		}

		// advance iteration number, reset cycle keff
		iteration_total++;
		keff_cycle = 0.0;

	}


	// print total transport running time
	runtime = get_time() - runtime;
	if(runtime>60.0){
		std::cout << "RUNTIME = " << runtime/60.0 << " minutes.\n";
	}
	else{
		std::cout << "RUNTIME = " << runtime << " seconds.\n";
	}

	write_results(runtime,keff,"w");
	fclose(statsfile);

}
void whistory::write_results(float runtime, float keff, std::string opentype){

	std::string resname = filename;
	resname.append(".results");
	FILE* f  = fopen(resname.c_str(),opentype.c_str());
	fprintf(f,"runtime = % 6.4E m    k_eff = % 6.4E     \ncycles total %u skipped %u scored %u \n %u source neutrons per cycle\n",runtime/60.0,keff,n_skip+n_cycles,n_skip,n_cycles,N);
	fclose(f);

}
void whistory::write_tally(unsigned tallynum){

	//tallynum is unused at this point
	float tally_err 	= 0;
	float tally_err_sq 	= 0;
	float this_square 	= 0;
	float this_mean 	= 0;
	float this_count 	= 0;
	float Emin 			= 1e-11;
	float Emax 			= 20.0;
	float edge 			= 0.0;

	// write tally values
	std::string tallyname = filename+".tally";
	FILE* tfile = fopen(tallyname.c_str(),"w");
	for (int k=0;k<n_tally;k++){
		this_count  	= (float) 	(N*n_cycles);//tally_count_total[k];
		this_mean 		= 			tally_score_total[k];
		this_square 	= 			tally_square_total[k];
		tally_err_sq 	= 			(1.0/((this_count - 1.0))) * ( (this_count*this_square)/(this_mean*this_mean) -  1.0 ) ;
		if(tally_err_sq>0.0){ 	
			tally_err = sqrtf(tally_err_sq);}
		else{					
			tally_err = 0.0;}
		fprintf(tfile,"%10.8E %10.8E %lu\n", this_mean/this_count, tally_err, tally_count_total[k]);
	}
	fclose(tfile);

	//write spacing, since there is n_tally+1 values, including it with the tally bins in a flat text messes with reading in as a matrix
	std::string binsname = filename+".tallybins";
	tfile = fopen(binsname.c_str(),"w");
	for (int k=0;k<=n_tally;k++){
		edge = Emin*powf((Emax/Emin), ((float)k)/ ((float)n_tally) );
		fprintf(tfile,"%10.8E\n",edge);
	}
	fclose(tfile);
}
void whistory::prep_secondaries(){

	// scan the yields to determine where the individual threads write into the done data
	res = cudppScan( scanplan_int, d_scanned,  d_yield,  Ndataset );
	if (res != CUDPP_SUCCESS){fprintf(stderr, "Error in scanning yield values\n");exit(-1);}

	// compact the done data to know where to write
	res = cudppCompact( compactplan, d_completed , (size_t*) d_num_completed , d_remap , d_done , Ndataset);
	if (res != CUDPP_SUCCESS){fprintf(stderr, "Error in compacting done values\n");exit(-1);}

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
unsigned whistory::map_active(){

	unsigned num_active=0;

	// flip done flag
	flip_done(NUM_THREADS, Ndataset, d_done);

	// remap to active
	res = cudppCompact( compactplan, d_active , (size_t*) d_num_active , d_remap , d_done , Ndataset);
	if (res != CUDPP_SUCCESS){fprintf(stderr, "Error in compacting done values\n");exit(-1);}
	cudaMemcpy(&num_active,d_num_active,1*sizeof(unsigned),cudaMemcpyDeviceToHost);

	// flip done flag back	
	flip_done(NUM_THREADS, Ndataset, d_done);

	return num_active;
}
void whistory::remap_active(unsigned* num_active, unsigned* escatter_N, unsigned* escatter_start, unsigned* iscatter_N, unsigned* iscatter_start, unsigned* cscatter_N, unsigned* cscatter_start, unsigned* fission_N, unsigned* fission_start){

	unsigned resamp_N = 0;
	unsigned resamp_start = 0;

	// sort key/value of rxn/tid
	//printf("N=%u\n",*num_active);
	if(*num_active>1){
		res = cudppRadixSort(radixplan, d_rxn, d_remap, *num_active );  //everything in 900s doesn't need to be sorted anymore
		if (res != CUDPP_SUCCESS){fprintf(stderr, "Error in sorting reactions\n");exit(-1);}
	}

	// launch edge detection kernel, writes mapped d_edges array
	reaction_edges(NUM_THREADS, *num_active, d_edges, d_rxn);

	// calculate lengths and starting indicies for the blocks, 0 indicates not found
	if (edges[0]){
		*escatter_N 	= (edges[1]-edges[0])+1;
		*escatter_start	= edges[0]-1;
	}
	else{
		*escatter_N 	= 0;
		*escatter_start	= 0;
	}
	if (edges[2]){
		*iscatter_N 	= (edges[3]-edges[2])+1;
		*iscatter_start	= edges[2]-1;
	}
	else{
		*iscatter_N 	= 0;
		*iscatter_start	= 0;
	}
	if (edges[4]){
		*cscatter_N 	= (edges[5]-edges[4])+1;
		*cscatter_start	= edges[4]-1;
	}
	else{
		*cscatter_N 	= 0;
		*cscatter_start	= 0;
	}
	if (edges[6]){
		   resamp_N 	= (edges[7]-edges[6])+1;
		   resamp_start	= edges[6]-1;
	}
	else{
		   resamp_N 	= 0;
		   resamp_start	= 0;
	}
	if (edges[8]){
		 *fission_N 	= (edges[9]-edges[8])+1;
		 *fission_start	= edges[8]-1;
	}
	else{
		 *fission_N 	= 0;
		 *fission_start	= 0;
	}

	//calculate total active, [10] is the starting index of >=811, so if==1, means that we are done!
	if(edges[8]==1){*num_active=0;}
	else           {
		*num_active=*escatter_N + *iscatter_N + *cscatter_N + resamp_N;
	}

	// debug
	//if(*num_active!=edges[8]-1){
	//	//print
	//	printf("num_active %u , edges[8] %u\n",*num_active,edges[8]-1);
	//	printf("nactive = %u, edges %u %u %u %u %u %u %u %u %u %u %u \n",*num_active,edges[0],edges[1],edges[2],edges[3],edges[4],edges[5],edges[6],edges[7],edges[8],edges[9],edges[10]);
	//	printf("escatter s %u n %u, iscatter s %u n %u, cscatter s %u n %u, resamp s %u n %u, fission s %u n %u \n\n",*escatter_start,*escatter_N,*iscatter_start,*iscatter_N,*cscatter_start,*cscatter_N,resamp_start,resamp_N, *fission_start, *fission_N);
	//	//dump
	//	write_to_file(d_remap, d_rxn, N,"remap","w");
	//	//exit
	//	assert(*num_active==edges[8]-1);
	//}

	// ensure order
	assert(*num_active<=edges[8]-1);
	if(*iscatter_N>0){ assert(*iscatter_start >= *escatter_start);}
	if(*cscatter_N>0){ assert(*cscatter_start >= *iscatter_start);}
	if(resamp_N>0){    assert(   resamp_start >= *cscatter_start);}

	// rezero edge vector (mapped, so is reflected on GPU)
	edges[0]  = 0; 
	edges[1]  = 0; 
	edges[2]  = 0;
	edges[3]  = 0;
	edges[4]  = 0;
	edges[5]  = 0;
	edges[6]  = 0;
	edges[7]  = 0;
	edges[8]  = 0;
	edges[9]  = 0;
	edges[10] = 0;

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

	tally_cell = cell;

}
void whistory::write_histories(unsigned iteration){

	//allocate
	unsigned*  done2;	
	unsigned*  cellnum2;
	unsigned*  matnum2;	
	unsigned*  isonum2;	
	unsigned*  yield2;
	unsigned*  rxn2;
	unsigned*  dex2;
	source_point* space2;
	float* E2;
	done2 		= new unsigned [N];
	cellnum2	= new unsigned [N];
	matnum2		= new unsigned [N];
	isonum2		= new unsigned [N];
	yield2 		= new unsigned [N];
	rxn2		= new unsigned [N];	
	dex2 		= new unsigned [N];
	space2 		= new source_point [N];
	E2 			= new float [N];

	// print actions
	std::string histfile_name = "history_file";
	char numstr[5];
	sprintf(numstr,"%u",iteration);
	histfile_name += numstr;
	printf("Writing to \"%s\"\n",histfile_name.c_str());

	// open file, appending to a new file.
	FILE* history_file = fopen(histfile_name.c_str(),"w");

	// write iteration delimiter
	fprintf(history_file,"==== ITERATION %u ====\n",iteration);

	// copy gemetrical (positions / directions)
	cudaMemcpy(space2,d_space,N*sizeof(source_point),cudaMemcpyDeviceToHost);

	// copy energies
	cudaMemcpy(E2,d_E,N*sizeof(float),cudaMemcpyDeviceToHost);

	// copy cell numbeer
	cudaMemcpy(cellnum2,d_cellnum,N*sizeof(unsigned),cudaMemcpyDeviceToHost);

	// copy material numbers
	cudaMemcpy(matnum2,d_matnum,N*sizeof(unsigned),cudaMemcpyDeviceToHost);

	// copy done flags
	cudaMemcpy(done2,d_done,N*sizeof(unsigned),cudaMemcpyDeviceToHost);

	// copy isotope numbers
	cudaMemcpy(isonum2,d_isonum,N*sizeof(unsigned),cudaMemcpyDeviceToHost);

	// copy yields
	cudaMemcpy(yield2,d_yield,N*sizeof(unsigned),cudaMemcpyDeviceToHost);

	// copy reaction numbers
	cudaMemcpy(rxn2,d_rxn,N*sizeof(unsigned),cudaMemcpyDeviceToHost);

	// copy index numbers
	cudaMemcpy(dex2,d_index,N*sizeof(unsigned),cudaMemcpyDeviceToHost);

	// sync device before write and return
	cudaDeviceSynchronize();

	// write history data to file
	for(unsigned k=0;k<N;k++){
		//fprintf(history_file,"tid %u (x,y,z) %8.6E %8.6E %8.6E (x,y,z)-hat %8.6E %8.6E %8.6E surf_dist %8.6E macro_t %8.6E enforce_BC %u E %8.6E cellnum %u matnum %u isonum %u rxn %u done %u yield %u\n",k,space2[k].x,space2[k].y,space2[k].z,space2[k].xhat,space2[k].yhat,space2[k].zhat,space2[k].surf_dist,space2[k].macro_t,space2[k].enforce_BC,E2[k],cellnum2[k],matnum2[k],isonum2[k],rxn2[k],done2[k],yield2[k]);
		fprintf(history_file,"a[%u,1:6]=[%8.6E,%8.6E,%8.6E,%8.6E,%u,%u,%u]\n",k,space2[k].x,space2[k].y,space2[k].z,E2[k],rxn2[k],yield2[k],dex2[k]);
	}
 	fclose(history_file);

 	//deallocate so can be alloaed again next time
 	delete done2 	;
	delete cellnum2	;
	delete matnum2	;
	delete isonum2	;
	delete yield2 	;
	delete rxn2		;
	delete space2 	;
	delete E2 		;
	delete dex2		;

}
void whistory::create_quad_tree(){

	std::cout << "\e[1;32m" << "Building quad tree for energy search... " << "\e[m \n";

	// node vectors
	std::vector<qnode_host>   nodes;
	std::vector<qnode_host>   nodes_next;

	// node variables
	qnode_host  this_qnode;
	qnode*      cuda_qnode;

	// build bottom-up from the unionized E vector
	unsigned k, depth;
	unsigned rows_by_four = MT_rows;
	for(k=0; k < rows_by_four ; k = k+4){
		this_qnode.node.values[0] = xs_data_main_E_grid[ k+0 ];
		this_qnode.node.values[1] = xs_data_main_E_grid[ k+1 ];
		this_qnode.node.values[2] = xs_data_main_E_grid[ k+2 ];
		this_qnode.node.values[3] = xs_data_main_E_grid[ k+3 ];
		this_qnode.node.leaves[0] = (qnode*) ((long unsigned)k + 0);  // recast grid index as the pointer for lowest nodes
		this_qnode.node.leaves[1] = (qnode*) ((long unsigned)k + 1);
		this_qnode.node.leaves[2] = (qnode*) ((long unsigned)k + 2);
		this_qnode.node.leaves[3] = (qnode*) ((long unsigned)k + 3);
		cudaMalloc(&cuda_qnode,sizeof(qnode));
		cudaMemcpy(cuda_qnode,&this_qnode.node,sizeof(qnode),cudaMemcpyHostToDevice);
		this_qnode.cuda_pointer = cuda_qnode;
		nodes.push_back(this_qnode);
	}
	//do the last values
	if(MT_rows%4){
		unsigned n=0;
		for(k=rows_by_four;k<MT_rows;k++){
			this_qnode.node.values[n] = xs_data_main_E_grid[ k ];
			this_qnode.node.leaves[n] = (qnode*) ((long unsigned)k); 
			n++;
		}
		//repeat the last values
		for(k=n;k<4;k++){
			this_qnode.node.values[k] = this_qnode.node.values[n-1];
			this_qnode.node.leaves[k] = this_qnode.node.leaves[n-1];
		}
		// device allocate and add to vector
		cudaMalloc(&cuda_qnode,sizeof(qnode));
		cudaMemcpy(cuda_qnode,&this_qnode.node,sizeof(qnode),cudaMemcpyHostToDevice);
		this_qnode.cuda_pointer = cuda_qnode;
		nodes.push_back(this_qnode);
	}


	//now build it up!  length will *always* need to be a multiple of 4.  this routine pads the end nodes with 
	unsigned this_width = nodes.size();
	unsigned lowest_length = this_width;
	unsigned mod4 		= this_width % 4;
	unsigned end_depth=  (logf(lowest_length)/logf(4))+1;
	unsigned num_repeats = 1;
	unsigned starting_index = 0;
	float inf = 1e45;
	//std::cout << "end_depth="<<end_depth<<"\n";
	for(unsigned depth=0;depth<end_depth;depth++){
		//for(unsigned copy_repeats=0;copy_repeats<num_repeats;copy_repeats++){
		//	starting_index=copy_repeats*this_width;
		//for(unsigned copy_iteration=0;copy_iteration<4;copy_iteration++){
			for( k=starting_index;k<(starting_index+this_width-mod4);k=k+4){
				//std::cout << "k=" << k << "\n";
				this_qnode.node.values[0] = nodes[ k+0 ].node.values[0]; // can use 0 since values overlap
				this_qnode.node.values[1] = nodes[ k+1 ].node.values[0];
				this_qnode.node.values[2] = nodes[ k+2 ].node.values[0];
				this_qnode.node.values[3] = nodes[ k+3 ].node.values[0];  
				this_qnode.node.leaves[0] = nodes[ k+0 ].cuda_pointer;  // set pointers as the cuda pointer of the children
				this_qnode.node.leaves[1] = nodes[ k+1 ].cuda_pointer;
				this_qnode.node.leaves[2] = nodes[ k+2 ].cuda_pointer;
				this_qnode.node.leaves[3] = nodes[ k+3 ].cuda_pointer;
				cudaMalloc(&cuda_qnode,sizeof(qnode));
				cudaMemcpy(cuda_qnode,&this_qnode.node,sizeof(qnode),cudaMemcpyHostToDevice);
				this_qnode.cuda_pointer = cuda_qnode;
				nodes_next.push_back(this_qnode);
			}
			if(mod4){
				//std::cout << "adding padded node at " << nodes_next.size()-1 << "\n";
				unsigned n=0;
				for( k=(starting_index+this_width-mod4) ; k<(starting_index+this_width) ; k++){
					//std::cout <<"n="<<n << " k="<<k<<"\n";
					this_qnode.node.values[n] = nodes[ k ].node.values[0];
					this_qnode.node.leaves[n] = nodes[ k ].cuda_pointer;
					n++;
				}
				for( n;n<4;n++){
					//std::cout <<"n="<<n << " k="<<k<<"\n";
					this_qnode.node.values[n] = inf;
					this_qnode.node.leaves[n] = NULL;
				}
				cudaMalloc(&cuda_qnode,sizeof(qnode));
				cudaMemcpy(cuda_qnode,&this_qnode.node,sizeof(qnode),cudaMemcpyHostToDevice);
				this_qnode.cuda_pointer = cuda_qnode;
				nodes_next.push_back(this_qnode);
			}
		//}
	//}
		if(mod4){
			this_width=(this_width)/4+1;
		}
		else{
			this_width=this_width/4;
		}
		mod4=this_width%4;
		nodes=nodes_next;
		nodes_next.clear();
		num_repeats=num_repeats*4;
		//std::cout << "--------------------------------------\n";
		//for(int g=0;g<nodes.size();g++){ //node vector check
		//	std::cout << "node " << g << " values " << nodes[g].node.values[0] << " " << nodes[g].node.values[1] << " "<< nodes[g].node.values[2] << " "<< nodes[g].node.values[3] << " "<< nodes[g].node.values[4] << " " <<"\n";
		//	std::cout << "node " << g << " leaves " << nodes[g].node.leaves[0] << " " << nodes[g].node.leaves[1] << " "<< nodes[g].node.leaves[2] << " "<< nodes[g].node.leaves[3] << " " <<"\n";
		//}
	}


	// copy size to object vars
	qnodes_depth = end_depth;
	qnodes_width = nodes.size();

	//copy root nodes vector to object variable
	//qnodes = new qnode[qnodes_width];
	//for(k=0;k<qnodes_width;k++){   //only need to copy heads, they have pointers to the rest in them
	//		qnodes[k].values[0] = nodes[k].node.values[0];
	//		qnodes[k].values[1] = nodes[k].node.values[1];
	//		qnodes[k].values[2] = nodes[k].node.values[2];
	//		qnodes[k].values[3] = nodes[k].node.values[3];
	//		qnodes[k].values[4] = nodes[k].node.values[4];
	//		qnodes[k].leaves[0] = nodes[k].node.leaves[0];
	//		qnodes[k].leaves[1] = nodes[k].node.leaves[1];
	//		qnodes[k].leaves[2] = nodes[k].node.leaves[2];
	//		qnodes[k].leaves[3] = nodes[k].node.leaves[3];
	//}

	// make and copy device data
	cudaMalloc(	&d_qnodes_root,				sizeof(qnode)	);
	cudaMemcpy(	 d_qnodes_root,	&nodes[0].node,		sizeof(qnode),	cudaMemcpyHostToDevice); 

	std::cout << "  Complete.  Depth of tree is "<< qnodes_depth << ", width is "<< qnodes_width <<".\n";

}
