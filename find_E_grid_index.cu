#include <cuda.h>
#include <stdio.h>
#include "datadef.h"

__global__ void find_E_grid_index_kernel(unsigned N, unsigned N_energies , unsigned* remap, float * main_E_grid, float* E , unsigned * index, unsigned* rxn){

	int tid_in = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid_in >= N){return;}

	// return if terminated
	unsigned this_rxn=rxn[tid_in];
	if (this_rxn>=900){return;}

	//remap
	int tid=remap[tid_in];

	// load data
	float value = E[tid];
	unsigned donesearching = 0;
	unsigned cnt  = 1;
	unsigned powtwo = 2;
	//unsigned olddex=index[tid];
	int dex  = (N_energies-1) / 2;  //N_energiesgth starts at 1, duh

	//printf("%p %d %10.4E\n",main_E_grid,dex,value);
	//int k;

	while(!donesearching){
		powtwo = powtwo * 2;
		if      ( 	main_E_grid[dex]   <= value && 
					main_E_grid[dex+1] >  value ) { donesearching = 1; }
		else if ( 	main_E_grid[dex]   >  value ) { dex  = dex - ((N_energies / powtwo) + 1) ; cnt++; }  // +1's are to do a ceiling instead of a floor on integer division
		else if ( 	main_E_grid[dex]   <  value ) { dex  = dex + ((N_energies / powtwo) + 1) ; cnt++; }

		if(cnt>30){
			donesearching=1;
			printf("binary search iteration overflow! %p %d % 10.8f tid=%u\n",main_E_grid,N_energies,value,tid);
			dex=0;
		}

		// edge checks... fix later???
		if(dex<0){
			//printf("binary search error! dex=%d, (ptr,N_energies,value) %p %d % 10.8f\n",dex,main_E_grid,N_energies,value);
			//for(k=0;k<N_energies;k++){printf("%10.8E\n",main_E_grid[k]);}
			dex=0;
			//donesearching=1;
		}
		if(dex>=N_energies){
			//printf("binary search error! dex=%d, (ptr,N_energies,value) %p %d % 10.8f\n",dex,main_E_grid,N_energies,value);
			//for(k=0;k<N_energies;k++){printf("%10.8E\n",main_E_grid[k]);}
			dex=N_energies-1;
			//donesearching=1;
		}
	}


	//write output index
	index[tid]=dex;
	//if(olddex!=dex){printf("E_i %6.4E E %6.4E E_i+1 %6.4E, dex %u quaddex %u\n",main_E_grid[dex],value,main_E_grid[dex+1],dex,olddex);}

}

void find_E_grid_index(unsigned NUM_THREADS, unsigned N, unsigned N_energies, unsigned* active, float * main_E_grid, float* E , unsigned * index , unsigned* rxn){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	find_E_grid_index_kernel <<< blks, NUM_THREADS >>> ( N, N_energies, active, main_E_grid,  E , index , rxn);
	cudaThreadSynchronize();

}

