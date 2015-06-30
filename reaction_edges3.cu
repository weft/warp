#include <cuda.h>
#include <stdio.h>

__global__ void reaction_edges_kernel( unsigned N, unsigned* edges, unsigned* rxn){

	int tid = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid >= N){return;}

	// the reaction vector has been sorted by this time.  data loads are expensive, computation is free, do as much as you can with 2 loads
	// need to find the lower AND upper bounds of the blocks.
	// edge array format is:
	// 0  = lower bound for 2 block 
	// 1  = upper bound for 2 block
	// 2  = lower bound for 51/90 block 
	// 3  = upper bound for 51/90 block
	// 4  = lower bound for 91 block 
	// 5  = upper bound for 91 block
	// 6  = lower bound for 800 block 
	// 7  = upper bound for 800 block
	// 8  = lower bound for 811/845 block 
	// 9  = upper bound for 811/845 block
	// 10 = lower bound for >900 block

	// init, load second if not last element
	unsigned rxn1, rxn2; 
	rxn1 = rxn[tid];
	if(tid<N-1){
		rxn2 = rxn[tid+1];
	}

	if(tid==0){
	// do first element
		if ( rxn1 == 2 ){
			edges[0] = 1;
		}
		else if( rxn1 >=51 && rxn1 <= 90 ){
			edges[2] = 1;
		}
		else if( rxn1 == 91 ){
			edges[4] = 1;
		} 
		else if( rxn1 == 800 | rxn1 == 801 ){
			edges[6] = 1;
		}
		else if( rxn1 >= 811 && rxn1 <=845 ){
			edges[8] = 1;
		}
		else if( rxn1 >= 900 ){
			edges[10] = 1;   // [10] = 1 is basically a done flag
		}
	}

	//check interior edges
	if (tid<N-1){
		// 2-edge
		if( rxn1 ==2 ){
			if( rxn2 >=51 && rxn2 <= 90){
				edges[1] = tid+1;
				edges[2] = tid+2;
			}
			else if( rxn2 == 91 ){
				edges[1] = tid+1;
				edges[4] = tid+2;
			}
			else if( rxn2 == 800 | rxn2==801 ){
				edges[1] = tid+1;
				edges[6] = tid+2;
			}
			else if( rxn2 >=811 && rxn2 <= 845 ){
				edges[1] = tid+1;
				edges[8] = tid+2;
			}
			else if( rxn2 >= 900 ){
				edges[1]  = tid+1;
				edges[10] = tid+2;
			}
		}

		// 51-90 edge
		if( rxn1 >=51 && rxn1 <= 90 ){
			if( rxn2 == 91 ){
				edges[3] = tid+1;
				edges[4] = tid+2;
			}
			else if( rxn2 == 800 | rxn2==801 ){
				edges[3] = tid+1;
				edges[6] = tid+2;
			}
			else if( rxn2 >=811 && rxn2 <= 845 ){
				edges[3] = tid+1;
				edges[8] = tid+2;
			}
			else if( rxn2 >= 900 ){
				edges[3]  = tid+1;
				edges[10] = tid+2;
			}
		}

		// 91 edge
		if( rxn1 == 91 ){
			if( rxn2 == 800 | rxn2==801 ){
				edges[5] = tid+1;
				edges[6] = tid+2;
			}
			else if( rxn2 >=811 && rxn2 <= 845 ){
				edges[5] = tid+1;
				edges[8] = tid+2;
			}
			else if( rxn2 >= 900 ){
				edges[5]  = tid+1;
				edges[10] = tid+2;
			}
		}

		// 800 edge
		if( rxn1 == 800 | rxn1==801 ){
			if( rxn2 >=811 && rxn2 <= 845 ){
				edges[7] = tid+1;
				edges[8] = tid+2;
			}
			else if( rxn2 >= 900 ){
				edges[7]  = tid+1;
				edges[10] = tid+2;
			}
		}

		// 811-845 edge
		if( rxn1 >=811 && rxn1 <= 845 ){
			if( rxn2 >= 900 ){
				edges[9]  = tid+1;
				edges[10] = tid+2;
			}
		}
	}
	else{

		// do last element
		if ( rxn1 == 2 ){
			edges[1] = N;
		}
		else if( rxn1 >=51 && rxn1 <= 90 ){
			edges[3] = N;
		}
		else if( rxn1 == 91 ){
			edges[5] = N;
		} 
		else if( rxn1 == 800 | rxn1==801 ){
			edges[7] = N;
		}
		else if( rxn1 >= 811 && rxn1 <=845 ){
			edges[9] = N;
		}
		//else if(rxn1 >= 900 ){
		//	no top edge for >=900
		//}
	}

}

void reaction_edges( unsigned NUM_THREADS,  unsigned N, unsigned* d_edges, unsigned* d_rxn){

	if(N<1){return;}
	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	reaction_edges_kernel <<< blks, NUM_THREADS >>> ( N, d_edges, d_rxn);
	cudaDeviceSynchronize();

}

