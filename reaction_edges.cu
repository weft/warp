#include <cuda.h>
#include <stdio.h>

__global__ void reaction_edges_kernel( unsigned N, unsigned* edges, unsigned* rxn){

	int tid = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid >= N){return;}

	// the reaction vector has been sorted by this time.  data loads are expensive, computation is free, do as much as you can with 2 loads
	// need to find the lower AND upper bounds of the blocks.
	// array structure is:
	// 0  = done flag
	// 1  = lower bound for 2 block 
	// 2  = upper bound for 2 block
	// 3  = lower bound for 51/90 block 
	// 4  = upper bound for 51/90 block
	// 5  = lower bound for 91 block 
	// 6  = upper bound for 91 block
	// 7  = lower bound for 800 block 
	// 8  = upper bound for 800 block
	// 9  = lower bound for 811/845 block 
	// 10 = upper bound for 811/845 block

	unsigned rxn1, rxn2; 
	rxn1 = rxn[tid];
	if(tid<N-1){
		rxn2 = rxn[tid+1];
	}
	else{  //check last element
		if(rxn1==2){
			edges[2] = tid+1;
		}
		else if(rxn1>=51 & rxn1<=90){
			edges[4] = tid+1;
		}
		else if(rxn1==91){
			edges[6] = tid+1;
		}
		else if(rxn1==800 & rxn1==800){
			edges[8] = tid+1;
		}
		else if(rxn1>=811 & rxn1<=845){
			edges[10] = tid+1;  //batch is done, set done flag
			//edges[0] = 1;
		}
		else if(rxn1>=900){
			//edges[0]=1;   //batch is definitely done, set done flag
		}
		return;
	}

	// first element doesn't have a preceeding
	if(tid==0){
		if(rxn1==2){
			edges[1] = tid;
		}
		else if(rxn1>=51 & rxn1<=90){
			edges[3] = tid;
		}
		else if(rxn1==91){
			edges[5] = tid;
		}
		else if(rxn1==800 & rxn1==800){
			edges[7] = tid;
		}
		else if(rxn1>=811 & rxn1<=845){
			edges[9] = tid;  //batch is done, set done flag
			edges[0] = 1;
		}
		else if(rxn1>=900){
			edges[10]=1;   //batch is definitely done, set done flag
		}
		return;
	}
	else{
		
		if(rxn1<2){   //can test start of 2   
			if(rxn2==2){
				edges[1]=tid+1;
				//return;
			}
		}
		else{         //can test the end of 2
			if(rxn1==2 & rxn2>2){
				edges[2]=tid+1;
				//return;
			}
		}

		if(rxn1<51){   //can test start of 51-90   
			if(rxn2>=51 & rxn2 <91){
				edges[3]=tid+1;
				//printf("found iscatter start %u\n",tid+1);
				//return;
			}
		}
		else{         //can test the end of 51-90
			if(rxn1>=51 & rxn1<=90 & rxn2>90){
				edges[4]=tid+1;
				//return;
			}
		}

		if(rxn1<91){   //can test start of 91  
			if(rxn2==91){
				edges[5]=tid+1;
				//return;
			}
		}
		else{         //can test the end of 91
			if(rxn1==91 & rxn2>91){
				edges[6]=tid+1;
				//return;
			}
		}

		if(rxn1<800){   //can test start of 800  
			if(rxn2==800){
				edges[7]=tid+1;
				//return;
			}
		}
		else{         //can test the end of 800
			if(rxn1==800 & rxn2>800){
				edges[8]=tid+1;
				//return;
			}
		}

		if(rxn1<811){   //can test start of 811-845  
			if(rxn2>=811 & rxn2 <=845){
				edges[9]=tid+1;
				//return;
			} 
		}
		else{         //can test the end of 811-845
			if(rxn1>=811 & rxn1<=845 & rxn2>845){
				edges[10]=tid+1;
				//return;
			}
		}

	}


}

void reaction_edges( unsigned NUM_THREADS,  unsigned N, unsigned* d_edges, unsigned* d_rxn){

	if(N<1){return;}
	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	reaction_edges_kernel <<< blks, NUM_THREADS >>> ( N, d_edges, d_rxn);
	cudaThreadSynchronize();

}

