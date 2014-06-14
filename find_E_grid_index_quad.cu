#include <cuda.h>
#include <stdio.h>
#include "datadef.h"

__global__ void find_E_grid_index_quad_kernel(unsigned N, unsigned depth, unsigned width, unsigned* active, qnode* qnode_root, float* E, unsigned * index, unsigned* done){

	int tid = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid >= N){return;}

	//remap to active
	//tid=active[tid];
	if(done[tid]){return;}

	//load data
	float this_E = E[tid];
	unsigned dex=0;

	//find initial node bin
	//unsigned 	bin  = tid*width/N;
	qnode 		node = *qnode_root;
	qnode*		next_ptr;

	//return;
	// traverse tree
	for(unsigned it=0; it<depth ;it++){
		//printf("tid=%d values=(%6.4E %6.4E %6.4E %6.4E)\n",tid,node.values[0],node.values[1],node.values[2],node.values[3]);
		if(this_E <= node.values[1]){
			next_ptr = node.leaves[0];
		}
		else if(this_E <= node.values[2]){
			next_ptr = node.leaves[1];
		}
		else if(this_E <= node.values[3]){
			next_ptr = node.leaves[2];
		}
		else{
			next_ptr = node.leaves[3];
		}
		memcpy(&node,next_ptr,sizeof(qnode));
	}

	// lowest level, copy index intead of pointer
	if(this_E <= node.values[1]){
		dex=(unsigned long)node.leaves[0];
	}
	else if(this_E <= node.values[2]){
		dex=(unsigned long)node.leaves[1];
	}
	else if(this_E <= node.values[3]){
		dex=(unsigned long)node.leaves[2];
	}
	else{
		dex=(unsigned long)node.leaves[3];
	}

	//write output index
	index[tid]=dex;

}

void find_E_grid_index_quad( unsigned NUM_THREADS, unsigned N, unsigned depth, unsigned width, unsigned* active, qnode* tree, float* E, unsigned * index, unsigned* done){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	find_E_grid_index_quad_kernel <<< blks, NUM_THREADS >>> (  N,  depth,  width, active, tree, E, index, done);
	cudaThreadSynchronize();

}

