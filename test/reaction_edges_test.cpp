#include "gtest/gtest.h"
#include <limits.h>
#include <cuda.h>
#include <stdio.h>
#include <vector> 
#include <iostream>
#include <sstream>
#include <cmath>
#include <assert.h>
#include <time.h>
#include <string.h>
#include "device_copies.h"

void reaction_edges( unsigned ,  unsigned , unsigned* , unsigned* );

namespace {

// The fixture for testing class Foo.
class ReactionEdgesTest : public ::testing::Test {

	protected:

	// data
	unsigned 	NUM_THREADS;
	unsigned 	N;
	unsigned 	n_edges;
	unsigned* 	edges;
	unsigned* 	edges_solution;
	unsigned* 	rxn;
	unsigned* 	d_edges;
	unsigned* 	d_rxn;	

	// edges array format is:
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


	ReactionEdgesTest() {
		// set defaults
		NUM_THREADS 	= 256;
		N 				= 1000;
		n_edges 		= 11;		
	}

	virtual ~ReactionEdgesTest() {
		// free mem
		delete rxn;
		delete edges;
		delete edges_solution;
		deallocate_on_device((void*)d_edges);
		deallocate_on_device((void*)d_rxn);
	}

	void allocate(){
		rxn						= new unsigned [N];
		edges_solution 			= new unsigned [n_edges];
		edges 					= new unsigned [n_edges];
		allocate_on_device((void**)&d_edges,n_edges*sizeof(unsigned));
		allocate_on_device((void**)&d_rxn,N*sizeof(unsigned));	
	}


	void compute_solution(){
	
		// do first element
		if ( rxn[0] == 2 ){
			edges_solution[0] = 1;
		}
		else if( rxn[0] >=51 && rxn[0] <= 90 ){
			edges_solution[2] = 1;
		}
		else if( rxn[0] == 91 ){
			edges_solution[4] = 1;
		} 
		else if( rxn[0] == 800 ){
			edges_solution[6] = 1;
		}
		else if( rxn[0] >= 811 && rxn[0] <=845 ){
			edges_solution[8] = 1;
		}
		else if(rxn[0] >= 900 ){
			edges_solution[10] = 1;   // [10] = 1 is basically a done flag
		}

		// scan array, mark edges
		for(unsigned i=0;i<N-1;i++){

			// 2-edge
			if( rxn[i]==2 ){
				if( rxn[i+1] >=51 && rxn[i+1] <= 90){
					edges_solution[1] = i+1;
					edges_solution[2] = i+2;
				}
				else if( rxn[i+1] == 91 ){
					edges_solution[1] = i+1;
					edges_solution[4] = i+2;
				}
				else if( rxn[i+1] == 800 ){
					edges_solution[1] = i+1;
					edges_solution[6] = i+2;
				}
				else if( rxn[i+1] >=811 && rxn[i+1] <= 845 ){
					edges_solution[1] = i+1;
					edges_solution[8] = i+2;
				}
				else if( rxn[i+1] >= 900 ){
					edges_solution[1]  = i+1;
					edges_solution[10] = i+2;
				}
			}

			// 51-90 edge
			if( rxn[i] >=51 && rxn[i] <= 90 ){
				if( rxn[i+1] == 91 ){
					edges_solution[3] = i+1;
					edges_solution[4] = i+2;
				}
				else if( rxn[i+1] == 800 ){
					edges_solution[3] = i+1;
					edges_solution[6] = i+2;
				}
				else if( rxn[i+1] >=811 && rxn[i+1] <= 845 ){
					edges_solution[3] = i+1;
					edges_solution[8] = i+2;
				}
				else if( rxn[i+1] >= 900 ){
					edges_solution[3]  = i+1;
					edges_solution[10] = i+2;
				}
			}

			// 91 edge
			if( rxn[i] == 91 ){
				if( rxn[i+1] == 800 ){
					edges_solution[5] = i+1;
					edges_solution[6] = i+2;
				}
				else if( rxn[i+1] >=811 && rxn[i+1] <= 845 ){
					edges_solution[5] = i+1;
					edges_solution[8] = i+2;
				}
				else if( rxn[i+1] >= 900 ){
					edges_solution[5]  = i+1;
					edges_solution[10] = i+2;
				}
			}

			// 800 edge
			if( rxn[i] == 800 ){
				if( rxn[i+1] >=811 && rxn[i+1] <= 845 ){
					edges_solution[7] = i+1;
					edges_solution[8] = i+2;
				}
				else if( rxn[i+1] >= 900 ){
					edges_solution[7]  = i+1;
					edges_solution[10] = i+2;
				}
			}

			// 811-845 edge
			if( rxn[i] >=811 && rxn[i] <= 845 ){
				if( rxn[i+1] >= 900 ){
					edges_solution[9]  = i+1;
					edges_solution[10] = i+2;
				}
			}
		}

		// do last element
		if ( rxn[N-1] == 2 ){
			edges_solution[1] = N;
		}
		else if( rxn[N-1] >=51 && rxn[N-1] <= 90 ){
			edges_solution[3] = N;
		}
		else if( rxn[N-1] == 91 ){
			edges_solution[5] = N;
		} 
		else if( rxn[N-1] == 800 ){
			edges_solution[7] = N;
		}
		else if( rxn[N-1] >= 811 && rxn[N-1] <=845 ){
			edges_solution[9] = N;
		}
		else if(rxn[N-1] >= 900 ){
			//edges_solution[10] = 1;   // [10] = 1 is basically a done flag
		}
	}

	void run_test(){

		// init to zeros
		for(unsigned i=0;i<11;i++){
			edges[i] = 0;
			edges_solution[i] = 0 ;
		}

		// copy set data
		copy_to_device( d_edges , edges , n_edges*sizeof(unsigned));
		copy_to_device( d_rxn   , rxn   ,       N*sizeof(unsigned));

		// compute solution serially
		compute_solution();

		// compute solution on device
		reaction_edges(  NUM_THREADS,   N,  d_edges,  d_rxn);
		
		// copy device solution back
		copy_from_device( edges , d_edges , n_edges*sizeof(unsigned));

		// print stuff just for debugging
		//printf("rxn\n");
		//for(unsigned i=0; i<N ;i++){ printf("%4u\n",rxn[i]);}
		//printf("device, solution\n");
		//for(unsigned i=0; i<n_edges ;i++){ printf("%6u %6u\n",edges[i],edges_solution[i]);}

		// check
		EXPECT_EQ(  edges_solution[0]  ,  edges[0]   );
		EXPECT_EQ(  edges_solution[1]  ,  edges[1]   );
		EXPECT_EQ(  edges_solution[2]  ,  edges[2]   );
		EXPECT_EQ(  edges_solution[3]  ,  edges[3]   );
		EXPECT_EQ(  edges_solution[4]  ,  edges[4]   );
		EXPECT_EQ(  edges_solution[5]  ,  edges[5]   );
		EXPECT_EQ(  edges_solution[6]  ,  edges[6]   );
		EXPECT_EQ(  edges_solution[7]  ,  edges[7]   );
		EXPECT_EQ(  edges_solution[8]  ,  edges[8]   );
		EXPECT_EQ(  edges_solution[9]  ,  edges[9]   );
		EXPECT_EQ(  edges_solution[10] ,  edges[10]  );

	}

};  // class

// TESTS GO HERE

//
//    TEST WITH ALL LENGTHS > 2
//

TEST_F(ReactionEdgesTest, LongAll){
	// this test does XXX
	
	// set size, allocate
	N = 1000;
	unsigned c;
	allocate();

	// set rxn for this case
	for(unsigned i=0; i<100; i++){
		rxn[i] = 2;
	}
	for(unsigned i=100; i<400; i++){
		c = (i-100)/7;
		if( c > 39 ){c=39;}
		rxn[i] = 51+c;
	}
	for(unsigned i=400; i<500; i++){
		rxn[i] = 91;
	}
	for(unsigned i=500; i<600; i++){
		rxn[i] = 800;
	}
	for(unsigned i=600; i<800; i++){
		c = (i-600)/5;
		if( c > 34 ){c=34;}
		rxn[i] = 811+c;
	}
	for(unsigned i=800; i<N; i++){
		rxn[i] = 900+(i-800);
	}

	// run the tests for this reaction vector
	run_test();
}   
TEST_F(ReactionEdgesTest, LongNoEscatter){
	// this test does XXX
	
	// set size, allocate
	N = 1000;
	unsigned c;
	allocate();

	// set rxn for this case
	for(unsigned i=0; i<400; i++){
		c = (i-0)/8;
		if( c > 39 ){c=39;}
		rxn[i] = 51+c;
	}
	for(unsigned i=400; i<500; i++){
		rxn[i] = 91;
	}
	for(unsigned i=500; i<600; i++){
		rxn[i] = 800;
	}
	for(unsigned i=600; i<800; i++){
		c = (i-600)/5;
		if( c > 34 ){c=34;}
		rxn[i] = 811+c;
	}
	for(unsigned i=800; i<N; i++){
		rxn[i] = 900+(i-800);
	}

	// run the tests for this reaction vector
	run_test();
}   
TEST_F(ReactionEdgesTest, LongNoIscatter){
	// this test does XXX
	
	// set size, allocate
	N = 1000;
	unsigned c;
	allocate();

	// set rxn for this case
	for(unsigned i=0; i<400; i++){
		rxn[i] = 2;
	}
	for(unsigned i=400; i<500; i++){
		rxn[i] = 91;
	}
	for(unsigned i=500; i<600; i++){
		rxn[i] = 800;
	}
	for(unsigned i=600; i<800; i++){
		c = (i-600)/5;
		if( c > 34 ){c=34;}
		rxn[i] = 811+c;
	}
	for(unsigned i=800; i<N; i++){
		rxn[i] = 900+(i-800);
	}

	// run the tests for this reaction vector
	run_test();
} 
TEST_F(ReactionEdgesTest, LongNoCscatter){
	// this test does XXX
	
	// set size, allocate
	N = 1000;
	unsigned c;
	allocate();

	// set rxn for this case
	for(unsigned i=0; i<100; i++){
		rxn[i] = 2;
	}
	for(unsigned i=100; i<400; i++){
		c = (i-100)/7;
		if( c > 39 ){c=39;}
		rxn[i] = 51+c;
	}
	for(unsigned i=400; i<600; i++){
		rxn[i] = 800;
	}
	for(unsigned i=600; i<800; i++){
		c = (i-600)/5;
		if( c > 34 ){c=34;}
		rxn[i] = 811+c;
	}
	for(unsigned i=800; i<N; i++){
		rxn[i] = 900+(i-800);
	}

	// run the tests for this reaction vector
	run_test();
}  
TEST_F(ReactionEdgesTest, LongNoResample){
	// this test does XXX
	
	// set size, allocate
	N = 1000;
	unsigned c;
	allocate();

	// set rxn for this case
	for(unsigned i=0; i<100; i++){
		rxn[i] = 2;
	}
	for(unsigned i=100; i<400; i++){
		c = (i-100)/7;
		if( c > 39 ){c=39;}
		rxn[i] = 51+c;
	}
	for(unsigned i=400; i<600; i++){
		rxn[i] = 91;
	}
	for(unsigned i=600; i<800; i++){
		c = (i-600)/5;
		if( c > 34 ){c=34;}
		rxn[i] = 811+c;
	}
	for(unsigned i=800; i<N; i++){
		rxn[i] = 900+(i-800);
	}

	// run the tests for this reaction vector
	run_test();
}      
TEST_F(ReactionEdgesTest, LongNoFission){
	// this test does XXX
	
	// set size, allocate
	N = 1000;
	unsigned c;
	allocate();

	// set rxn for this case
	for(unsigned i=0; i<100; i++){
		rxn[i] = 2;
	}
	for(unsigned i=100; i<400; i++){
		c = (i-100)/7;
		if( c > 39 ){c=39;}
		rxn[i] = 51+c;
	}
	for(unsigned i=400; i<500; i++){
		rxn[i] = 91;
	}
	for(unsigned i=500; i<800; i++){
		rxn[i] = 800;
	}
	for(unsigned i=800; i<N; i++){
		rxn[i] = 900+(i-800);
	}

	// run the tests for this reaction vector
	run_test();
}   
TEST_F(ReactionEdgesTest, LongNoDone){
	// this test does XXX
	
	// set size, allocate
	N = 1000;
	unsigned c;
	allocate();

	// set rxn for this case
	for(unsigned i=0; i<100; i++){
		rxn[i] = 2;
	}
	for(unsigned i=100; i<400; i++){
		c = (i-100)/7;
		if( c > 39 ){c=39;}
		rxn[i] = 51+c;
	}
	for(unsigned i=400; i<500; i++){
		rxn[i] = 91;
	}
	for(unsigned i=500; i<600; i++){
		rxn[i] = 800;
	}
	for(unsigned i=600; i<N; i++){
		c = (i-600)/5;
		if( c > 34 ){c=34;}
		rxn[i] = 811+c;
	}

	// run the tests for this reaction vector
	run_test();
}  

TEST_F(ReactionEdgesTest, LongOnlyEscatter){
	// this test does XXX
	
	// set size, allocate
	N = 1000;
	allocate();

	// set rxn for this case
	for(unsigned i=0; i<N; i++){
		rxn[i] = 2;
	}

	// run the tests for this reaction vector
	run_test();
}    
TEST_F(ReactionEdgesTest, LongOnlyIscatter){
	// this test does XXX
	
	// set size, allocate
	N = 1000;
	unsigned c;
	allocate();

	// set rxn for this case
	for(unsigned i=0; i<N; i++){
		c = i/25;
		if( c > 39 ){c=39;}
		rxn[i] = 51+c;
	}
	
	// run the tests for this reaction vector
	run_test();
}   
TEST_F(ReactionEdgesTest, LongOnlyCscatter){
	// this test does XXX
	
	// set size, allocate
	N = 1000;
	allocate();

	// set rxn for this case
	for(unsigned i=0; i<N; i++){
		rxn[i] = 91;
	}
	
	// run the tests for this reaction vector
	run_test();
}  
TEST_F(ReactionEdgesTest, LongOnlyResample){
	// this test does XXX
	
	// set size, allocate
	N = 1000;
	allocate();

	// set rxn for this case
	for(unsigned i=0; i<N; i++){
		rxn[i] = 800;
	}
	
	// run the tests for this reaction vector
	run_test();
}    
TEST_F(ReactionEdgesTest, LongOnlyFission){
	// this test does XXX
	
	// set size, allocate
	N = 1000;
	unsigned c;
	allocate();

	// set rxn for this case
	for(unsigned i=0; i<N; i++){
		c = i/29;
		if( c > 34 ){c=34;}
		rxn[i] = 811+c;
	}
	
	// run the tests for this reaction vector
	run_test();
}   
TEST_F(ReactionEdgesTest, LongOnlyDone){
	// this test does XXX
	
	// set size, allocate
	N = 1000;
	unsigned c;
	allocate();

	// set rxn for this case
	for(unsigned i=0; i<N; i++){
		c = i/5;
		rxn[i] = 900+c;
	}
	
	// run the tests for this reaction vector
	run_test();
}   
//
//    Tests where the block len == 2
//

TEST_F(ReactionEdgesTest, TwoAll){
	// this test does XXX
	
	// set size, allocate
	N = 12;
	allocate();

	// set rxn for this case
	rxn[0] = 2;
	rxn[1] = 2;
	rxn[2] = 51;
	rxn[3] = 58;
	rxn[4] = 91;
	rxn[5] = 91;
	rxn[6] = 800;
	rxn[7] = 800;
	rxn[8] = 811;
	rxn[9] = 816;
	rxn[10] = 900;
	rxn[11] = 1102;

	// run the tests for this reaction vector
	run_test();
}   
TEST_F(ReactionEdgesTest, TwoNoEscatter){
	// this test does XXX
	
	// set size, allocate
	N = 10;
	allocate();

	// set rxn for this case
	rxn[0] = 51;
	rxn[1] = 58;
	rxn[2] = 91;
	rxn[3] = 91;
	rxn[4] = 800;
	rxn[5] = 800;
	rxn[6] = 811;
	rxn[7] = 816;
	rxn[8] = 900;
	rxn[9] = 1102;
	
	// run the tests for this reaction vector
	run_test();
}  
TEST_F(ReactionEdgesTest, TwoNoIscatter){
	// this test does XXX
	
	// set size, allocate
	N = 10;
	allocate();

	// set rxn for this case
	rxn[0] = 2;
	rxn[1] = 2;
	rxn[2] = 91;
	rxn[3] = 91;
	rxn[4] = 800;
	rxn[5] = 800;
	rxn[6] = 811;
	rxn[7] = 816;
	rxn[8] = 900;
	rxn[9] = 1102;
	
	// run the tests for this reaction vector
	run_test();
} 
TEST_F(ReactionEdgesTest, TwoNoCscatter){
	// this test does XXX
	
	// set size, allocate
	N = 10;
	allocate();

	// set rxn for this case
	rxn[0] = 2;
	rxn[1] = 2;
	rxn[2] = 51;
	rxn[3] = 90;
	rxn[4] = 800;
	rxn[5] = 800;
	rxn[6] = 811;
	rxn[7] = 816;
	rxn[8] = 900;
	rxn[9] = 1102;
	
	// run the tests for this reaction vector
	run_test();
}  
TEST_F(ReactionEdgesTest, TwoNoResample){
	// this test does XXX
	
	// set size, allocate
	N = 10;
	allocate();

	// set rxn for this case
	rxn[0] = 2;
	rxn[1] = 2;
	rxn[2] = 51;
	rxn[3] = 90;
	rxn[4] = 91;
	rxn[5] = 91;
	rxn[6] = 811;
	rxn[7] = 816;
	rxn[8] = 900;
	rxn[9] = 1102;
	
	// run the tests for this reaction vector
	run_test();
}       
TEST_F(ReactionEdgesTest, TwoNoFission){
	// this test does XXX
	
	// set size, allocate
	N = 10;
	allocate();

	// set rxn for this case
	rxn[0] = 2;
	rxn[1] = 2;
	rxn[2] = 51;
	rxn[3] = 90;
	rxn[4] = 91;
	rxn[5] = 91;
	rxn[6] = 800;
	rxn[7] = 800;
	rxn[8] = 900;
	rxn[9] = 1102;
	
	// run the tests for this reaction vector
	run_test();
}   
TEST_F(ReactionEdgesTest, TwoNoDone){
	// this test does XXX
	
	// set size, allocate
	N = 10;
	allocate();

	// set rxn for this case
	rxn[0] = 2;
	rxn[1] = 2;
	rxn[2] = 51;
	rxn[3] = 90;
	rxn[4] = 91;
	rxn[5] = 91;
	rxn[6] = 800;
	rxn[7] = 800;
	rxn[8] = 811;
	rxn[9] = 840;
	
	// run the tests for this reaction vector
	run_test();
} 

//
//    Tests where the block len == 1
//

TEST_F(ReactionEdgesTest, OneAll){
	// this test does XXX
	
	// set size, allocate
	N = 6;
	allocate();

	// set rxn for this case
	rxn[0] = 2;
	rxn[1] = 51;
	rxn[2] = 91;
	rxn[3] = 800;
	rxn[4] = 811;
	rxn[5] = 1102;

	// run the tests for this reaction vector
	run_test();
}   
TEST_F(ReactionEdgesTest, OneNoEscatter){
	// this test does XXX
	
	// set size, allocate
	N = 5;
	allocate();

	// set rxn for this case
	rxn[0] = 51;
	rxn[1] = 91;
	rxn[2] = 800;
	rxn[3] = 811;
	rxn[4] = 1102;
	
	// run the tests for this reaction vector
	run_test();
}  
TEST_F(ReactionEdgesTest, OneNoIscatter){
	// this test does XXX
	
	// set size, allocate
	N = 5;
	allocate();

	// set rxn for this case
	rxn[0] = 2;
	rxn[1] = 91;
	rxn[2] = 800;
	rxn[3] = 811;
	rxn[4] = 1102;
	
	// run the tests for this reaction vector
	run_test();
} 
TEST_F(ReactionEdgesTest, OneNoCscatter){
	// this test does XXX
	
	// set size, allocate
	N = 5;
	allocate();

	// set rxn for this case
	rxn[0] = 2;
	rxn[1] = 51;
	rxn[2] = 800;
	rxn[3] = 811;
	rxn[4] = 1102;
	
	// run the tests for this reaction vector
	run_test();
}  
TEST_F(ReactionEdgesTest, OneNoResample){
	// this test does XXX
	
	// set size, allocate
	N = 5;
	allocate();

	// set rxn for this case
	rxn[0] = 2;
	rxn[1] = 51;
	rxn[2] = 91;
	rxn[3] = 816;
	rxn[4] = 1102;
	
	// run the tests for this reaction vector
	run_test();
}       
TEST_F(ReactionEdgesTest, OneNoFission){
	// this test does XXX
	
	// set size, allocate
	N = 5;
	allocate();

	// set rxn for this case
	rxn[0] = 2;
	rxn[1] = 51;
	rxn[2] = 91;
	rxn[3] = 800;
	rxn[4] = 1102;
	
	// run the tests for this reaction vector
	run_test();
}   
TEST_F(ReactionEdgesTest, OneNoDone){
	// this test does XXX
	
	// set size, allocate
	N = 5;
	allocate();

	// set rxn for this case
	rxn[0] = 2;
	rxn[1] = 51;
	rxn[2] = 91;
	rxn[3] = 800;
	rxn[4] = 816;
	
	// run the tests for this reaction vector
	run_test();
} 

//
//   short tests
//

TEST_F(ReactionEdgesTest, ShortEscatter){
	// this test does XXX
	
	// set size, allocate
	N = 1;
	allocate();

	// set rxn for this case
	rxn[0] = 2;
	
	// run the tests for this reaction vector
	run_test();
} 
TEST_F(ReactionEdgesTest, ShortIscatter){
	// this test does XXX
	
	// set size, allocate
	N = 1;
	allocate();

	// set rxn for this case
	rxn[0] = 51;
	
	// run the tests for this reaction vector
	run_test();
} 
TEST_F(ReactionEdgesTest, ShortCscatter){
	// this test does XXX
	
	// set size, allocate
	N = 1;
	allocate();

	// set rxn for this case
	rxn[0] = 91;
	
	// run the tests for this reaction vector
	run_test();
} 
TEST_F(ReactionEdgesTest, ShortResample){
	// this test does XXX
	
	// set size, allocate
	N = 1;
	allocate();

	// set rxn for this case
	rxn[0] = 800;
	
	// run the tests for this reaction vector
	run_test();
} 
TEST_F(ReactionEdgesTest, ShortFission){
	// this test does XXX
	
	// set size, allocate
	N = 1;
	allocate();

	// set rxn for this case
	rxn[0] = 816;
	
	// run the tests for this reaction vector
	run_test();
} 
TEST_F(ReactionEdgesTest, ShortDone){
	// this test does XXX
	
	// set size, allocate
	N = 1;
	allocate();

	// set rxn for this case
	rxn[0] = 919;
	
	// run the tests for this reaction vector
	run_test();
} 

//
//   Pathological tests?
//   Things like list being unsorted should not happen.  Should be tested seperately from edge detection.
//

}  // namespace

 int main(int argc, char **argv) {

  	::testing::InitGoogleTest(&argc, argv);
  	return RUN_ALL_TESTS();

}