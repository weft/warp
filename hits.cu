#include <optix.h>
#include <optixu/optixu_math_namespace.h>
#include "datadef.h"

using namespace optix;

rtDeclareVariable(intersection_point, payload, rtPayload, ); 
rtDeclareVariable(float, int_dist, rtIntersectionDistance, );
rtDeclareVariable(optix::Ray, ray, rtCurrentRay, );
rtDeclareVariable(unsigned,  cellnum, , );
rtDeclareVariable(unsigned,  cellmat, , );
rtDeclareVariable(unsigned,  cellfissile, , );
rtDeclareVariable(unsigned,  outer_cell, , );

RT_PROGRAM void closest_hit()
{
	unsigned j = 0;
	unsigned k = 0;
	unsigned notfound = 1;

	// write this cell's info into buffer element
	hit_buffer this_buff;
	this_buff.cell = cellnum;
	this_buff.mat  = cellmat;
	this_buff.fiss = cellfissile;

	// stop ray iterations if outer cell is hit
	if(cellnum==outer_cell){
		payload.cont=0;
	}

	// make end element
	hit_buffer end_buff;
	end_buff.cell = -1; 
	end_buff.mat  = -1; 
	end_buff.fiss = -1; 

	//rtPrintf("cellnum,matnum,isfiss %d %d %d\n",this_buff.cell,this_buff.mat,this_buff.fiss);

	// always update current position and intersection distance, camera takes care of recording the first one
	payload.x=int_dist*ray.direction.x+ray.origin.x;
	payload.y=int_dist*ray.direction.y+ray.origin.y;
	payload.z=int_dist*ray.direction.z+ray.origin.z;
	payload.surf_dist = int_dist;

	// scan for this cell
	for(j=0;j<payload.buff_index;j++){

		if(payload.hitbuff[j].cell == cellnum){   //found cell, remove this value and rpeserve order by shifting all down by one and appending -1's

			for(k=j ; k<9 ; k++){
				payload.hitbuff[k] 	= payload.hitbuff[k+1];
			}

			payload.hitbuff[9] = end_buff;
			payload.buff_index--;  //decrement buff_index
			notfound = 0;

			break;

		}
	
	}

	if(notfound){ 

		//check if overrrun
		if(payload.buff_index>9){

			rtPrintf("HIT BUFFER OVERRUN, NEED LARGE BUFFER.\n");
			rtThrow(RT_EXCEPTION_USER + 0);

		}
		//append at buff_index
		else{

			payload.hitbuff[payload.buff_index] = this_buff;
			payload.buff_index++;  //increment index

		}

	}


}
