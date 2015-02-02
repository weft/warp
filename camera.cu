#include "optix.h"
#include <optix_world.h>
#include "datadef.h"

using namespace optix;

rtBuffer<source_point,1>            positions_buffer;
rtBuffer<unsigned,1>                rxn_buffer;
rtBuffer<unsigned,1>                remap_buffer;
rtBuffer<unsigned,1>                done_buffer;
rtBuffer<unsigned,1>                cellnum_buffer;
rtBuffer<unsigned,1>                matnum_buffer;
rtDeclareVariable(rtObject,      top_object, , );
rtDeclareVariable(uint, launch_index_in, rtLaunchIndex, );
rtDeclareVariable(uint, launch_dim,   rtLaunchDim, );
rtDeclareVariable(unsigned,  outer_cell, , );
rtDeclareVariable(unsigned,  trace_type, , );
rtDeclareVariable(unsigned,  boundary_condition, , );

RT_PROGRAM void camera()
{
	//skip done particles

	//remap if 2
	unsigned launch_index;
	if(trace_type==2){
		launch_index=remap_buffer[launch_index_in];
		if(rxn_buffer[launch_index_in]>900){return;}
	}
	else{
		launch_index = launch_index_in;
	}

	// declare important stuff
	int                 cnt;
	float               epsilon=1e-5; 	
	intersection_point  payload;
	
	// init payload flags
	payload.cont=1;
	payload.buff_index=0;
	
	// init ray
	float3 ray_direction  = make_float3(positions_buffer[launch_index].xhat, positions_buffer[launch_index].yhat, positions_buffer[launch_index].zhat);
	float3 ray_origin     = make_float3(positions_buffer[launch_index].x,    positions_buffer[launch_index].y,    positions_buffer[launch_index].z);
	optix::Ray ray        = optix::make_Ray( ray_origin, ray_direction, 0, epsilon, RT_DEFAULT_MAX );

	// init hitbuff
	for(cnt=0;cnt<10;cnt++){
		payload.hitbuff[cnt].cell = -1;
		payload.hitbuff[cnt].mat  = -1;
		payload.hitbuff[cnt].fiss = -1;
	}

	// first trace to find closest hit, set bc flag
	rtTrace(top_object, ray, payload);
	//rtPrintf("did first trace, type %u\n",trace_type);
	positions_buffer[launch_index].surf_dist = payload.surf_dist; 
	if(trace_type==2){
		if(payload.hitbuff[0].cell == outer_cell){
			positions_buffer[launch_index].enforce_BC=1;
			//rtPrintf("should enforce BC\n");
		}
		else{
			positions_buffer[launch_index].enforce_BC=0;
			//rtPrintf("DO NOT enforce BC\n");
		}
	}

	payload.hitbuff[0].cell = -1;
	payload.hitbuff[0].mat  = -1;
	payload.hitbuff[0].fiss = -1;
	payload.buff_index=0;
	payload.cont=1;

	// find entering cell otherwise, trace will write, use downward z 
	ray_direction  = make_float3(0, 0, -1.0);
	ray_origin     = make_float3(positions_buffer[launch_index].x,    positions_buffer[launch_index].y,    positions_buffer[launch_index].z);
	ray = optix::make_Ray( ray_origin, ray_direction, 0, epsilon, RT_DEFAULT_MAX );
	rtTrace(top_object, ray, payload); 
	while(payload.cont){
		ray_origin = make_float3(payload.x+1.0*epsilon,payload.y+1.0*epsilon,payload.z+1.0*epsilon);
		ray = optix::make_Ray( ray_origin, ray_direction, 0, epsilon, RT_DEFAULT_MAX );
		rtTrace(top_object, ray, payload);      
	}
	if(trace_type == 2){ //write material to buffer normally, write surface distance
		matnum_buffer[launch_index] 				= payload.hitbuff[0].mat;
		cellnum_buffer[launch_index] 				= payload.hitbuff[0].cell;
	}
	else if(trace_type == 3){  //write fissile flag if fissile query
		matnum_buffer[launch_index] 				= payload.hitbuff[0].fiss;
		cellnum_buffer[launch_index] 				= payload.hitbuff[0].cell;
		rxn_buffer[launch_index_in] 				= 818;
	}

}

RT_PROGRAM void exception()
{
	const unsigned int code = rtGetExceptionCode();
	rtPrintf( "Caught exception 0x%X at launch index (%d)\n", code, launch_index_in);
	rtPrintExceptionDetails();
}
