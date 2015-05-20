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
<<<<<<< HEAD
	int                 sense = 0;
=======
	int                 sense;
>>>>>>> 3118af2c2a4c1900ecdbcb4a58f9cbf9a12d48e7
	float               epsilon=5.0e-4; 	
	intersection_point  payload;
	
	// init payload flags
	payload.sense = 0;
	payload.cell  = 999999;
	payload.mat   = 999999;
	payload.cell  = 999999;
	payload.fiss  = 0;
	
	// init ray
	float3 ray_direction  = make_float3(positions_buffer[launch_index].xhat, positions_buffer[launch_index].yhat, positions_buffer[launch_index].zhat);
	float3 ray_origin     = make_float3(positions_buffer[launch_index].x,    positions_buffer[launch_index].y,    positions_buffer[launch_index].z);
	optix::Ray ray        = optix::make_Ray( ray_origin, ray_direction, 0, epsilon, RT_DEFAULT_MAX );

	// first trace to find closest hit, set norm/distance, set bc flag
	rtTrace(top_object, ray, payload);
	sense = payload.sense;
	if(trace_type==2){
		positions_buffer[launch_index].surf_dist = payload.surf_dist; 
		positions_buffer[launch_index].norm[0]   = payload.norm[0];
		positions_buffer[launch_index].norm[1]   = payload.norm[1];
		positions_buffer[launch_index].norm[2]   = payload.norm[2];
		// write bc flag if first hit is outer cell
		if(payload.cell == outer_cell){
			positions_buffer[launch_index].enforce_BC = boundary_condition;
		}
		else{
			positions_buffer[launch_index].enforce_BC = 0;
		}
	}

	// find entering cell otherwise, trace will write, use downward z 
	ray_direction  = make_float3(0.0,0.0,-1.0);
	while(sense>=0 && payload.cell != outer_cell){
		ray_origin = make_float3(payload.x,payload.y,payload.z);
		ray = optix::make_Ray( ray_origin, ray_direction, 0, epsilon, RT_DEFAULT_MAX );
		rtTrace(top_object, ray, payload);
		sense += payload.sense;
	}

	// write cell/material numbers to buffer
	if(trace_type == 2){ //write material to buffer normally, write surface distance
		matnum_buffer[launch_index] 				= payload.mat;
		cellnum_buffer[launch_index] 				= payload.cell;
	}
	else if(trace_type == 3){  //write fissile flag if fissile query
		matnum_buffer[launch_index] 				= payload.fiss;
		cellnum_buffer[launch_index] 				= payload.cell;
		rxn_buffer[launch_index_in] 				= 818;
	}

}

RT_PROGRAM void exception()
{
	const unsigned int code = rtGetExceptionCode();
	rtPrintf( "Caught exception 0x%X at launch index (%d)\n", code, launch_index_in);
	rtPrintExceptionDetails();
}
