#include "optix.h"
#include <optix_world.h>
#include "datadef.h"

using namespace optix;

rtBuffer<spatial_data,1>            positions_buffer;
rtBuffer<unsigned,1>                rxn_buffer;
rtBuffer<unsigned,1>                remap_buffer;
rtBuffer<unsigned,1>                cellnum_buffer;
rtBuffer<unsigned,1>                matnum_buffer;
rtBuffer<unsigned,1>                talnum_buffer;
rtDeclareVariable(rtObject,		top_object, 		, 				);
rtDeclareVariable(uint, 		launch_index_in, 	rtLaunchIndex, 	);
rtDeclareVariable(uint, 		launch_dim,   		rtLaunchDim, 	);
rtDeclareVariable(unsigned,  	outer_cell, 		, 				);
rtDeclareVariable(unsigned,  	trace_type, 		, 				);
rtDeclareVariable(unsigned,  	boundary_condition, , 				);

RT_PROGRAM void camera()
{
	//skip done particles
//	rtPrintf("top of camera\n");

	//remap if 2
	unsigned launch_index;
	if(trace_type==2){
		launch_index=remap_buffer[launch_index_in];

//		rtPrintf("launch index %u rxn %u\n",launch_index,rxn_buffer[launch_index_in]);

		if(rxn_buffer[launch_index_in]>=900){return;}
	}
	else{
		launch_index = launch_index_in;
	}

	// declare important stuff
	int                 sense = 0;
	float               epsilon=5.0e-6; 	
	float3 ray_direction, ray_origin;
	optix::Ray ray;
	intersection_point  payload;

	payload.launch_dex = launch_index;
	payload.cont	   = 1;
	payload.buff_index = 0;

	// null rxn, miss will set it if there is a miss
	// rxn_buffer[launch_index_in] = 0.0;

	//rtPrintf("ray %u rxn %u xyz-hat % 10.8E % 10.8E % 10.8E\n",launch_index,rxn_buffer[launch_index_in],positions_buffer[launch_index].xhat,positions_buffer[launch_index].yhat,positions_buffer[launch_index].zhat);

	// find nearest surface if and BC if type 2 
	if(trace_type==2){
		// init payload flags
		payload.sense 		= 0;
		payload.mat   		= 999999;
		payload.cell  		= 999999;
		payload.fiss  		= 0;
		payload.x 		= 0.0;
		payload.y 		= 0.0;
		payload.z 		= 0.0; 
		payload.surf_dist 	= 50000;  
		payload.norm[0]   	= 0.0; 
		payload.norm[1]   	= 0.0; 
		payload.norm[2]   	= 0.0;    
//		payload.sense     	= 0.0;   
//		if(launch_index==1){rtPrintf("outer cell: %u\n",outer_cell);}

		// init ray
		ray_direction  	= make_float3(positions_buffer[launch_index].xhat, positions_buffer[launch_index].yhat, positions_buffer[launch_index].zhat);
		ray_origin     	= make_float3(positions_buffer[launch_index].x,    positions_buffer[launch_index].y,    positions_buffer[launch_index].z);
		ray        = optix::make_Ray( ray_origin, ray_direction, 0, epsilon, RT_DEFAULT_MAX );
		// first trace to find closest hit, set norm/distance, set bc flag
		rtTrace(top_object, ray, payload);
		positions_buffer[launch_index].surf_dist = payload.surf_dist; 
		positions_buffer[launch_index].norm[0]   = payload.norm[0];
		positions_buffer[launch_index].norm[1]   = payload.norm[1];
		positions_buffer[launch_index].norm[2]   = payload.norm[2];
		// write bc flag if first hit is outer cell
//		if(payload.cell == outer_cell){
			positions_buffer[launch_index].enforce_BC = boundary_condition;
//		}
//		else{
//			positions_buffer[launch_index].enforce_BC = 0;
//		}
	}

	// re-init sense, payload, ray
	sense			= 0;
	payload.sense		= 0;
	payload.mat		= 999999;
	payload.cell		= 999999;
	payload.fiss		= 0;
	payload.x		= 0.0;
	payload.y		= 0.0;
	payload.z		= 0.0; 
	payload.surf_dist	= 50000;  
	payload.norm[0]		= 0.0; 
	payload.norm[1]		= 0.0; 
	payload.norm[2]		= 0.0;    
//	payload.sense		= 0.0;   
	
	ray_direction	= make_float3(0,0,-1);
	ray_origin	= make_float3(positions_buffer[launch_index].x,    
			positions_buffer[launch_index].y,    positions_buffer[launch_index].z);
//	ray_direction  	= make_float3(positions_buffer[launch_index].xhat, 
//			positions_buffer[launch_index].yhat, positions_buffer[launch_index].zhat);
	ray		= optix::make_Ray( ray_origin, ray_direction, 0, epsilon, RT_DEFAULT_MAX );
	const float	push_value		= 2.0;
	float dotp=0.0;

	// then find entering cell, use downward z to make problems with high x-y density faster
	rtTrace(top_object, ray, payload);
	sense = payload.sense;
	while(sense>=0){// & (outer_cell!=payload.cell)){
//	while((sense>=0) || (payload.cont)){// & (outer_cell!=payload.cell)){
//	while(payload.cont){
		dotp = 	payload.norm[0]*ray_direction.x +
			payload.norm[1]*ray_direction.y +
			payload.norm[2]*ray_direction.z;
		ray_origin = make_float3(
				payload.x+copysignf(1.0,dotp)*push_value*epsilon*payload.norm[0],
				payload.y+copysignf(1.0,dotp)*push_value*epsilon*payload.norm[1],
				payload.z+copysignf(1.0,dotp)*push_value*epsilon*payload.norm[2]);
		ray = optix::make_Ray( ray_origin, ray_direction, 0, epsilon, RT_DEFAULT_MAX );
		rtTrace(top_object, ray, payload);
		sense = sense + payload.sense;
	}

	matnum_buffer[launch_index] = payload.mat;
	cellnum_buffer[launch_index] = payload.cell;
	talnum_buffer[launch_index] = payload.tally_index;

	if(trace_type == 3)
	{
		rxn_buffer[launch_index_in] = 0;
		matnum_buffer[launch_index] = payload.fiss;
	}

	if(trace_type == 2)
	{
		for(int i = 0; i < 10; i++)
		{
			positions_buffer[launch_index].cell[i] = -1;
			positions_buffer[launch_index].dist[i] = -1.0;
			positions_buffer[launch_index].mat[i] = -1;
		}
		payload.buff_index	= 0;
		payload.cont		= 1;
		ray_direction  	= make_float3(positions_buffer[launch_index].xhat, 
				positions_buffer[launch_index].yhat, positions_buffer[launch_index].zhat);
		ray_origin     	= make_float3(positions_buffer[launch_index].x,    
				positions_buffer[launch_index].y,    positions_buffer[launch_index].z);
		ray             = optix::make_Ray( ray_origin, ray_direction, 0, epsilon, RT_DEFAULT_MAX);
		rtTrace(top_object, ray, payload);
		while(payload.cont)
		{
			dotp = 	payload.norm[0]*ray_direction.x +
				payload.norm[1]*ray_direction.y +
				payload.norm[2]*ray_direction.z;
			ray_origin = make_float3(
					payload.x+copysignf(1.0,dotp)*push_value*epsilon*payload.norm[0],
					payload.y+copysignf(1.0,dotp)*push_value*epsilon*payload.norm[1],
					payload.z+copysignf(1.0,dotp)*push_value*epsilon*payload.norm[2]);
			ray = optix::make_Ray( ray_origin, ray_direction, 0, epsilon, RT_DEFAULT_MAX );
			rtTrace(top_object, ray, payload);
		}
	}

}

RT_PROGRAM void exception()
{
	const unsigned int code = rtGetExceptionCode();
	rtPrintf( "Caught exception 0x%X at launch index (%d)\n", code, launch_index_in);
	rtPrintExceptionDetails();
}
