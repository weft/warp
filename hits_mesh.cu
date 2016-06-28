#include <optix.h>
#include <optixu/optixu_math_namespace.h>
#include "datadef.h"

using namespace optix;

rtDeclareVariable(intersection_point, payload, rtPayload, ); 
rtDeclareVariable(float, int_dist, rtIntersectionDistance, );
rtDeclareVariable(optix::Ray, ray, rtCurrentRay, );
rtDeclareVariable(unsigned,  cellnum,     attribute cell_num, );
rtDeclareVariable(int,       celltal,     attribute cell_tal, );
rtDeclareVariable(unsigned,  cellmat,     attribute cell_mat, );
rtDeclareVariable(unsigned,  cellfissile, attribute cell_fis, );
rtDeclareVariable(unsigned,  sense,       attribute cell_sense, );
rtDeclareVariable(float3,    normal,      attribute normal, );
rtDeclareVariable(unsigned,  outer_cell,  , );
rtDeclareVariable(uint, launch_index, rtLaunchIndex, );
rtBuffer<spatial_data,1> positions_buffer;

RT_PROGRAM void closest_hit()
{

	// always update current position and intersection distance, camera takes care of recording the first one
	payload.x         = int_dist * ray.direction.x + ray.origin.x;
	payload.y         = int_dist * ray.direction.y + ray.origin.y;
	payload.z         = int_dist * ray.direction.z + ray.origin.z;
	payload.surf_dist = int_dist;
	
	//write normals
	payload.norm[0] = normal.x;
	payload.norm[1] = normal.y;
	payload.norm[2] = normal.z;

	// update sense
	payload.sense = sense;
	if (sense == 0){rtPrintf("sense of closest_hit is 0!\n");}

	if(cellnum == outer_cell){ payload.cont = 0; }

	//update mat, cell, fiss
	payload.mat		= cellmat;
	payload.cell		= cellnum;
	payload.tally_index	= celltal;
	payload.fiss		= cellfissile;

	positions_buffer[payload.launch_dex].cell[payload.buff_index] = cellnum;
	positions_buffer[payload.launch_dex].dist[payload.buff_index] = int_dist;
	positions_buffer[payload.launch_dex].mat[payload.buff_index] = cellmat;

	payload.buff_index++;
}
