#include <optix.h>
#include <optixu/optixu_math_namespace.h>
#include "datadef.h"

using namespace optix;

rtDeclareVariable(intersection_point, payload, rtPayload, ); 
rtDeclareVariable(float, int_dist, rtIntersectionDistance, );
rtDeclareVariable(optix::Ray, ray, rtCurrentRay, );
rtDeclareVariable(unsigned,  cellnum,     attribute cell_num, );
rtDeclareVariable(unsigned,  cellmat,     attribute cell_mat, );
rtDeclareVariable(unsigned,  cellfissile, attribute cell_fis, );
rtDeclareVariable(unsigned,  sense,       attribute cell_sense, );
rtDeclareVariable(float3,    normal,      attribute normal, );

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

	//update mat, cell, fiss
	payload.mat  = cellmat;
	payload.cell = cellnum;
	payload.fiss = cellfissile;


}
