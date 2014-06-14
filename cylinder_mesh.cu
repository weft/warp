#include <optix.h>
#include <optixu/optixu_math_namespace.h>
#include <optixu/optixu_matrix_namespace.h>
#include <optixu/optixu_aabb_namespace.h>
#include "datadef.h"

using namespace optix;

rtBuffer<geom_data,1>               dims;
rtDeclareVariable(optix::Ray, ray, rtCurrentRay, );
rtDeclareVariable(unsigned,  cellnum,     attribute cell_num, );
rtDeclareVariable(unsigned,  cellmat,     attribute cell_mat, );
rtDeclareVariable(unsigned,  cellfissile, attribute cell_fis, );

RT_PROGRAM void intersect(int object_dex)
{
	float t1, t2, z1, z2, sdisc;

	float3 loc  = make_float3(dims[object_dex].loc[0],dims[object_dex].loc[1],dims[object_dex].loc[2]);
	float3 xformed_origin = ray.origin - loc;

	float radius    = dims[object_dex].max[0];
	float zmin 		= dims[object_dex].min[2]-loc.z;
	float zmax 		= dims[object_dex].max[2]-loc.z;

	float a =        ( ray.direction.x  * ray.direction.x  ) + ( ray.direction.y  * ray.direction.y  );
    float b = 2.0 * (( ray.direction.x  * xformed_origin.x ) + ( ray.direction.y  * xformed_origin.y ));
    float c =        ( xformed_origin.x * xformed_origin.x ) + ( xformed_origin.y * xformed_origin.y ) - (radius * radius);
    float disc = (b*b)-(4*a*c);

	bool report = false;
	bool check_second = true;

	if (disc > 0.0f){  // the line intersects the circle

		report = true;

		sdisc = sqrt(disc);
        t1 = (-b-sdisc)/(2*a);
        t2 = (-b+sdisc)/(2*a);

		z1 = ray.direction.z * t1 + xformed_origin.z;
		z2 = ray.direction.z * t2 + xformed_origin.z;

		//rtPrintf("zmin %6.4E zmax %6.4E t1 %6.4E t2 %6.4E z1 %6.4E z2 %6.4E report %u\n",zmin,zmax,t1,t2,z1,z2,report);

		if( ((z1 > zmax) & (z2 > zmax)) | ((z1 < zmin) & (z2 < zmin)) ){  //miss in corners
			report=false;
		}
		else{ 

			if (z1 > zmax ){  //  top intersection z1
				t1 = (zmax - xformed_origin.z) / ray.direction.z;
			}
			else if(z1 < zmin ) { // bottom intersection z1
				t1 = (zmin - xformed_origin.z) / ray.direction.z;
			}

			if (z2 > zmax){  //  top intersection z2
				t2 = (zmax - xformed_origin.z) / ray.direction.z;
			}
			else if(z2 < zmin) { // bottom intersection z2
				t2 = (zmin - xformed_origin.z) / ray.direction.z;
			}
		}
	    
	}
	else if( (xformed_origin.x*xformed_origin.x+xformed_origin.y*xformed_origin.y)<(radius*radius) ) {  // exactly perpendicular

		report = true;

		t1 = fminf((zmax - xformed_origin.z) / ray.direction.z , (zmin - xformed_origin.z) / ray.direction.z);
		t2 = fmaxf((zmax - xformed_origin.z) / ray.direction.z , (zmin - xformed_origin.z) / ray.direction.z);	

	}

	//rtPrintf("zmin %6.4E zmax %6.4E t1 %6.4E t2 %6.4E z1 %6.4E z2 %6.4E report %u\n",zmin,zmax,t1,t2,z1,z2,report);

	if (report){
		if(t1>0){
			if (rtPotentialIntersection(t1) ) {
				cellnum     = dims[object_dex].cellnum;
				cellmat     = dims[object_dex].matnum;
				cellfissile = dims[object_dex].is_fissile;
				if(rtReportIntersection(0)){
					check_second=false;
				}
			}
		}
		if(check_second & t2>0){
			if (rtPotentialIntersection(t2) ) {
				cellnum     = dims[object_dex].cellnum;
				cellmat     = dims[object_dex].matnum;
				cellfissile = dims[object_dex].is_fissile;
				rtReportIntersection(0);
			}
		}
	}

}

RT_PROGRAM void bounds (int object_dex, float result[6])
{
	float3 mins = make_float3(-dims[object_dex].max[0],-dims[object_dex].max[0], dims[object_dex].min[2]);  //set all to the radius
  	float3 maxs = make_float3( dims[object_dex].max[0], dims[object_dex].max[0], dims[object_dex].max[2]);
  	float3 loc  = make_float3( dims[object_dex].loc[0], dims[object_dex].loc[1], dims[object_dex].loc[2]);
	
  	optix::Aabb* aabb = (optix::Aabb*)result;
  	aabb->set(mins+loc, maxs+loc);
}
