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
rtDeclareVariable(unsigned,  sense      , attribute cell_sense, );
rtDeclareVariable(float3,    normal,      attribute normal,   );

RT_PROGRAM void intersect(int object_dex)
{
	float t1, t2, sdisc, sgn;
	float3 this_norm1, this_norm2, int1, int2;

	float3 loc  = make_float3(dims[object_dex].loc[0],dims[object_dex].loc[1],dims[object_dex].loc[2]);
	float3 xformed_origin = ray.origin - loc;

	float radius    = dims[object_dex].max[0];
	float zmin 		= dims[object_dex].min[2]-loc.z;
	float zmax 		= dims[object_dex].max[2]-loc.z;

	float a =        ( ray.direction.x  * ray.direction.x  ) + ( ray.direction.y  * ray.direction.y  );
	float b = 2.0 * (( ray.direction.x  * xformed_origin.x ) + ( ray.direction.y  * xformed_origin.y ));
	float c =        ( xformed_origin.x * xformed_origin.x ) + ( xformed_origin.y * xformed_origin.y ) - (radius * radius);
	float disc = (b*b)-(4.0*a*c);

	bool report = false;
	bool check_second = true;
	sgn = 0.0;

	if (disc > 0.0f){  // the line intersects the circle

		report = true;

		sdisc = sqrt(disc);
		t1 = (-b-sdisc)/(2.0*a);
		t2 = (-b+sdisc)/(2.0*a);

		// intersection points
		int1 = ray.direction * t1 + xformed_origin;
		int2 = ray.direction * t2 + xformed_origin;

		// normals point to radius
		this_norm1 = make_float3(  int1.x/radius ,  int1.y/radius , 0 );
		this_norm2 = make_float3(  int2.x/radius ,  int2.y/radius , 0 );

		//miss in corners
		if( ((int1.z > zmax) & (int2.z > zmax)) | ((int1.z < zmin) & (int2.z < zmin)) ){  
			report=false;
		}
		else{   // t1 always smaller 

			if (int1.z > zmax ){  //  top intersection z1
				t1 = (zmax - xformed_origin.z) / ray.direction.z;
				this_norm1 = make_float3(0, 0, 1);
			}
			else if(int1.z < zmin ) { // bottom intersection z1
				t1 = (zmin - xformed_origin.z) / ray.direction.z;
				this_norm1 = make_float3(0, 0, -1);
			}

			if (int2.z > zmax){  //  top intersection z2
				t2 = (zmax - xformed_origin.z) / ray.direction.z;
				this_norm2 = make_float3(0, 0, 1);
			}
			else if(int2.z < zmin) { // bottom intersection z2
				t2 = (zmin - xformed_origin.z) / ray.direction.z;
				this_norm2 = make_float3(0, 0, -1);
			}
		}
	    
	}
	else if( (xformed_origin.x*xformed_origin.x+xformed_origin.y*xformed_origin.y)<(radius*radius) ) {  // exactly perpendicular

		report = true;

		t1 = (zmax - xformed_origin.z) / ray.direction.z;
		t2 = (zmin - xformed_origin.z) / ray.direction.z;

		// sense
		if (t1*t2 < 0.0 ){ // neg means inside
			sgn = -1.0;
		}
		else{
			sgn =  1.0;
		}
		
		// intersection points
		int1 = ray.direction * t1 + xformed_origin;
		int2 = ray.direction * t2 + xformed_origin;	

		// bottom always -1, top always +1, sense used to flip if point is inside
		this_norm1 = make_float3(0,0,1);
		this_norm2 = make_float3(0,0,-1);
	}

	// sense
	if (t1*t2 < 0.0 ){ // neg means inside
		sgn = -1.0;
	}
	else{
		sgn =  1.0;
	}

	// report
	if (report){
		if(t1>0){
			if (rtPotentialIntersection(t1)) {
				cellnum     = dims[object_dex].cellnum;
				cellmat     = dims[object_dex].matnum;
				cellfissile = dims[object_dex].is_fissile;
				normal 		= sgn*this_norm1;
				sense       = int(sgn);
				if(rtReportIntersection(0)){
					check_second=false;
				}
			}
		}
		if(check_second & t2>0){
			if (rtPotentialIntersection(t2)) {
				cellnum     = dims[object_dex].cellnum;
				cellmat     = dims[object_dex].matnum;
				cellfissile = dims[object_dex].is_fissile;
				normal 		= sgn*this_norm2;
				sense       = int(sgn);
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
