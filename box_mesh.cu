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
rtDeclareVariable(float3, normal, attribute normal, );


static __device__ float3 boxnormal(float t, float3 t0, float3 t1)
{
  float3 neg = make_float3(t==t0.x?1:0, t==t0.y?1:0, t==t0.z?1:0);
  float3 pos = make_float3(t==t1.x?1:0, t==t1.y?1:0, t==t1.z?1:0);
  return pos-neg;
}

RT_PROGRAM void intersect(int object_dex)
{
  float3 mins = make_float3(dims[object_dex].min[0],dims[object_dex].min[1],dims[object_dex].min[2]);
  float3 maxs = make_float3(dims[object_dex].max[0],dims[object_dex].max[1],dims[object_dex].max[2]);
  float3 loc  = make_float3(dims[object_dex].loc[0],dims[object_dex].loc[1],dims[object_dex].loc[2]);
  float3 xformed_origin = ray.origin - loc;

  //rtPrintf("%6.4E %6.4E %6.4E %6.4E %6.4E %6.4E %6.4E %6.4E %6.4E\n",xformed_origin.x,xformed_origin.y,xformed_origin.z,mins.x,mins.y,mins.z,maxs.x,maxs.y,maxs.z);

  float3 t0 = (mins - xformed_origin)/ray.direction;
  float3 t1 = (maxs - xformed_origin)/ray.direction;
  float3 near = fminf(t0, t1);
  float3 far  = fmaxf(t0, t1);
  float tmin  = fmaxf( near );
  float tmax  = fminf( far );
  float sgn   = 0;

    // compute sign, if true, points should be outside.  if one product is positive, point lies on same side of two parallel planes, meaning it must be outside the box
  if ( (t0.x*t1.x > 0) | (t0.y*t1.y > 0) | (t0.z*t1.z > 0) ){
    sgn =  1.0;
  }
  else{
    sgn = -1.0;
  }

  // report intersection
  if(tmin <= tmax) {
    bool check_second = true;
    if( rtPotentialIntersection( tmin ) ) {
        cellnum     = dims[object_dex].cellnum;
        cellmat     = dims[object_dex].matnum;
        cellfissile = dims[object_dex].is_fissile;
        normal      = boxnormal( tmin , t0 , t1 );
        normal      = sgn*normal;
        sense       = int(sgn);
        if(rtReportIntersection(0))
         check_second = false;
    } 
    if(check_second) {
      if( rtPotentialIntersection( tmax ) ) {
         cellnum     = dims[object_dex].cellnum;
         cellmat     = dims[object_dex].matnum;
         cellfissile = dims[object_dex].is_fissile;
         normal      = boxnormal( tmax , t0 , t1 );
         normal      = sgn*normal;
         sense       = int(sgn);
         rtReportIntersection(0);
      }
    }
  }
}

RT_PROGRAM void bounds (int object_dex, float result[6])
{
  float3 mins = make_float3(dims[object_dex].min[0],dims[object_dex].min[1],dims[object_dex].min[2]);
  float3 maxs = make_float3(dims[object_dex].max[0],dims[object_dex].max[1],dims[object_dex].max[2]);
  float3 loc  = make_float3(dims[object_dex].loc[0],dims[object_dex].loc[1],dims[object_dex].loc[2]);

  optix::Aabb* aabb = (optix::Aabb*)result;
  aabb->set(mins+loc, maxs+loc);
}
