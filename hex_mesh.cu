#include <optix.h>
#include <optixu/optixu_math_namespace.h>
#include <optixu/optixu_matrix_namespace.h>
#include <optixu/optixu_aabb_namespace.h>
#include "datadef.h"

using namespace optix;

rtBuffer<geom_data,1>               dims;
rtDeclareVariable(optix::Ray, ray, rtCurrentRay, );
rtDeclareVariable(unsigned,  cellnum,     attribute cell_num, );
rtDeclareVariable(int,       celltal,     attribute cell_tal, );
rtDeclareVariable(unsigned,  cellmat,     attribute cell_mat, );
rtDeclareVariable(unsigned,  cellfissile, attribute cell_fis, );
rtDeclareVariable(unsigned,  sense      , attribute cell_sense, );
rtDeclareVariable(float3,    normal,      attribute normal,   );
rtDeclareVariable(uint, launch_index_in, rtLaunchIndex, );

static __device__ bool accept_point(float3 pnt, float a, float x1, float x2, float zmin, float zmax)
{

  // point accepted if within the hex region
  float x = fabsf(pnt.x);
  float y = fabsf(pnt.y);
  float line = -a*(x-x2)/(x2-x1); 
  float tol = 1e-7;

  if( (pnt.z==zmax|pnt.z==zmin) & x > x1 & y > (line+fabsf(tol*line))){
    return false;
  }
  else{
    return true;
  }

}

static __device__ float get_t(float3 hat, float3 dir, float3 diff_points){

  float ndD;

  ndD = dot(hat,dir);

  if (ndD!=0.0){
    return dot(hat,diff_points) / ndD;
  }
  else{
    return 1e99;   // something big just to put it out of range
  }

}


RT_PROGRAM void intersect(int object_dex)
{

  // maxs.x is apothem
  float3 mins = make_float3(dims[object_dex].min[0],dims[object_dex].min[1],dims[object_dex].min[2]);
  float3 maxs = make_float3(dims[object_dex].max[0],dims[object_dex].max[1],dims[object_dex].max[2]);
  float3 loc  = make_float3(dims[object_dex].loc[0],dims[object_dex].loc[1],dims[object_dex].loc[2]);
  float3 xformed_origin = ray.origin - loc;

  // init
  float  t0=-1e34, t1=1e34, sgn=1.0, this_t, new_t=0.0, closest_xy_t=0.0, int_c_z=0.0;
  int    d0=0, d1=0;
  float3 norms[8], pts[8];
  bool report=true, check_second=true;

  // box/line region delimiters
  float x1 = maxs.x/sqrtf(3.0);
  float x2 = 2.0*x1;
  
  // normal vectors
  norms[0] = make_float3( 0.0            , 0.0     , 1.0 );  //  z
  norms[1] = make_float3( 0.0            , 0.0     ,-1.0 );  // -z
  norms[2] = make_float3( 0.0            , 1.0     , 0.0 );  //  y
  norms[3] = make_float3( 0.0            ,-1.0     , 0.0 );  // -y
  norms[4] = make_float3( sqrtf(3.0)/2.0 , 1.0/2.0 , 0.0 );  //  r
  norms[5] = make_float3(-sqrtf(3.0)/2.0 ,-1.0/2.0 , 0.0 );  // -r
  norms[6] = make_float3(-sqrtf(3.0)/2.0 , 1.0/2.0 , 0.0 );  //  l
  norms[7] = make_float3( sqrtf(3.0)/2.0 ,-1.0/2.0 , 0.0 );  // -l

  // points that define positions of all planes
  pts[0] = make_float3(  x2,  0.0         , maxs.z );   // 4
  pts[1] = make_float3(  x2,  0.0         , mins.z );   // 1
  pts[2] = make_float3( -x1,  maxs.x      , maxs.z );   // 2
  pts[3] = make_float3( -x1, -maxs.x      , maxs.z );   // 3
  pts[4] = make_float3(  x2,  0.0         , maxs.z );   // 1
  pts[5] = make_float3( -x1, -maxs.x      , maxs.z );   // 3
  pts[6] = make_float3( -x1,  maxs.x      , maxs.z );   // 2
  pts[7] = make_float3(  x2,  0.0         , maxs.z );   // 1

  // calculate t for closest point in xy
  if ( (xformed_origin.x*ray.direction.x + xformed_origin.y*ray.direction.y) == 0.0 | (ray.direction.x*ray.direction.x + ray.direction.y*ray.direction.y) == 0.0 ){ // z-only, use z
      closest_xy_t = -dot(xformed_origin,ray.direction) / dot(ray.direction,ray.direction);
  }
  else{
      closest_xy_t = -(xformed_origin.x*ray.direction.x + xformed_origin.y*ray.direction.y) / (ray.direction.x*ray.direction.x + ray.direction.y*ray.direction.y);
  }

  // calculate z point of the closest x-y point
  int_c_z = ray.direction.z * closest_xy_t + xformed_origin.z;

  // correct if the closest intersection point is outside z bounds
  if (int_c_z>maxs.z){
    closest_xy_t = (1.0 + copysignf(1e-6,-ray.direction.z))*get_t(norms[0],ray.direction,(pts[0]-xformed_origin)); // make t value slightly inside that of the top z plane
  }
  else if (int_c_z<mins.z){
    closest_xy_t = (1.0 + copysignf(1e-6,ray.direction.z))*get_t(norms[1],ray.direction,(pts[1]-xformed_origin));  // make t value slightly inside that of the bottom z plane
  }

  // get two xy points that are closest to origin
  // t1 is nearest positive w.r.t. closest_xy_t, t0 is nearest negative
  for (int i=0; i<8; i++) {
    // calculate intersection t value
    this_t = get_t(norms[i],ray.direction,(pts[i]-xformed_origin));
    // find the points to closest_xy_t
    new_t = this_t - closest_xy_t;
    if( new_t < 0.0) { 
      if(new_t > t0){
        t0=new_t;
        d0=i;
      }
    }
    else{ 
      if(new_t < t1){
        t1=new_t;
        d1=i;
      }
    }
  }

  // shift back to normal t
  t0 = t0 + closest_xy_t;
  t1 = t1 + closest_xy_t;

  // check these points for corner misses
  report =          accept_point( (xformed_origin + t0*ray.direction) , maxs.x, x1, x2, mins.z, maxs.z);
  report = report & accept_point( (xformed_origin + t1*ray.direction) , maxs.x, x1, x2, mins.z, maxs.z);
  if (!report){
    rtPrintf("CORNER MISS \no=numpy.array([% 10.8E,% 10.8E,% 10.8E])\ndir=numpy.array([% 10.8E,% 10.8E,% 10.8E])\n",xformed_origin.x,xformed_origin.y,xformed_origin.z,ray.direction.x,ray.direction.y,ray.direction.z);
  }

  // sense
  if (t0*t1 < 0.0){ // neg means inside
    sgn = -1.0;
  }
  else{
    sgn =  1.0;
  }

  // report intersection
  if(report) {
    if( rtPotentialIntersection( t0 ) ) {
        cellnum     = dims[object_dex].cellnum;
        celltal     = dims[object_dex].talnum;
        cellmat     = dims[object_dex].matnum;
        cellfissile = dims[object_dex].is_fissile;
        normal      = sgn*norms[d0];
        sense       = int(sgn);
       if(rtReportIntersection(0))
         check_second = false;
    }
    if(check_second) {
      if( rtPotentialIntersection( t1 ) ) {
         cellnum     = dims[object_dex].cellnum;
         celltal     = dims[object_dex].talnum;
         cellmat     = dims[object_dex].matnum;
         cellfissile = dims[object_dex].is_fissile;
         normal      = sgn*norms[d1];
         sense       = int(sgn);
        rtReportIntersection(0);
      }
    }
  }

}

RT_PROGRAM void bounds (int object_dex, float result[6])
{
  // max.x is apothem
  float3 mins = make_float3(dims[object_dex].min[0],dims[object_dex].min[1],dims[object_dex].min[2]);
  float3 maxs = make_float3(dims[object_dex].max[0],dims[object_dex].max[1],dims[object_dex].max[2]);
  float3 loc  = make_float3(dims[object_dex].loc[0],dims[object_dex].loc[1],dims[object_dex].loc[2]);

  result[0] = -2.0*maxs.x/sqrt(3.0)     + loc.x;
  result[1] =     -maxs.x               + loc.y;
  result[2] =      mins.z               + loc.z;
  result[3] =  2.0*maxs.x/sqrt(3.0)     + loc.x;
  result[4] =      maxs.x               + loc.y;
  result[5] =      maxs.z               + loc.z;
}
