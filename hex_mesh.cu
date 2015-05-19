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
rtDeclareVariable(float3,    normal,      attribute normal,   );

static __device__ bool accept_z(float3 pnt, float a, float x1, float x2, float zmin, float zmax)
{

  float x = fabsf(pnt.x);
  float y = fabsf(pnt.y);
  float line = -a*(x-x2)/(x2-x1); 

  if(x <= x1){

    if (y <= a){
      return true;
    }
    else{
      return false;
    }

  }
  else if(x <= x2){

    if (y <= line){
      return true;
    }
    else{
      return false;
    }

  }
  else{
    return false;
  }

}

static __device__ bool accept_y(float3 pnt, float a, float x1, float x2, float zmin, float zmax)
{

  float x = fabsf(pnt.x);
  float z =       pnt.z; 

  if(z <= zmax && z >= zmin){

    if(x <= x1){
      return true;
    }
    else{
      return false;
    }

  }
  else{
    return false;
  }

}

static __device__ bool accept_line(float3 pnt, float a, float x1, float x2, float zmin, float zmax)
{

  float x = fabsf(pnt.x);
  float z =       pnt.z;

  if(z <= zmax && z >= zmin){

    if(x >= x1 && x <= x2){
      return true;
    }
    else{
      return false;
    }
    
  }
  else{
    return false;
  }

}


RT_PROGRAM void intersect(int object_dex)
{

  // max.x is apothem
  float3 mins = make_float3(dims[object_dex].min[0],dims[object_dex].min[1],dims[object_dex].min[2]);
  float3 maxs = make_float3(dims[object_dex].max[0],dims[object_dex].max[1],dims[object_dex].max[2]);
  float3 loc  = make_float3(dims[object_dex].loc[0],dims[object_dex].loc[1],dims[object_dex].loc[2]);
  float3 xformed_origin = ray.origin - loc;

  float   tmin, tmax, t0, t1, ndD, sense, sgn;
  float3  norm_max, norm_min;

  //init
  tmin       =  9999999999999999999.0;
  tmax       = -9999999999999999999.0;
  sense      = -1.0;
  norm_max.x = 0.0;
  norm_max.y = 0.0;
  norm_max.z = 0.0;
  norm_min.x = 0.0;
  norm_min.y = 0.0;
  norm_min.z = 0.0;

  // box/line region delimiters
  float x1 = maxs.x/sqrtf(3.0);
  float x2 = 2.0*x1;
  
  // normal vectors
  float3 z_hat = make_float3( 0.0            , 0.0     , 1.0 );
  float3 y_hat = make_float3( 0.0            , 1.0     , 0.0 );
  float3 r_hat = make_float3( sqrtf(3.0)/2.0 , 1.0/2.0 , 0.0 );
  float3 l_hat = make_float3(-sqrtf(3.0)/2.0 , 1.0/2.0 , 0.0 );

  // points that define all planes
  float3 p4 = make_float3(  x2,  0.0         , mins.z );
  float3 p1 = make_float3(  x2,  0.0         , maxs.z );
  float3 p2 = make_float3( -x1,  maxs.x      , maxs.z );
  float3 p3 = make_float3( -x1, -maxs.x      , maxs.z );

  // z plane intersections
  ndD=dot(z_hat,ray.direction);
  if (ndD!=0.0){
    t0 = dot(  z_hat , ( p1 - xformed_origin ) ) / ndD;
    t1 = dot(  z_hat , ( p4 - xformed_origin ) ) / ndD;
    if (sense < 0.0){
      sense = t0*t1;
    }
    if(t0>0.0 && accept_z(xformed_origin+t0*ray.direction,maxs.x,x1,x2,mins.z,maxs.z)){
      tmin = fminf(t0,tmin);
      tmax = fmaxf(t0,tmax);
      if (tmin == t0){
        norm_min = z_hat;
      }
      if (tmax == t0){
        norm_max = z_hat;
      }
    }
    if(t1>0.0 && accept_z(xformed_origin+t1*ray.direction,maxs.x,x1,x2,mins.z,maxs.z)){
      tmin = fminf(t1,tmin);
      tmax = fmaxf(t1,tmax);
      if (tmin == t1){
        norm_min = -z_hat;
      }
      if (tmax == t1){
        norm_max = -z_hat;
      }
    }
  }

  // y line intersections
  ndD=dot(y_hat,ray.direction);
  if (ndD!=0.0){
    t0 = dot(  y_hat , ( p2 - xformed_origin ) ) / ndD;
    t1 = dot(  y_hat , ( p3 - xformed_origin ) ) / ndD;
    if (sense < 0.0){
      sense = t0*t1;
    }
    if(t0>0.0 && accept_y(xformed_origin+t0*ray.direction,maxs.x,x1,x2,mins.z,maxs.z)){
      tmin = fminf(t0,tmin);
      tmax = fmaxf(t0,tmax);
      if (tmin == t0){
        norm_min = y_hat;
      }
      if (tmax == t0){
        norm_max = y_hat;
      }
    }
    if(t1>0.0 && accept_y(xformed_origin+t1*ray.direction,maxs.x,x1,x2,mins.z,maxs.z)){
      tmin = fminf(t1,tmin);
      tmax = fmaxf(t1,tmax);
      if (tmin == t1){
        norm_min = -y_hat;
      }
      if (tmax == t1){
        norm_max = -y_hat;
      }
    }
  }

  // left line intersections
  ndD=dot(l_hat,ray.direction);
  if (ndD!=0.0){
    t0 = dot(  l_hat , ( p2 - xformed_origin ) ) / ndD;
    t1 = dot(  l_hat , ( p1 - xformed_origin ) ) / ndD;
    if (sense < 0.0){
      sense = t0*t1;
    }
    if(t0>0.0 && accept_line(xformed_origin+t0*ray.direction,maxs.x,x1,x2,mins.z,maxs.z)){
      tmin = fminf(t0,tmin);
      tmax = fmaxf(t0,tmax);
      if (tmin == t0){
        norm_min = l_hat;
      }
      if (tmax == t0){
        norm_max = l_hat;
      }
    }
    if(t1>0.0 && accept_line(xformed_origin+t1*ray.direction,maxs.x,x1,x2,mins.z,maxs.z)){
      tmin = fminf(t1,tmin);
      tmax = fmaxf(t1,tmax);
      if (tmin == t1){
        norm_min = -l_hat;
      }
      if (tmax == t1){
        norm_max = -l_hat;
      }
    }
  }

  // right line intersections
  ndD=dot(r_hat,ray.direction);
  if (ndD!=0.0){
    t0 = dot(  r_hat , ( p1 - xformed_origin ) ) / ndD;
    t1 = dot(  r_hat , ( p3 - xformed_origin ) ) / ndD;
    if (sense < 0.0){
      sense = t0*t1;
    }
    if(t0>0.0 && accept_line(xformed_origin+t0*ray.direction,maxs.x,x1,x2,mins.z,maxs.z)){
      tmin = fminf(t0,tmin);
      tmax = fmaxf(t0,tmax);
      if (tmin == t0){
        norm_min = r_hat;
      }
      if (tmax == t0){
        norm_max = r_hat;
      }
    }
    if(t1>0.0 && accept_line(xformed_origin+t1*ray.direction,maxs.x,x1,x2,mins.z,maxs.z)){
      tmin = fminf(t1,tmin);
      tmax = fmaxf(t1,tmax);
      if (tmin == t1){
        norm_min = -r_hat;
      }
      if (tmax == t1){
        norm_max = -r_hat;
      }
    }
  }

  // if a single sense is positive, point lies on same side of two parallel planes and it is outside the hexagon
  if ( sense < 0.0){
    sgn =  1.0;   
  }
  else{
    sgn = -1.0;   // switch sign of normal to inward
  }

  // report intersection
  if(tmin <= tmax) {
    bool check_second = true;
    if( rtPotentialIntersection( tmin ) ) {
        cellnum     = dims[object_dex].cellnum;
        cellmat     = dims[object_dex].matnum;
        cellfissile = dims[object_dex].is_fissile;
        normal      = sgn*norm_min;
       if(rtReportIntersection(0))
         check_second = false;
    } 
    if(check_second) {
      if( rtPotentialIntersection( tmax ) ) {
         cellnum     = dims[object_dex].cellnum;
         cellmat     = dims[object_dex].matnum;
         cellfissile = dims[object_dex].is_fissile;
         normal      = sgn*norm_max;
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
