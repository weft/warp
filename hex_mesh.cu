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
rtDeclareVariable(float3, normal, attribute normal, );


RT_PROGRAM void intersect(int object_dex)
{

  float3 mins = make_float3(dims[object_dex].min[0],dims[object_dex].min[1],dims[object_dex].min[2]);
  float3 maxs = make_float3(dims[object_dex].max[0],dims[object_dex].max[1],dims[object_dex].max[2]);
  float3 loc  = make_float3(dims[object_dex].loc[0],dims[object_dex].loc[1],dims[object_dex].loc[2]);
  float3 xformed_origin = ray.origin - loc;

  float   tvec[8], tmin, ndD;
  float   ax, ay, az, sgn;
  float3  pvec[8];
  float3  norms[8];
  int     k;
  bool    report;

  // box/line region delimiters
  float x1 = maxs.y/sqrt(3.0);
  float x2 = 2*x1;
  
  // normal vectors
  float3 this_norm;
  float3 z_hat = make_float3( 0.0 , 0.0         , 1.0);
  float3 y_hat = make_float3( 0.0 , 1.0         , 0.0 );
  float3 r_hat = make_float3( sqrt(3.0)/2 , 1/2 , 0.0 );
  float3 l_hat = make_float3(-sqrt(3.0)/2 , 1/2 , 0.0 );

  // points that define all planes
  float3 p4 = make_float3(  x2,  0           , mins.x );
  float3 p1 = make_float3(  x2,  0           , maxs.x );
  float3 p2 = make_float3( -x1,  maxs.y      , maxs.x );
  float3 p3 = make_float3( -x1, -maxs.y      , maxs.x );

  // find all plane intersections
  ndD=dot(z_hat,ray.direction);
  if (ndD!=0.0){
    tvec[0]=dot(z_hat,(p1-xformed_origin))/ndD;
    tvec[1]=dot(z_hat,(p4-xformed_origin))/ndD;
    pvec[0] = xformed_origin+tvec[0]*ray.direction;
    pvec[1] = xformed_origin+tvec[1]*ray.direction;
    norms[0]= -z_hat;
    norms[1]=  z_hat;
  }
  ndD=dot(y_hat,ray.direction);
  if (ndD!=0.0){
    tvec[2]=dot(  y_hat , ( p2 - xformed_origin ) ) / ndD;
    tvec[3]=dot(  y_hat , ( p3 - xformed_origin ) ) / ndD;
    pvec[2] = xformed_origin+tvec[2]*ray.direction;
    pvec[3] = xformed_origin+tvec[3]*ray.direction;
    norms[2]=  y_hat;
    norms[3]= -y_hat;
  }
  ndD=dot(l_hat,ray.direction);
  if (ndD!=0.0){
    tvec[4]=dot(  l_hat , ( p2 - xformed_origin ) ) / ndD;
    tvec[5]=dot(  l_hat , ( p1 - xformed_origin ) ) / ndD;
    pvec[4] = xformed_origin+tvec[4]*ray.direction;
    pvec[5] = xformed_origin+tvec[5]*ray.direction;
    norms[4]=  l_hat;
    norms[5]= -l_hat;
  }
  ndD=dot(r_hat,ray.direction);
  if (ndD!=0.0){
    tvec[6]=dot(  r_hat , ( p1 - xformed_origin ) ) / ndD;
    tvec[7]=dot(  r_hat , ( p3 - xformed_origin ) ) / ndD;
    pvec[6] = xformed_origin+tvec[6]*ray.direction;
    pvec[7] = xformed_origin+tvec[7]*ray.direction;
    norms[6]=  r_hat;
    norms[7]= -r_hat;
  }

  // compute sign, if true, points should be outside.  if product is positive, point lies on same side of two parallel planes
  if ( (tvec[0]*tvec[1] > 0) | (tvec[2]*tvec[3] > 0) | (tvec[4]*tvec[5] > 0) | (tvec[6]*tvec[7] > 0) ){
    sgn = -1.0;   // switch sign of normal to inward
  }
  else{
    sgn =  1.0;
  }

  // get hits that are in-bounds (should only be 2, report one with smallest t)
    report=false;
    tmin=1.0/0.0;
    for (k=0;k<8;k++) {
        ax = fabsf(pvec[k].x);
        ay = fabsf(pvec[k].y);
        az = fabsf(pvec[k].z);
        // is in box region
        if (ax<=x1 && ay<=maxs.y) {
            if (az>=mins.x && az<=maxs.x) {
                if (tvec[k] >= 1e-8 && tvec[k] <= tmin) {
                    tmin=tvec[k];
                    report=true;
                    this_norm = sgn*norms[k];
                }
            }
        }
        // is in line region
        else if (ax>x1 && ax<=x2 && (ay-(maxs.y*(2-ax/x1)))<=1e-6){
            if (az>=mins.x && az<=maxs.x) {
                if (tvec[k] >= 1e-8 && tvec[k] <= tmin) {
                    tmin=tvec[k];
                    this_norm = sgn*norms[k];
                    report=true;
                }
            }
        }
    }

  // report t value of first intersection
  if(report) {
    if( rtPotentialIntersection( tmin ) ) {
      cellnum     = dims[object_dex].cellnum;
      cellmat     = dims[object_dex].matnum;
      cellfissile = dims[object_dex].is_fissile;
      normal      = this_norm;
      rtReportIntersection(0);
    }
  }

}

RT_PROGRAM void bounds (int object_dex, float result[6])
{
  float3 mins = make_float3(dims[object_dex].min[0],dims[object_dex].min[1],dims[object_dex].min[2]);
  float3 maxs = make_float3(dims[object_dex].max[0],dims[object_dex].max[1],dims[object_dex].max[2]);
  float3 loc  = make_float3(dims[object_dex].loc[0],dims[object_dex].loc[1],dims[object_dex].loc[2]);

  result[0] = -2*maxs.y/sqrt(3.0) + loc.x;
  result[1] =   -maxs.y           + loc.y;
  result[2] =    mins.x           + loc.z;
  result[3] =  2*maxs.y/sqrt(3.0) + loc.x;
  result[4] =    maxs.y           + loc.y;
  result[5] =    maxs.x           + loc.z;
}
