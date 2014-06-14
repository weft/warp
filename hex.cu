#include <optix.h>
#include <optixu/optixu_math_namespace.h>
#include <optixu/optixu_matrix_namespace.h>
#include <optixu/optixu_aabb_namespace.h>

using namespace optix;

rtDeclareVariable(float3, mins, , );
rtDeclareVariable(float3, maxs, , );
rtDeclareVariable(optix::Ray, ray, rtCurrentRay, );

RT_PROGRAM void intersect(int)
{
  float   tvec[8], tmin, ndD;
  float   ax, ay, az;
  float3  pvec[8];
  int     k;
  bool    report;

  // box/line region delimiters
  float x1 = maxs.y/sqrt(3.0);
  float x2 = 2*x1;
  
  // normal vectors
  float3 z_hat = make_float3( 0.0 , 0.0         , 1.0);
  float3 y_hat = make_float3( 0.0 , 1.0         , 0.0 );
  float3 r_hat = make_float3( 1.0 , 1/sqrt(3.0) , 0.0 );
  float3 l_hat = make_float3(-1.0 , 1/sqrt(3.0) , 0.0 );

  // points that define all planes
  float3 p4 = make_float3(  x2,  0           , mins.x );
  float3 p1 = make_float3(  x2,  0           , maxs.x );
  float3 p2 = make_float3( -x1,  maxs.y , maxs.x );
  float3 p3 = make_float3( -x1, -maxs.y , maxs.x );

  // find all plane intersections
  ndD=dot(z_hat,ray.direction);
  if (ndD!=0.0){
    tvec[0]=dot(z_hat,(p1-ray.origin))/ndD;
    tvec[1]=dot(z_hat,(p4-ray.origin))/ndD;
    pvec[0] = ray.origin+tvec[0]*ray.direction;
    pvec[1] = ray.origin+tvec[1]*ray.direction;
  }
  ndD=dot(y_hat,ray.direction);
  if (ndD!=0.0){
    tvec[2]=dot(  y_hat , ( p2 - ray.origin ) ) / ndD;
    tvec[3]=dot(  y_hat , ( p3 - ray.origin ) ) / ndD;
    pvec[2] = ray.origin+tvec[2]*ray.direction;
    pvec[3] = ray.origin+tvec[3]*ray.direction;
  }
  ndD=dot(l_hat,ray.direction);
  if (ndD!=0.0){
    tvec[4]=dot(  l_hat , ( p2 - ray.origin ) ) / ndD;
    tvec[5]=dot(  l_hat , ( p1 - ray.origin ) ) / ndD;
    pvec[4] = ray.origin+tvec[4]*ray.direction;
    pvec[5] = ray.origin+tvec[5]*ray.direction;
  }
  ndD=dot(r_hat,ray.direction);
  if (ndD!=0.0){
    tvec[6]=dot(  r_hat , ( p1 - ray.origin ) ) / ndD;
    tvec[7]=dot(  r_hat , ( p3 - ray.origin ) ) / ndD;
    pvec[6] = ray.origin+tvec[6]*ray.direction;
    pvec[7] = ray.origin+tvec[7]*ray.direction;
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
                }
            }
        }
        // is in line region
        else if (ax>x1 && ax<=x2 && (ay-(maxs.y*(2-ax/x1)))<=1e-6){
            if (az>=mins.x && az<=maxs.x) {
                if (tvec[k] >= 1e-8 && tvec[k] <= tmin) {
                    tmin=tvec[k];
                    report=true;
                }
            }
        }
    }

  // report t value of first intersection
  if(report) {
    if( rtPotentialIntersection( tmin ) ) {
       rtReportIntersection(0);
    }
  }

}

RT_PROGRAM void bounds (int, float result[6])
{
  result[0] = -2*maxs.y/sqrt(3.0);
  result[1] = -maxs.y;
  result[2] = mins.x;
  result[3] = 2*maxs.y/sqrt(3.0);
  result[4] = maxs.y;
  result[5] = maxs.x;
}
