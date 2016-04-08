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
  float tol = 0.0;//1e-6; 

  // check z
  if( pnt.z > (zmax+fabsf(tol*zmax)) | pnt.z < (zmin-fabsf(tol*zmin)) ){
    return false;
  }

  // check xy
  if (     x > (x2+fabsf(tol*x2)) ){
    return false;
  }
  else if( x > x1){
    if(    y <= (line+fabsf(tol*line)) ){return true;}
    else{                                return false;}
  }
  else{
    if(    y <= (a+fabsf(tol*a))       ){return true;}
    else{                                return false;}
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
  //float  max_diff = 2.0*sqrtf(2.0*maxs.x*maxs.x+maxs.z*maxs.z); // BB chord, maxium difference possible in t values
  float  t0=1e34, t1=1e34, sgn=1.0, this_t, t[8];
  int    d[8];
  float3  this_int, norm0, norm1, norms[8], pts[8];
  bool report=true, check_second=true;

  // box/line region delimiters
  float x1 = maxs.x/sqrtf(3.0);
  float x2 = 2.0*x1;

  // init
  t[0] = 0.0;
  t[1] = 0.0;
  t[2] = 0.0;
  t[3] = 0.0;
  t[4] = 0.0;
  t[5] = 0.0;
  t[6] = 0.0;
  t[7] = 0.0;
  d[0] = 0;
  d[1] = 0;
  d[2] = 0;
  d[3] = 0;
  d[4] = 0;
  d[5] = 0;
  d[6] = 0;
  d[7] = 0;
  
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

  //  do xy-perpendicular planes first
  int k=0;
  for (int i=0; i<8; i++) {
    // calculate intersection t value
    this_t = get_t(norms[i],ray.direction,(pts[i]-xformed_origin));
    // calculate intersection point from t value
    this_int = ray.direction * this_t + xformed_origin;
    // accept if within maximum radius and within z values
    if( accept_point(this_int, maxs.x, x1, x2, mins.z, maxs.z) ){
      t[k]=this_t;
      d[k]=i;
      k++;
    }
    // break if a double is flagged
    if(k==3){break;}
  }

  // now find any missing points or determine if its a corner miss
  if(k==0){
    report = false;
  }
  else if(k==1){
    // not good
    t0=t[0];
    norm0=norms[d[0]];
    report = true;
    check_second = false;
    rtPrintf("k==1! t(%u,:)=[%10.8E, %10.8E, %10.8E, %10.8E, %10.8E, %10.8E, %10.8E, %10.8E];\n",launch_index_in,t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7]);
    rtPrintf("k==1! o(%u,:)=[%10.8E, %10.8E, %10.8E];dir(%u,:)=[%10.8E, %10.8E, %10.8E];\n",launch_index_in,xformed_origin.x,xformed_origin.y,xformed_origin.z,tid,ray.direction.x,ray.direction.y,ray.direction.z);
  }
  else if(k==2){
    // good
    t0=t[0];
    t1=t[1];
    norm0=norms[d[0]];
    norm1=norms[d[1]];
    report = true;
    check_second = true;
    //rtPrintf("k==2! t=[%10.8E, %10.8E];o=[%10.8E, %10.8E, %10.8E];d=[%10.8E, %10.8E, %10.8E];\n",t0,t1,xformed_origin.x,xformed_origin.y,xformed_origin.z,ray.direction.x,ray.direction.y,ray.direction.z);
  }
  else{
    // also not good
    t0=t[0];
    t1=t[1];
    norm0=norms[d[0]];
    norm1=norms[d[1]];
    report = true;
    check_second = true;
    rtPrintf("k==%d! t(%u,:)=[%10.8E, %10.8E, %10.8E, %10.8E, %10.8E, %10.8E, %10.8E, %10.8E];\n",k,launch_index_in,t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7]);
    rtPrintf("k==%d! o(%u,:)=[%10.8E, %10.8E, %10.8E];dir(%u,:)=[%10.8E, %10.8E, %10.8E];\n",k,launch_index_in,xformed_origin.x,xformed_origin.y,xformed_origin.z,tid,ray.direction.x,ray.direction.y,ray.direction.z);
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
        normal      = sgn*norm0;
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
         normal      = sgn*norm1;
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
