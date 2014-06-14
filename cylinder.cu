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
    float t1, t2, z1, z2, sdisc;

    float radius    = maxs.y;
    float zmin      = mins.x;
    float zmax      = maxs.x;

    float a =        ( ray.direction.x  * ray.direction.x  ) + ( ray.direction.y  * ray.direction.y  );
    float b = 2.0 * (( ray.direction.x  * ray.origin.x ) + ( ray.direction.y  * ray.origin.y ));
    float c =        ( ray.origin.x * ray.origin.x ) + ( ray.origin.y * ray.origin.y ) - (radius * radius);
    float disc = (b*b)-(4*a*c);

    bool report = false;
    bool check_second = true;

    if (disc > 0.0f){  // the line intersects the circle

        report = true;

        sdisc = sqrt(disc);
        t1 = (-b-sdisc)/(2*a);
        t2 = (-b+sdisc)/(2*a);

        z1 = ray.direction.z * t1 + ray.origin.z;
        z2 = ray.direction.z * t2 + ray.origin.z;

        //rtPrintf("zmin %6.4E zmax %6.4E t1 %6.4E t2 %6.4E z1 %6.4E z2 %6.4E report %u\n",zmin,zmax,t1,t2,z1,z2,report);

        if( ((z1 > zmax) & (z2 > zmax)) | ((z1 < zmin) & (z2 < zmin)) ){  //miss in corners
            report=false;
        }
        else{ 

            if (z1 > zmax ){  //  top intersection z1
                t1 = (zmax - ray.origin.z) / ray.direction.z;
            }
            else if(z1 < zmin ) { // bottom intersection z1
                t1 = (zmin - ray.origin.z) / ray.direction.z;
            }

            if (z2 > zmax){  //  top intersection z2
                t2 = (zmax - ray.origin.z) / ray.direction.z;
            }
            else if(z2 < zmin) { // bottom intersection z2
                t2 = (zmin - ray.origin.z) / ray.direction.z;
            }
        }
        
    }
    else if( (ray.origin.x*ray.origin.x+ray.origin.y*ray.origin.y)<(radius*radius) ) {  // exactly perpendicular

        report = true;

        t1 = fminf((zmax - ray.origin.z) / ray.direction.z , (zmin - ray.origin.z) / ray.direction.z);
        t2 = fmaxf((zmax - ray.origin.z) / ray.direction.z , (zmin - ray.origin.z) / ray.direction.z);  

    }

    //rtPrintf("zmin %6.4E zmax %6.4E t1 %6.4E t2 %6.4E z1 %6.4E z2 %6.4E report %u\n",zmin,zmax,t1,t2,z1,z2,report);

    if (report){
        if(t1>0){
            if (rtPotentialIntersection(t1) ) {
                if(rtReportIntersection(0)){
                    check_second=false;
                }
            }
        }
        if(check_second & t2>0){
            if (rtPotentialIntersection(t2) ) {
                rtReportIntersection(0);
            }
        }
    }

}

RT_PROGRAM void bounds (int, float result[6])
{
    float r    = maxs.y;
    float zmin = mins.x;
    float zmax = maxs.x;

    result[0] = -r;
    result[1] = -r;
    result[2] = zmin;
    result[3] = r;
    result[4] = r;
    result[5] = zmax;
}
