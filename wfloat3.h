#ifndef WFLOAT3_H
#define WFLOAT3_H

//class definitions for device vector operations
class wfloat3{
public:
	float x,y,z;
	inline __device__ 			wfloat3();
	inline __device__ 			wfloat3(float);
	inline __device__ 			wfloat3(float,float,float);
	inline __device__ wfloat3 	operator+ (wfloat3); 
	inline __device__ wfloat3 	operator- (wfloat3);
	inline __device__ wfloat3 	operator* (wfloat3);
	inline __device__ wfloat3 	operator+ (float); 
	inline __device__ wfloat3 	operator- (float);
	inline __device__ wfloat3 	operator* (float);
	inline __device__ wfloat3 	operator/ (float);
	inline __device__ wfloat3 	cross(wfloat3);
	inline __device__ float   	dot(wfloat3);
	inline __device__ void  	rodrigues_rotation(wfloat3,float);
	inline __device__ wfloat3  	rotate(float,float);
	inline __device__ float 	norm2();
};
__device__ wfloat3::wfloat3(){x=0;y=0;z=0;};
__device__ wfloat3::wfloat3(float a){x=a;y=a;z=a;};
__device__ wfloat3::wfloat3(float a,float b,float c){x=a;y=b;z=c;};
__device__ wfloat3 wfloat3::operator+ (wfloat3 arg){
	wfloat3 result(x+arg.x,y+arg.y,z+arg.z);
	return result;
}; 
__device__ wfloat3 wfloat3::operator- (wfloat3 arg){
	wfloat3 result(x-arg.x,y-arg.y,z-arg.z);
	return result;
};
__device__ wfloat3 wfloat3::operator* (wfloat3 arg){
	wfloat3 result(x*arg.x,y*arg.y,z*arg.z);
	return result;
};
__device__ wfloat3 wfloat3::operator+ (float arg){
	wfloat3 result(x+arg,y+arg,z+arg);
	return result;
}; 
__device__ wfloat3 wfloat3::operator- (float arg){
	wfloat3 result(x-arg,y-arg,z-arg);
	return result;
};
__device__ wfloat3 wfloat3::operator* (float arg){
	wfloat3 result(x*arg,y*arg,z*arg);
	return result;
};
__device__ wfloat3 wfloat3::operator/ (float arg){
	wfloat3 result(x/arg,y/arg,z/arg);
	return result;
};
__device__ wfloat3 wfloat3::cross(wfloat3 arg){
	wfloat3 result;
	result.x =  y*arg.z - arg.y*z ;
	result.y = -x*arg.z + arg.x*z ;
	result.z =  x*arg.y - arg.x*y ;
	return result;
};
__device__ float wfloat3::dot(wfloat3 arg){
	return (x*arg.x + y*arg.y + z*arg.z);
};
__device__ void wfloat3::rodrigues_rotation(wfloat3 k, float theta){
	*this = (*this)*cosf(theta) - (k.cross(*this))*sinf(theta) + k*(k.dot(*this))*(1.0-cosf(theta));
};
__device__ float wfloat3::norm2(){
	return sqrtf( x*x + y*y + z*z );
};
__device__ wfloat3 wfloat3::rotate(float mu, float rn){
	// borrowed from OpenMC
	wfloat3 out;
	float phi = 6.28318530718 * rn;
    float a = sqrtf(max(0.0, 1.0 - mu*mu));
    float b = sqrtf(max(0.0, 1.0 - this[0].z*this[0].z));
    float cosphi = cosf(phi);
    float sinphi = sinf(phi);
    // Need to treat special case where sqrt(1 - w**2) is close to zero by
    // expanding about the v component rather than the w component
    if (b > 1e-10) {
      out.x = mu*this[0].x + a*(this[0].x*this[0].z*cosphi - this[0].y*sinphi)/b;
      out.y = mu*this[0].y + a*(this[0].y*this[0].z*cosphi + this[0].x*sinphi)/b;
      out.z = mu*this[0].z - a*b*cosphi;
  	}
    else{
      b = sqrtf(1.0 - this[0].y*this[0].y);
      out.x = mu*this[0].x + a*(this[0].x*this[0].y*cosphi + this[0].z*sinphi)/b;
      out.y = mu*this[0].y - a*b*cosphi;
      out.z = mu*this[0].z + a*(this[0].y*this[0].z*cosphi - this[0].x*sinphi)/b;
    }
    return out;
};

#endif
