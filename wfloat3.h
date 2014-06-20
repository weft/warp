#ifndef WFLOAT3_H
#define WFLOAT3_H

/**
 * \class wfloat3 wfloat3.h
 * \brief class definitions for device vector operations
 */

class wfloat3{
public:
	float x; /**< x-coordinate */
	float y; /**< y-coordinate */
	float z; /**< z-coordinate */
	/**
	 * \brief sets x,y,z to 0,0,0
	 */
	inline __device__ 		wfloat3();
	/**
	 * \brief sets x,y,z to a,a,a
	 * @param[in] a - point to set
	 */ 
	inline __device__ 		wfloat3(float);
	/**
	 * \brief sets x,y,z to a,b,c
	 * @param[in] a,b,c - points to set
	 */
	inline __device__ 		wfloat3(float,float,float);
	/**
	 * \brief vector addition operator
	 * \details adds x and x-component of input wfloat3, etc.
	 * @param[in] arg - wfloat3 coordinates to be added
	 * \returns result - resultant wfloat3
	 */
	inline __device__ wfloat3 	operator+ (wfloat3); 
	/**
	 * \brief vector subtraction operator
	 * \details subtracts x-component of input wfloat3 from x, etc.
	 * @param[in] arg - wfloat3 coordinates to be subtracted
	 * \returns result - resultant wfloat3
	 */
	inline __device__ wfloat3 	operator- (wfloat3);
	/**
	 * \brief vector multiplication operator
	 * \details multiplies x and x-component of input wfloat3, etc.
	 * @param[in] arg - wfloat3 coordinates to be multiplied
	 * \returns result - resultant wfloat3
	 */
	inline __device__ wfloat3 	operator* (wfloat3);
	/**
	 * \brief scalar addition operator
	 * \details adds x and arg, etc.
	 * @param[in] arg - number to be added
	 * \returns result - resultant wfloat3
	 */
	inline __device__ wfloat3 	operator+ (float); 
	/**
	 * \brief scalar subtraction operator
	 * \details subtracts arg from x, etc.
	 * @param[in] arg - number to be subtracted
	 * \returns result - resultant wfloat3
	 */
	inline __device__ wfloat3 	operator- (float);
	/**
	 * \brief scalar multiplication operator
	 * \details multiplies x and arg, etc.
	 * @param[in] arg - number by which to multiply
	 * \returns result - resultant wfloat3
	 */
	inline __device__ wfloat3 	operator* (float);
	/**
	 * \brief scalar divison operator
	 * \details divides x by arg, etc.
	 * @param[in] arg - number by which to divide
	 * \returns result - resultant wfloat3
	 */
	inline __device__ wfloat3 	operator/ (float);
	/**
	 * \brief cross product operator
	 * \details returns the cross product of the vector and arg
	 * @param[in] arg - vector to cross
	 * \returns result - resultant wfloat3
	 */
	inline __device__ wfloat3 	cross(wfloat3);
	/**
	 * \brief dot product operator
	 * \details returns the dot product of the vector and arg
	 * @param[in] arg - vector to dot
	 */
	inline __device__ float   	dot(wfloat3);
	/**
	 * \brief Rodrigues' rotation operator
	 * \details rotates a vector in space, given axis and angle of rotation
	 * @param[in]  k - unit vector describing axis of rotation about which to rotate
	 * @param[in] theta - angle by which to rotate
	 */
	inline __device__ void  	rodrigues_rotation(wfloat3,float);
	/**
	 * \brief rotation about random cosine
	 * \details borrowed from OpenMC
	 * @param[in] mu - random cos(theta)
	 * @param[in] rn - random number
	 */ 
	inline __device__ wfloat3  	rotate(float,float);
	/**
	 * \brief returns square root of sum of squares of coordinates
	 */
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
