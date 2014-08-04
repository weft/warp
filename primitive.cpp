#include <vector> 
#include <iostream>
#include <sstream>
#include <cmath>
#include <assert.h>
#include <time.h>
#include "datadef.h"
#include "primitive.h"

int primitive::num_primitives=0;
primitive::primitive(){
	//box default constructor
	min[0]=0;min[1]=0;min[2]=0;
	max[0]=0;max[1]=0;max[2]=0;
	location[0]=0;location[1]=0;location[2]=0;
	type=0;
	material=0;
	primitive_id=num_primitives;
	num_primitives++;
	n_transforms=0;

}
primitive::primitive(int ptype, unsigned cellmat , std::vector<float> mins, std::vector<float> maxs, std::vector<float> locs){
	//box valued constructor
	min[0]=mins[0];
	min[1]=mins[1];
	min[2]=mins[2];
	max[0]=maxs[0];
	max[1]=maxs[1];
	max[2]=maxs[2];
	location[0]=locs[0];
	location[1]=locs[1];
	location[2]=locs[2];
	type=ptype;
	material=cellmat;
	primitive_id=num_primitives;
	num_primitives++;
	n_transforms=0;
}
primitive::~primitive(){
	//delete min;
	//delete max;
	//delete location;
	num_primitives--;
	//delete &type;
	//delete &primitive_id;
	//delete &n_transforms;
	//delete &material;
	//transforms.~vector();
}
	/**
*/
unsigned primitive::add_transform(){
	wtransform this_transform;
	if(transforms.empty()){
		this_transform.cellnum = 0;
	}
	else{
	 	this_transform = transforms.back();
	}
	this_transform.cellmat = material;
	this_transform.dx      = 0;
	this_transform.dy      = 0;
	this_transform.dz      = 0;
	this_transform.theta   = 0;
	this_transform.phi     = 0;
	transforms.push_back(this_transform);
	n_transforms=transforms.size();
	return (n_transforms-1);
}
/**
*/
unsigned primitive::add_transform(unsigned cellnum , float dx , float dy , float dz , float theta , float phi ){
	wtransform this_transform;
	this_transform.cellnum = cellnum;
	this_transform.cellmat = material;
	this_transform.dx      = dx;
	this_transform.dy      = dy;
	this_transform.dz      = dz;
	this_transform.theta   = theta;
	this_transform.phi     = phi;
	transforms.push_back(this_transform);
	n_transforms=transforms.size();
	return (n_transforms-1);
}
/**
*/
unsigned primitive::add_transform(unsigned cellnum ,unsigned cellmat, float dx , float dy , float dz , float theta , float phi ){
	wtransform this_transform;
	this_transform.cellnum = cellnum;
	this_transform.cellmat = cellmat;
	this_transform.dx      = dx;
	this_transform.dy      = dy;
	this_transform.dz      = dz;
	this_transform.theta   = theta;
	this_transform.phi     = phi;
	transforms.push_back(this_transform);
	n_transforms=transforms.size();
	return (n_transforms-1);
}
void primitive::print_transform(){
	std::cout << "--- primitive id  = " << primitive_id << " ---" << "\n";
	std::cout << "   min,max  = (" << min[0] << " , " << min[1] << " , " << min[2] << "),(" << max[0] << " , " << max[1] << " , " << max[2] << ")" << "\n";
	std::cout << "   location = (" << location[0] << " , " << location[1] << " , " << location[2] << ")" << "\n";
	std::cout << "   type =  " << type << "\n";
	std::cout << "   material = " << material << "\n";
	for (int tnum=0;tnum<n_transforms;tnum++){
		std::cout << "   ************ " << "\n";
		std::cout << "   transform= " << tnum << "\n";
		std::cout << "   cellnum  = " << transforms[tnum].cellnum << "\n";
		std::cout << "   cellmat  = " << transforms[tnum].cellmat << "\n";
		std::cout << "   dx       = " << transforms[tnum].dx << "\n";
		std::cout << "   dy       = " << transforms[tnum].dy << "\n";
		std::cout << "   dz       = " << transforms[tnum].dz << "\n";
		std::cout << "   theta    = " << transforms[tnum].theta << "\n";
		std::cout << "   phi      = " << transforms[tnum].phi << "\n";
	}
}
void primitive::print_transform(int tnum){
	std::cout << "--- primitive id  = " << primitive_id << " ---" << "\n";
	std::cout << "   min,max  = (" << min[0] << " , " << min[1] << " , " << min[2] << "),(" << max[0] << " , " << max[1] << " , " << max[2] << ")" << "\n";
	std::cout << "   location = (" << location[0] << " , " << location[1] << " , " << location[2] << ")" << "\n";
	std::cout << "   type =  " << type << "\n";
	std::cout << "   material = " << material << "\n";
	std::cout << "   ************ " << "\n";
	std::cout << "   transform= " << tnum << "\n";
	std::cout << "   cellnum  = " << transforms[tnum].cellnum << "\n";
	std::cout << "   cellmat  = " << transforms[tnum].cellmat << "\n";
	std::cout << "   dx       = " << transforms[tnum].dx << "\n";
	std::cout << "   dy       = " << transforms[tnum].dy << "\n";
	std::cout << "   dz       = " << transforms[tnum].dz << "\n";
	std::cout << "   theta    = " << transforms[tnum].theta << "\n";
	std::cout << "   phi      = " << transforms[tnum].phi << "\n";
}
void primitive::make_hex_array(int n, float x, float y, float PD_ratio, unsigned starting_index){

	wtransform this_transform;

	int k, j, num, cnt;
	float offsetx, offsety, fnum, lattr; //PD_ratio;

	// get strting cell number as the one set for the last
	//wtransform this_transform = transforms.back();
	//unsigned starting_index = this_transform.cellnum + 1;

	// add num of transforms to whatever is there
    n_transforms += 3*n*(n-1)+1;

	num=n;
	cnt=0;
	//PD_ratio=1.164;
	lattr = PD_ratio * (2.0 * max[1]) / sqrt(3.0);
	offsety=(n-1)*lattr*1.5;

	//row
	for (k=0;k<(2*n-1);k++){
		fnum=num-1;
		offsetx=-(sqrt(3.0)*lattr)*(fnum/2.0);
		//column
		for(j=0;j<num;j++){

    		this_transform.cellnum=starting_index+cnt;
    		this_transform.cellmat=material;
    		this_transform.dx=offsetx;
    		this_transform.dy=offsety;
    		this_transform.dz=0;
    		this_transform.theta=0;
    		this_transform.phi=0;
    		transforms.push_back(this_transform);
    		cnt++;
    		offsetx+=sqrt(3.0)*lattr;
		}

		if ( k < n-1 ){
			num++;   }
		else{
		 	num--;   }

		offsety-=lattr*1.5;

	}

}
