#ifndef PRIMITIVE_H
#define PRIMITIVE_H

class primitive
{	
	public:
	 primitive();
	 primitive(int,unsigned,unsigned,float,float,float,float,float,float,float,float,float);
	~primitive();
	void add_transform();
	void add_transform(unsigned,float,float,float,float,float);  //defaults to primitive material
	void add_transform(unsigned,unsigned,float,float,float,float,float); //set own material
	void print_transform();
	void print_transform(int);
	void make_hex_array(int,float,float,float,unsigned);
	void make_hex_array(int,float,float,float,float,unsigned);
	float       min[3];
	float       max[3];
	float       location[3];
	static int  num_primitives;
	int			type;      // 0=box, 1=cyl, 2=hex
	int 		primitive_id;
	int         n_transforms;
	int         material;
	std::vector<wtransform>   transforms;
};

#endif
