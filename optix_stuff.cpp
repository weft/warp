#include <vector> 
#include <iostream>
#include <sstream>
#include <cmath>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <png++/png.hpp>
#include "datadef.h"
#include "primitive.h"
#include "wgeometry.h"
#include "optix_stuff.h"
#include "device_copies.h"

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//							OptiX stuff
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////


optix_stuff::optix_stuff(){}
optix_stuff::optix_stuff(unsigned Nin,unsigned mult){
	//set stack size multiplier
	stack_size_multiplier = mult;
	//set main N
	N=Nin;
}
optix_stuff::~optix_stuff(){
	//try {
	//	context->destroy();	
	//} 
	//catch( optix::Exception &e ){
	//	std::cout << "OptiX Error on destroy: " << e.getErrorString().c_str() << "\n";
	//	exit(1);
	//}
}
void optix_stuff::init_internal(wgeometry problem_geom, unsigned compute_device_in, std::string accel_type_in){

	using namespace optix;

	// local variables
	char	 path_to_ptx[512];
	Program  ray_gen_program;
	Program  exception_program;  
	Program  miss_program;
	Buffer 	 positions_buffer;
	Buffer 	 rxn_buffer;
	Buffer 	 done_buffer;
	Buffer 	 cellnum_buffer;
	Buffer 	 matnum_buffer;
	Buffer 	 remap_buffer;
	Variable positions_var;
	Variable rxn_var;
	Variable done_var;
	Variable cellnum_var;
	Variable matnum_var;
	Variable remap_var;
	Variable outer_cell_var;
	Variable boundary_condition_var;
	Variable trace_type_var;
	RTsize   stack_size;
	RTsize	 printf_size;

	//set vars
	compute_device 	= compute_device_in;
	optix_device 	= 0;
	accel_type     	= accel_type_in;   
	     if(accel_type.compare("Sbvh")==0){traverse_type="Bvh";}
	else if(accel_type.compare("Bvh")==0) {traverse_type="Bvh";}
	else if(accel_type.compare("MedianBvh")==0) {traverse_type="Bvh";}
	else if(accel_type.compare("Lbvh")==0) {traverse_type="Bvh";}
	else if(accel_type.compare("TriangleKdTree")==0){traverse_type="KdTree";}

	// set geom type  0=primitive instancing, 1=transform instancing, 2=transform instancing with common primitives
	GEOM_FLAG = 0;

	//set image type
	image_type = "cell";

  	//int computeCaps[2];
  	//if (RTresult code = rtDeviceGetAttribute(deviceId, RT_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY, sizeof(computeCaps), &computeCaps))
  	//	throw Exception::makeException(code, 0);
  	//for(unsigned int index = 1; index < Context::getDeviceCount(); ++index) {
	//	int computeCapsB[2];
	//	if (RTresult code = rtDeviceGetAttribute(index, RT_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY, sizeof(computeCaps), &computeCapsB))
	//		throw Exception::makeException(code, 0);
	//	if (computeCaps[0] == computeCapsB[0] && computeCaps[1] == computeCapsB[1]) {
	//		deviceId = index;
	//		break;
	//	}
  	//}

	// create context
	context = Context::create();
	//get device info
	unsigned deviceId = compute_device;
	unsigned deviceCount = context -> getDeviceCount();
	assert(deviceId<deviceCount);
	context -> setDevices(&deviceId, &deviceId+1);  // iterator_start, iterator_end
	unsigned enabled_count = context -> getEnabledDeviceCount();
	std::vector<int> enabled_ids = context -> getEnabledDevices();
	printf("OptiX using device ");
	for(unsigned h=0;h<enabled_count;h++){printf("%u ",enabled_ids[h]); optix_device=enabled_ids[h];}
	printf("\n");
	//set up scene info
  	context->setRayTypeCount( 1u );
  	context->setEntryPointCount( 1u );
  	context["radiance_ray_type"]->setUint( 0u );
  	context["scene_epsilon"]->setFloat( 1.e-4f );
	context->setPrintEnabled( 1);
	printf_size = context->getPrintBufferSize();
	context->setPrintBufferSize(printf_size*10);
	context->setExceptionEnabled( RT_EXCEPTION_ALL, 1);

	// set stack size
	stack_size = context->getStackSize();
	stack_size = stack_size_multiplier*stack_size;
	context->setStackSize( stack_size );
	
	// Render particle buffer and attach to variable, get pointer for CUDA
	positions_buffer = context->createBuffer(RT_BUFFER_INPUT_OUTPUT,RT_FORMAT_USER,N);
	positions_buffer -> setElementSize( sizeof(source_point) );
	positions_buffer -> getDevicePointer(optix_device,&positions_ptr);  // 0 is optix device
	positions_var = context["positions_buffer"];
	positions_var -> set(positions_buffer);

	// Render reaction buffer and attach to variable, get pointer for CUDA
	rxn_buffer = context->createBuffer(RT_BUFFER_INPUT_OUTPUT,RT_FORMAT_USER,N);
	rxn_buffer -> setElementSize( sizeof(unsigned) );
	rxn_buffer -> getDevicePointer(optix_device,&rxn_ptr);
	rxn_var = context["rxn_buffer"];
	rxn_var -> set(rxn_buffer);

	// Render done buffer and attach to variable, get pointer for CUDA
	done_buffer = context->createBuffer(RT_BUFFER_INPUT_OUTPUT,RT_FORMAT_USER,N);
	done_buffer -> setElementSize( sizeof(unsigned) );
	done_buffer -> getDevicePointer(optix_device,&done_ptr);
	done_var = context["done_buffer"];
	done_var -> set(done_buffer);

	// Render cellnum buffer and attach to variable, get pointer for CUDA
	cellnum_buffer = context->createBuffer(RT_BUFFER_INPUT_OUTPUT,RT_FORMAT_USER,N);
	cellnum_buffer -> setElementSize( sizeof(unsigned) );
	cellnum_buffer -> getDevicePointer(optix_device,&cellnum_ptr);
	cellnum_var = context["cellnum_buffer"];
	cellnum_var -> set(cellnum_buffer);

	// Render matnum buffer and attach to variable, get pointer for CUDA
	matnum_buffer = context->createBuffer(RT_BUFFER_INPUT_OUTPUT,RT_FORMAT_USER,N);
	matnum_buffer -> setElementSize( sizeof(unsigned) );
	matnum_buffer -> getDevicePointer(optix_device,&matnum_ptr);
	matnum_var = context["matnum_buffer"];
	matnum_var -> set(matnum_buffer);

	// Render remap buffer and attach to variable, get pointer for CUDA
	remap_buffer = context->createBuffer(RT_BUFFER_INPUT_OUTPUT,RT_FORMAT_USER,N);
	remap_buffer -> setElementSize( sizeof(unsigned) );
	remap_buffer -> getDevicePointer(optix_device,&remap_ptr);
	remap_var = context["remap_buffer"];
	remap_var -> set(remap_buffer);

	// Ray generation program 
	sprintf( path_to_ptx, "%s", "camera.ptx" );
	ray_gen_program = context->createProgramFromPTXFile( path_to_ptx, "camera" );
	context->setRayGenerationProgram( 0, ray_gen_program );
	
	// Exception program
	exception_program = context->createProgramFromPTXFile( path_to_ptx, "exception" );
	context->setExceptionProgram( 0, exception_program ); 
	
	// Miss program 
	sprintf( path_to_ptx, "%s", "miss.ptx" );
	miss_program = context->createProgramFromPTXFile( path_to_ptx, "miss" );
	context->setMissProgram( 0, miss_program ); 

	//set boundary condition for outer cell
	context["boundary_condition"]->setUint(0);

	//set trace type, 1=transport (writes intersection point and next cell), 2=fission (writes origin and current cell)
	context["trace_type"]->setUint(1);

	//set outer cell adn get its dimensions
	context["outer_cell"]->setUint(outer_cell);
	outer_cell_type = problem_geom.get_outer_cell_dims(outer_cell_dims);

	// make all geometry instances
	if(GEOM_FLAG==1){
		make_geom_xform(problem_geom);
	}
	else if(GEOM_FLAG==2){
		make_geom_xform_common(problem_geom);
	} 
	else{
		make_geom_prim(problem_geom);
	}

	//set exceptions on
	context->setExceptionEnabled( RT_EXCEPTION_ALL, 1 );

	//validate and compile
	context->validate();
	context->compile();

}
void optix_stuff::init(wgeometry problem_geom, unsigned compute_device_in, std::string accel_type_in){
	// set min and max cell numbers
	mincell = problem_geom.get_minimum_cell();
	maxcell = problem_geom.get_maximum_cell();
	outer_cell = problem_geom.get_outer_cell();
	outer_cell_type = problem_geom.get_outer_cell_type();
	// get material numbers
	n_materials = problem_geom.get_material_count();
	// try to init optix
	try {
		init_internal(problem_geom, compute_device_in, accel_type_in);	
	} 
	catch( optix::Exception &e ){
		std::cout << "OptiX Error in init:" << e.getErrorString().c_str() <<"\n";
		exit(1);
	}
}
void optix_stuff::set_trace_type(unsigned trace_type){
	context["trace_type"]->setUint(trace_type);
}
void optix_stuff::trace(unsigned trace_type){
	context["trace_type"]->setUint(trace_type);
	context -> launch( 0 , N );
}
void optix_stuff::trace(unsigned trace_type, unsigned n_active){
	context["trace_type"]->setUint(trace_type);
	context -> launch( 0 , n_active );
}
void optix_stuff::trace(){
	context -> launch( 0 , N );
}
void optix_stuff::make_geom_xform(wgeometry problem_geom){

	using namespace optix;

	Group 			top_level_group;
	Variable 		top_object;
	Acceleration 		top_level_acceleration;
	Acceleration 		this_accel;
	GeometryGroup 		this_geom_group;
	Variable 		this_geom_min;
	Variable 		this_geom_max;
	Geometry 		this_geom;
	GeometryInstance 	ginst;
	Material 		material;
	Program  		intersection_program;
	Program  		bounding_box_program;
	Program  		closest_hit_program;
	Transform 		this_transform;
	Acceleration  		acceleration;
	Variable  		cellnum_var;
	Variable  		cellmat_var;
	Variable 		cellfissile_var;
	char 			path_to_ptx[512];
	unsigned 		cellnum,cellmat;
	float 			dx,dy,dz,theta,phi;
	float 			m[16];
	unsigned 		uniqueindex = 0;
	unsigned 		is_fissile = 0;

	// Make top level group/accel as children of the top level object
	this_accel 	= context -> createAcceleration(accel_type.c_str(),traverse_type.c_str());
	this_accel 	-> markDirty();
	top_level_group = context->createGroup();
	top_level_group ->setChildCount(problem_geom.get_transform_count());   // every primitive has at least 1 transform, so the total number of transforms is the number of instances
	top_level_group -> setAcceleration( this_accel );
	context["top_object"]-> set( top_level_group );

	for(int j=0;j<problem_geom.get_primitive_count();j++){

		//create this geometry type
		this_geom = context->createGeometry();
		this_geom -> setPrimitiveCount(1u);

		//set intersection and BB programs
		if      (problem_geom.primitives[j].type == 0)	{sprintf( path_to_ptx, "%s", "box.ptx" );}
		else if (problem_geom.primitives[j].type == 1)	{sprintf( path_to_ptx, "%s", "cylinder.ptx" );}
		else if (problem_geom.primitives[j].type == 2)	{sprintf( path_to_ptx, "%s", "hex.ptx" );}
		bounding_box_program = context->createProgramFromPTXFile( path_to_ptx, "bounds" );
		intersection_program = context->createProgramFromPTXFile( path_to_ptx, "intersect" );
		this_geom -> setBoundingBoxProgram ( bounding_box_program );
		this_geom -> setIntersectionProgram( intersection_program );

		//set hit programs to material
		sprintf( path_to_ptx, "%s", "hits.ptx" );
		closest_hit_program = context->createProgramFromPTXFile( path_to_ptx, "closest_hit" );
		material = context -> createMaterial();
		material -> setClosestHitProgram( 0, closest_hit_program );

		//set program variables for this instance
    		this_geom_min = this_geom["mins"];
    		this_geom_max = this_geom["maxs"];
    		this_geom_min -> set3fv( problem_geom.primitives[j].min );
    		this_geom_max -> set3fv( problem_geom.primitives[j].max );

		for (int k=0;k<problem_geom.primitives[j].n_transforms;k++){

			dx =          problem_geom.primitives[j].transforms[k].dx;
			dy =          problem_geom.primitives[j].transforms[k].dy;
			dz =          problem_geom.primitives[j].transforms[k].dz;
			theta =       problem_geom.primitives[j].transforms[k].theta;
			phi =         problem_geom.primitives[j].transforms[k].phi;
			cellnum =     problem_geom.primitives[j].transforms[k].cellnum;
			cellmat =     problem_geom.primitives[j].transforms[k].cellmat;
			for(int z=0;z<problem_geom.get_material_count();z++){
				if (cellmat == problem_geom.materials[z].matnum){
					is_fissile = problem_geom.materials[z].is_fissile;   // set fissile flag
					cellmat    = problem_geom.materials[z].id;            // hash the material number to the ID, which is the matrix index, not that user-set number
					break;
				}
			}
			//std::cout << "cellnum " << cellnum << " matnum " << cellmat << " is fissile " << is_fissile << "\n";

			//create instances
			ginst = context -> createGeometryInstance();
			ginst -> setGeometry( this_geom );
			ginst -> setMaterialCount( 1u );
			ginst -> setMaterial( 0, material );

			//set cell-specific variables
			cellnum_var 	= ginst["cellnum"];
			cellmat_var 	= ginst["cellmat"];
			cellfissile_var = ginst["cellfissile"];
			cellnum_var 	-> setUint(cellnum);
			cellmat_var 	-> setUint(cellmat);
			cellfissile_var -> setUint(is_fissile);
			//std::cout << "cellnum,matnum,isfiss " << cellnum << " " << cellmat << " " << is_fissile << "\n";

			// make geometry group for this primitive (to attach acceleration to)
			this_accel = context->createAcceleration(accel_type.c_str(),traverse_type.c_str());
			this_accel -> markDirty();
			this_geom_group = context -> createGeometryGroup();
			this_geom_group -> setChildCount( 1u );
			this_geom_group -> setAcceleration( this_accel );
        
			//put geom instance into geomgroup
			this_geom_group -> setChild( 0, ginst );
    
        	//make transforms as necessary and attach to root node
        	//printf("cell %d: applying transform %d -  dx=%f dy=%f dz=%f theta=%f phi=%f\n",cellnum,k,dx,dy,dz,theta,phi);

			m[ 0] = cos(theta)*cos(phi);    m[ 1] = -cos(theta)*sin(phi);   m[ 2] = sin(theta);     m[ 3] = dx;
			m[ 4] = sin(phi);               m[ 5] = cos(phi);               m[ 6] = 0.0f;           m[ 7] = dy;
			m[ 8] = -sin(theta)*cos(phi);   m[ 9] = sin(theta)*sin(phi);    m[10] = cos(theta);     m[11] = dz;
			m[12] = 0.0f;                   m[13] = 0.0f;                   m[14] = 0.0f;           m[15] = 1.0f;
  
			this_transform = context -> createTransform();
			this_transform -> setChild(this_geom_group);
			this_transform -> setMatrix( 0, m, 0 );
			top_level_group -> setChild( uniqueindex , this_transform );
			uniqueindex++;

		}

	}

}
void optix_stuff::make_geom_xform_common(wgeometry problem_geom){

	using namespace optix;

	Group 			top_level_group;
	Variable 		top_object;
	Acceleration 		top_level_acceleration;
	Acceleration 		this_accel;
	GeometryGroup 		this_geom_group;
	Variable 		this_geom_min;
	Variable 		this_geom_max;
	Geometry 		this_geom;
	GeometryInstance 	ginst;
	Material 		material;
	Program  		intersection_program;
	Program  		bounding_box_program;
	Program  		closest_hit_program;
	Transform 		this_transform;
	Acceleration  		acceleration;
	Variable  		cellnum_var;
	Variable  		cellmat_var;
	Variable 		cellfissile_var;
	char 			path_to_ptx[512];
	unsigned 		cellnum,cellmat;
	float 			dx,dy,dz,theta,phi;
	float 			m[16];
	unsigned 		uniqueindex = 0;
	unsigned 		is_fissile = 0;

	// Make top level group/accel as children of the top level object
	this_accel 	= context -> createAcceleration(accel_type.c_str(),traverse_type.c_str());
	this_accel 	-> markDirty();
	top_level_group = context->createGroup();
	top_level_group ->setChildCount(problem_geom.get_transform_count());   // every primitive has at least 1 transform, so the total number of transforms is the number of instances
	top_level_group -> setAcceleration( this_accel );
	context["top_object"]-> set( top_level_group );

	for(int j=0;j<problem_geom.get_primitive_count();j++){

		//create this geometry type
		this_geom = context->createGeometry();
		this_geom -> setPrimitiveCount(1u);

		//set intersection and BB programs
		if      (problem_geom.primitives[j].type == 0)	{sprintf( path_to_ptx, "%s", "box.ptx" );}
		else if (problem_geom.primitives[j].type == 1)	{sprintf( path_to_ptx, "%s", "cylinder.ptx" );}
		else if (problem_geom.primitives[j].type == 2)	{sprintf( path_to_ptx, "%s", "hex.ptx" );}
		bounding_box_program = context->createProgramFromPTXFile( path_to_ptx, "bounds" );
		intersection_program = context->createProgramFromPTXFile( path_to_ptx, "intersect" );
		this_geom -> setBoundingBoxProgram ( bounding_box_program );
		this_geom -> setIntersectionProgram( intersection_program );

		//set hit programs to material
		sprintf( path_to_ptx, "%s", "hits.ptx" );
		closest_hit_program = context->createProgramFromPTXFile( path_to_ptx, "closest_hit" );
		material = context -> createMaterial();
		material -> setClosestHitProgram( 0, closest_hit_program );

		//set program variables for this instance
    		this_geom_min = this_geom["mins"];
    		this_geom_max = this_geom["maxs"];
    		this_geom_min -> set3fv( problem_geom.primitives[j].min );
    		this_geom_max -> set3fv( problem_geom.primitives[j].max );

    		//create instances
		ginst = context -> createGeometryInstance();
		ginst -> setGeometry( this_geom );
		ginst -> setMaterialCount( 1u );
		ginst -> setMaterial( 0, material );

		//set cell-specific variables
		cellnum_var 	= ginst["cellnum"];
		cellmat_var 	= ginst["cellmat"];
		cellfissile_var = ginst["cellfissile"];
		cellnum_var 	-> setUint(cellnum);
		cellmat_var 	-> setUint(cellmat);
		cellfissile_var -> setUint(is_fissile);
		//std::cout << "cellnum,matnum,isfiss " << cellnum << " " << cellmat << " " << is_fissile << "\n";

		// make geometry group for this primitive (to attach acceleration to)
		this_accel = context->createAcceleration(accel_type.c_str(),traverse_type.c_str());
		this_accel -> markDirty();
		this_geom_group = context -> createGeometryGroup();
		this_geom_group -> setChildCount( 1u );
		this_geom_group -> setAcceleration( this_accel );

		//put geom instance into geomgroup
		this_geom_group -> setChild( 0, ginst );

		for (int k=0;k<problem_geom.primitives[j].n_transforms;k++){

			dx =          problem_geom.primitives[j].transforms[k].dx;
			dy =          problem_geom.primitives[j].transforms[k].dy;
			dz =          problem_geom.primitives[j].transforms[k].dz;
			theta =       problem_geom.primitives[j].transforms[k].theta;
			phi =         problem_geom.primitives[j].transforms[k].phi;
			cellnum =     problem_geom.primitives[j].transforms[k].cellnum;
			cellmat =     problem_geom.primitives[j].transforms[k].cellmat;
			for(int z=0;z<problem_geom.get_material_count();z++){
				if (cellmat == problem_geom.materials[z].matnum){
					is_fissile =  problem_geom.materials[z].is_fissile;   // set fissile flag
					cellmat    = problem_geom.materials[z].id;            // hash the material number to the ID, which is the matrix index, not that user-set number
					break;
				}
			}
			//std::cout << "cellnum " << cellnum << " matnum " << cellmat << " is fissile " << is_fissile << "\n";
    
        	//make transforms as necessary and attach to root node
        	//printf("cell %d: applying transform %d -  dx=%f dy=%f dz=%f theta=%f phi=%f\n",cellnum,k,dx,dy,dz,theta,phi);

			m[ 0] = cos(theta)*cos(phi);    m[ 1] = -cos(theta)*sin(phi);   m[ 2] = sin(theta);     m[ 3] = dx;
			m[ 4] = sin(phi);               m[ 5] = cos(phi);               m[ 6] = 0.0f;           m[ 7] = dy;
			m[ 8] = -sin(theta)*cos(phi);   m[ 9] = sin(theta)*sin(phi);    m[10] = cos(theta);     m[11] = dz;
			m[12] = 0.0f;                   m[13] = 0.0f;                   m[14] = 0.0f;           m[15] = 1.0f;
  
			this_transform = context -> createTransform();
			this_transform -> setChild(this_geom_group);
			this_transform -> setMatrix( 0, m, 0 );
			top_level_group -> setChild( uniqueindex , this_transform );
			uniqueindex++;

		}

	}

}
void optix_stuff::make_geom_prim(wgeometry problem_geom){

	using namespace optix;

	Group 			top_level_group;
	Variable 		top_object;
	Acceleration 		this_accel;
	Buffer 			geom_buffer;
	geom_data* 		geom_buffer_host;
	geom_data*		geom_buffer_ptr;
	GeometryGroup 		this_geom_group;
	Geometry 		this_geom;
	GeometryInstance 	ginst;
	Material 		material;
	Program  		intersection_program;
	Program  		bounding_box_program;
	Program  		closest_hit_program;
	Transform 		this_transform;
	Acceleration  		acceleration;
	char 			path_to_ptx[512];
	unsigned 		cellnum,cellmat;
	float 			dx,dy,dz,theta,phi;
	float 			m[16];
	unsigned 		uniqueindex = 0;
	unsigned 		is_fissile = 0;

	// Make top level group/accel as children of the top level object
	this_accel = context->createAcceleration(accel_type.c_str(),traverse_type.c_str());
	this_accel -> markDirty();
	this_geom_group = context -> createGeometryGroup();
	this_geom_group -> setChildCount( problem_geom.get_primitive_count() );
	this_geom_group -> setAcceleration( this_accel );
	context["top_object"]-> set( this_geom_group );

	for(int j=0 ; j<problem_geom.get_primitive_count() ; j++){

		//create this geometry type
		this_geom = context->createGeometry();
		this_geom -> setPrimitiveCount( problem_geom.primitives[j].n_transforms );

		//set intersection and BB programs
		if      (problem_geom.primitives[j].type == 0)	{sprintf( path_to_ptx, "%s", "box_mesh.ptx" );}
		else if (problem_geom.primitives[j].type == 1)	{sprintf( path_to_ptx, "%s", "cylinder_mesh.ptx" );}
		else if (problem_geom.primitives[j].type == 2)	{sprintf( path_to_ptx, "%s", "hex_mesh.ptx" );}
		else if (problem_geom.primitives[j].type == 3)	{sprintf( path_to_ptx, "%s", "sphere_mesh.ptx" );}
		bounding_box_program = context->createProgramFromPTXFile( path_to_ptx, "bounds" );
		intersection_program = context->createProgramFromPTXFile( path_to_ptx, "intersect" );
		this_geom -> setBoundingBoxProgram ( bounding_box_program );
		this_geom -> setIntersectionProgram( intersection_program );

		//set hit programs to material
		sprintf( path_to_ptx, "%s", "hits_mesh.ptx" );
		closest_hit_program = context->createProgramFromPTXFile( path_to_ptx, "closest_hit" );
		material = context -> createMaterial();
		material -> setClosestHitProgram( 0, closest_hit_program );

		//create instance
		ginst = context -> createGeometryInstance();
		ginst -> setGeometry( this_geom );
		ginst -> setMaterialCount( 1u );
		ginst -> setMaterial( 0, material );

		// create internal optix buffer and copy "transform" data, as well as cellnum, matnum, and isfiss into the buffer stuct
		geom_buffer = context->createBuffer(RT_BUFFER_INPUT,RT_FORMAT_USER,problem_geom.primitives[j].n_transforms);
		geom_buffer -> setElementSize( sizeof(geom_data) );
		this_geom["geom_data"] -> set( geom_buffer );
		geom_buffer_ptr = (geom_data*) geom_buffer->map();
		for(int k=0 ; k<problem_geom.primitives[j].n_transforms ; k++){
			// apply transforms to base primitive
			geom_buffer_ptr[k].min[0]		= problem_geom.primitives[j].min[0];
			geom_buffer_ptr[k].min[1]		= problem_geom.primitives[j].min[1];
			geom_buffer_ptr[k].min[2]		= problem_geom.primitives[j].min[2];
			geom_buffer_ptr[k].max[0]		= problem_geom.primitives[j].max[0];
			geom_buffer_ptr[k].max[1]		= problem_geom.primitives[j].max[1];
			geom_buffer_ptr[k].max[2]		= problem_geom.primitives[j].max[2];
			geom_buffer_ptr[k].loc[0]		= problem_geom.primitives[j].transforms[k].dx;
			geom_buffer_ptr[k].loc[1]		= problem_geom.primitives[j].transforms[k].dy;
			geom_buffer_ptr[k].loc[2]		= problem_geom.primitives[j].transforms[k].dz;
			// cell properties
			geom_buffer_ptr[k].cellnum		= problem_geom.primitives[j].transforms[k].cellnum;
			geom_buffer_ptr[k].matnum		= problem_geom.primitives[j].transforms[k].cellmat;
			for(int z=0;z<problem_geom.get_material_count();z++){  // resolve the material number (user input) to material ID (index) to be used in the phyics routines
				if (geom_buffer_ptr[k].matnum == problem_geom.materials[z].matnum){
					geom_buffer_ptr[k].is_fissile =  problem_geom.materials[z].is_fissile;   // set fissile flag
					geom_buffer_ptr[k].matnum    =  problem_geom.materials[z].id;            // hash the material number to the ID, which is the matrix index, not that user-set number
					break;
				}
			}
			//std::cout << "Prim/xform " << j << "," << k << " " << geom_buffer_ptr[k].min[0] << " "	<< geom_buffer_ptr[k].min[1]	<< " "	<< geom_buffer_ptr[k].min[2]	<< " "	<< geom_buffer_ptr[k].max[0]	<< " "	<< geom_buffer_ptr[k].max[1]	<< " "	<< geom_buffer_ptr[k].max[2]	<< " "	<< geom_buffer_ptr[k].cellnum	<< " "	<< geom_buffer_ptr[k].matnum	<< " "	<< geom_buffer_ptr[k].is_fissile << "\n";
		}
		geom_buffer->unmap();

		//set program variables for this object type
    	this_geom["dims"] -> setBuffer( geom_buffer );

    	//set instance as geomgroup child
    	this_geom_group -> setChild(j,ginst);

	}

}
void optix_stuff::trace_geometry(unsigned width_in,unsigned height_in,std::string type,std::string filename){

//	std::cout << "\e[1;32m" << "Plotting Geometry... " << "\e[m \n";
//
//	using namespace optix;
//
//	//get aspect ratio and make N-compatible corresponding heights and widths
//	float aspect = width_in / height_in;
//	float mu, theta;
//	float pi=3.14159;
//	unsigned width  = sqrtf(N*aspect); 
//	unsigned height = sqrtf(N/aspect);
//	std::cout << "width  = " << width << "\n";
//	std::cout << "height = " << height << "\n";
//
//	// init the starting points to be across the z=0 plane and pointing downwards or isotropically random, should produce the same results
//	//FILE* positionsfile = fopen("positionsfile","w");
//	source_point * positions_local = new source_point[width*height];
//	float dx = (42.0-(-42.0))/width;
//	float dy = (42.0-(-42.0))/height;
//	unsigned index;
//	for(int j=0;j<height;j++){
//		for(int k=0;k<width;k++){
//			mu = 2.0*rand()-1.0;
//			theta = 2.0*pi*rand();
//			index = j * width + k;
//			positions_local[index].x = -42.0 + dx/2 + k*dx;
//			positions_local[index].y = -42.0 + dy/2 + j*dy;
//			positions_local[index].z = 0.0;
//			positions_local[index].xhat = sqrtf(1-mu*mu) * cosf( theta ); //0.0;
//			positions_local[index].yhat = sqrtf(1-mu*mu) * sinf( theta ); //0.0;
//			positions_local[index].zhat =       mu; //-1.0;
//			positions_local[index].samp_dist = 50000.0; 
//		}
//	}
//	//fclose(positionsfile);
//
//	// copy starting positions data to pointer
//	cudaMemcpy((void*)positions_ptr,positions_local,width*height*sizeof(source_point),cudaMemcpyHostToDevice);
//	
//	// trace with whereami?
//	context["trace_type"]->setUint(2);
//	context->launch(0,width*height);
//	
//	//copy to local buffer
//	unsigned * image_local = new unsigned[width*height];
//	cudaMemcpy(image_local,(void*)cellnum_ptr,width*height*sizeof(unsigned),cudaMemcpyDeviceToHost);
//
//	// make image
//	png::image< png::rgb_pixel > image(height, width);
//	float * colormap = new float[3];
//	for (size_t y = 0; y < image.get_height(); ++y)
//	{
//	    for (size_t x = 0; x < image.get_width(); ++x)
//	    {
//	    	//mincell=0;
//	    	//maxcell=3;
//	    	make_color(colormap,image_local[y*width+x],mincell,maxcell);
//	    	//printf("%u %u %6.3f %6.3f %6.3f\n",mincell,maxcell,colormap[0],colormap[1],colormap[2]);
//	        image[y][x] = png::rgb_pixel(colormap[0],colormap[1],colormap[2]);
//	    }
//	}
//
//	image.write(filename);
//
//	std::cout << "Done.  Written to " << filename << "\n";
//
//	delete image_local;
//	delete colormap;
//	delete positions_local;
//
}
float optix_stuff::trace_test(){

	float mu, theta;
	float pi=3.14159;

	//FILE* positionsfile = fopen("positionsfile","w");
	source_point * positions_local = new source_point[N];
	unsigned index;
	float x_min = outer_cell_dims[0];
	float y_min = outer_cell_dims[1];
	float z_min = outer_cell_dims[2];
	float x_max = outer_cell_dims[3];
	float y_max = outer_cell_dims[4];
	float z_max = outer_cell_dims[5];

	//y_min = y_max = 0.0;

	// make distribution random now
	int height = (int) sqrtf(N);
	int width  = (int) sqrtf(N);
	printf("image w/h %dx%d\n",width,height);
	float dx = (x_max-x_min)/width;
	float dy = (x_max-x_min)/height;
	float dz = (z_max-z_min)/height;
	for(int j=0;j<height;j++){
		for(int k=0;k<width;k++){
			mu = 2.0*get_rand()-1.0;
			theta = 2.0*pi*get_rand();
			index = j * width + k;
			positions_local[index].x = x_min + dx/2 + k*dx;
			positions_local[index].y = y_min + dy/2 + j*dy;
			positions_local[index].z = 0.0;//z_min + dz/2 + j*dz;
			positions_local[index].xhat = 	sqrtf(1-mu*mu) * cosf( theta ); 
			positions_local[index].yhat = 	sqrtf(1-mu*mu) * sinf( theta ); 
			positions_local[index].zhat = 	mu; 
			positions_local[index].surf_dist = 50000.0; 
			//printf("%6.4E %6.4E %6.4E %6.4E %6.4E %6.4E\n",positions_local[index].x,positions_local[index].y,positions_local[index].z,positions_local[index].xhat,positions_local[index].yhat,positions_local[index].zhat);
		}
	}
	for(index;index<N;index++){  //these are the end bits that don't fit into a square
			positions_local[index].x = 0.0;
			positions_local[index].y = 0.0;
			positions_local[index].z = 0.0;
			positions_local[index].xhat = 0.0;
			positions_local[index].yhat = 0.0;
			positions_local[index].zhat =-1.0;
			positions_local[index].surf_dist = 50000.0; 
	}

	// copy starting positions data to pointer
	copy_to_device((void*)positions_ptr,positions_local,N*sizeof(source_point));
	
	// trace with with plane to generate image
	std::cout << "\e[1;32m" << "Tracing to image to build accel struct and to prove whereami..." << "\e[m \n";
	trace(2);

	//copy to local buffer
	unsigned * image_local = new unsigned[width*height];
	if(image_type.compare("cell")==0){			copy_from_device(image_local,(void*)cellnum_ptr,width*height*sizeof(unsigned));}
	else if (image_type.compare("material")==0){	copy_from_device(image_local,(void*)matnum_ptr,width*height*sizeof(unsigned));}

	// make image
	png::image< png::rgb_pixel > image(height, width);
	float * colormap = new float[3];
	for (size_t y = 0; y < image.get_height(); ++y)
	{
	    for (size_t x = 0; x < image.get_width(); ++x)
	    {
	    	if(image_type.compare("cell")==0){			make_color(colormap,image_local[y*width+x],mincell,maxcell);}
	    	else if (image_type.compare("material")==0){	make_color(colormap,image_local[y*width+x],0  ,n_materials);}
	        image[y][x] = png::rgb_pixel(colormap[0],colormap[1],colormap[2]);
	    }
	}

	image.write("geom_test.png");

	// make distribution random now
	for(index=0;index<N;index++){
			mu 				 = 2.0*get_rand()-1.0;
			theta				 = 2.0*pi*get_rand();
			positions_local[index].surf_dist =  500000;   
			positions_local[index].x         =     0.9 * ( ( x_max - x_min ) * get_rand() + x_min );  
			positions_local[index].y         =     0.9 * ( ( y_max - y_min ) * get_rand() + y_min );  
			positions_local[index].z         =     0.9 * ( ( z_max - z_min ) * get_rand() + z_min ); 
			positions_local[index].xhat      =     sqrtf(1-mu*mu) * cosf( theta );
			positions_local[index].yhat      =     sqrtf(1-mu*mu) * sinf( theta );
			positions_local[index].zhat      =     mu;
			//printf("%6.4E %6.4E %6.4E %6.4E %6.4E %6.4E\n",positions_local[index].x,positions_local[index].y,positions_local[index].z,positions_local[index].xhat,positions_local[index].yhat,positions_local[index].zhat);
	}

	// copy starting positions data to pointer
	copy_to_device((void*)positions_ptr,positions_local,N*sizeof(source_point));

	//trace and time
	std::cout << "\e[1;32m" << "Timing trace " << N << " particles... " ;
	float time_out = ((float)clock())/((float)CLOCKS_PER_SEC);
	trace(2);
	time_out = ((float)clock())/((float)CLOCKS_PER_SEC) - time_out; 
	std::cout << "\e[1;32m" << " done in " << time_out << " seconds \e[m \n";

	// print out intersection points


	return time_out;

}
void optix_stuff::print(){
	std::string instancing;
	if(GEOM_FLAG==1){instancing="transform";}
	else if(GEOM_FLAG==2){instancing="common node transform";}
	else         {instancing="primitive";}
	std::cout << "\e[1;32m" << "--- OptiX SUMMARY ---" << "\e[m \n";
	std::cout << "  Using \e[1;31m"<< instancing <<"\e[m-based instancing\n";
	std::cout << "  Image type set to \e[1;31m"<< image_type <<"\e[m\n";	
	std::cout << "  Device set to "<<compute_device<<"\n";
	std::cout << "  Acceleration set to "<<accel_type<<"/"<<traverse_type<<"\n";
	std::cout << "  stack  size = " << context->getStackSize() << " bytes\n";
	std::cout << "  printf size = " << context->getPrintBufferSize() << " bytes\n";
}
void optix_stuff::make_color(float* color, unsigned x, unsigned min, unsigned max){
	// red linear blue linear green sin colormap
	if (x==4294967295){ //miss value
		color[0]=0.0;
		color[1]=0.0;
		color[2]=0.0;
	}
	else{
		float normed_value = (float) (x-min+1)/(max+2-min);
		color[0] = normed_value;              //red
		color[1] = sin(normed_value*3.14159); //green
		color[2] = 1.0-normed_value;          //blue
	}
	
	//bring up to 256?
	color[0]=color[0]*256;
	color[1]=color[1]*256;
	color[2]=color[2]*256;

}
float optix_stuff::get_rand(){

	return static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

}
void optix_stuff::set_image_type(std::string string_in){

	if(string_in.compare("material") & string_in.compare("cell") ){
		std::cout << "\"" << string_in << "\" is not a valid option for image type, must be \"cell\" or \"material\"\n";
	}
	else{
		//std::cout << "Image type set to \"" << string_in << "\"\n";
		image_type = string_in;
	}

}
unsigned optix_stuff::get_outer_cell(){

	return outer_cell;

}

unsigned optix_stuff::get_outer_cell_type(){

	return outer_cell_type;

}

