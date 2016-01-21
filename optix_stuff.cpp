#include <vector> 
#include <iostream>
#include <sstream>
#include <cmath>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <png++/png.hpp>
#include "datadef.h"
#include "wprimitive.h"
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
	Buffer 	 cellnum_buffer;
	Buffer 	 matnum_buffer;
	Buffer 	 talnum_buffer;
	Buffer 	 remap_buffer;
	Variable positions_var;
	Variable rxn_var;
	Variable cellnum_var;
	Variable matnum_var;
	Variable talnum_var;
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
	optix_device = 0;  // always the first enabled device, 0!
	//for(unsigned h=0;h<enabled_count;h++){printf("%u ",enabled_ids[h]); optix_device=enabled_ids[h];}
	//printf("\n");
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
	positions_buffer -> setElementSize( sizeof(spatial_data) );
	positions_buffer -> getDevicePointer(optix_device,&positions_ptr);  // 0 is optix device
	positions_var = context["positions_buffer"];
	positions_var -> set(positions_buffer);

	// Render reaction buffer and attach to variable, get pointer for CUDA
	rxn_buffer = context->createBuffer(RT_BUFFER_INPUT_OUTPUT,RT_FORMAT_USER,N);
	rxn_buffer -> setElementSize( sizeof(unsigned) );
	rxn_buffer -> getDevicePointer(optix_device,&rxn_ptr);
	rxn_var = context["rxn_buffer"];
	rxn_var -> set(rxn_buffer);

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

	// Render talnum buffer and attach to variable, get pointer for CUDA
	talnum_buffer = context->createBuffer(RT_BUFFER_INPUT_OUTPUT,RT_FORMAT_USER,N);
	talnum_buffer -> setElementSize( sizeof(int) );
	talnum_buffer -> getDevicePointer(optix_device,&talnum_ptr);
	talnum_var = context["talnum_buffer"];
	talnum_var -> set(talnum_buffer);

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
	context["boundary_condition"]->setUint(boundary_condition);

	//set trace type, 2=transport (finds nearest surface, normal, writes cell number and material number), 3=fissile query(writes fissile flag into material number, writes cell number), 4=geometry plot(same as 2, but misses are squelched, no normals/intersection distances reported)
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
	boundary_condition = problem_geom.get_boundary_condition();
	if(boundary_condition == 0) {
	 	printf("BC of outer cell not set!  Setting to 1 (black)...\n");
	 	boundary_condition = 1;
	 }
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

	printf("TRANSFORM-BASED INSTANTIATION ROUTINES ARE DEPRECATED AND WILL NOT WORK RIGHT.  ONLY KEPT FOR POSTERITY.\n");
	exit(0);

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
			geom_buffer_ptr[k].talnum		= problem_geom.primitives[j].transforms[k].tally_index;
			for(int z=0;z<problem_geom.get_material_count();z++){  // resolve the material number (user input) to material ID (index) to be used in the phyics routines
				if (geom_buffer_ptr[k].matnum == problem_geom.materials[z].matnum){
					geom_buffer_ptr[k].is_fissile =  problem_geom.materials[z].is_fissile;    // set fissile flag
					geom_buffer_ptr[k].matnum     =  problem_geom.materials[z].id;            // hash the material number to the ID, which is the matrix index, not that user-set number
					break;
				}
			}
			//std::cout << "Prim/xform " << j << "," << k << " " << geom_buffer_ptr[k].min[0] << " "	<< geom_buffer_ptr[k].min[1]	<< " "	<< geom_buffer_ptr[k].min[2]	<< " "	<< geom_buffer_ptr[k].max[0]	<< " "	<< geom_buffer_ptr[k].max[1]	<< " "	<< geom_buffer_ptr[k].max[2]	<< " cellnm "	<< geom_buffer_ptr[k].cellnum	<< " matnum "	<< geom_buffer_ptr[k].matnum	<< " talnum " << geom_buffer_ptr[k].talnum << " fissile "	<< geom_buffer_ptr[k].is_fissile << "\n";
		}
		geom_buffer->unmap();

		//set program variables for this object type
    	this_geom["dims"] -> setBuffer( geom_buffer );

    	//set instance as geomgroup child
    	this_geom_group -> setChild(j,ginst);

	}

}
void optix_stuff::print(){
	std::string instancing;
	if(GEOM_FLAG==1){instancing="transform";}
	else if(GEOM_FLAG==2){instancing="common node transform";}
	else         {instancing="primitive";}
	std::vector<int> enabled_ids = context -> getEnabledDevices();
	std::string enabled_name = context ->getDeviceName(enabled_ids[0]);

	std::cout << "\e[1;32m" << "--- OptiX SUMMARY ---" << "\e[m \n";
	printf(      "  Using device %u: \e[1;31m%s\e[m\n",enabled_ids[0],enabled_name.c_str());	
	std::cout << "  Using \e[1;31m"<< instancing <<"\e[m-based instancing\n";
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
unsigned optix_stuff::get_outer_cell(){

	return outer_cell;

}

unsigned optix_stuff::get_outer_cell_type(){

	return outer_cell_type;

}

