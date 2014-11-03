#include <limits.h>
#include "optix_stuff_test.h"
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

/**
 * \brief constructor test
 *
 * \details tests that the number of histories and the stack size multiplier are both equal to the 
 * user-specified values.
 */
TEST_F(optixStuffTest, Constructor)
{
	N = 10000;
	mult = 4;
	optix_stuff constructor_test(N,mult);
	mult = constructor_test.stack_size_multiplier;
	N = constructor_test.N;
	EXPECT_EQ(4u, mult);
	EXPECT_EQ(10000u, N);
}

/**
 * \brief initialization test
 *
 * \details creates a wgeometry and adds a material, a primitive, and a transform to it. sets the outer 
 * cell of the wgeometry and updates and checks it. sets the number of histories, stack size multiplier, 
 * and image type of the OptiX stuff object. tests that no exceptions are thrown by the initialization.
 */
TEST_F(optixStuffTest, Init)
{
	optix_stuff init_test;
	wgeometry geom;
	compute_device = 0;
	accel_type = "Sbvh";
	material = 1;
	is_fissile = 1;
	n_isotopes = 1;
	density = 19.816;
	isotopes.push_back(94239);
	mat_fracs.push_back(1);
	geom.add_material(material,is_fissile,n_isotopes,density,isotopes,mat_fracs);
	type = 3;
	min[0] = -5.1; min[1] = -5.1; min[2] = -5.1;
	max[0] =  5.1; max[1] =  5.1; max[2] =  5.1;
	origin[0] = 0; origin[1] = 0; origin[2] = 0;
	geom.add_primitive(type, material, min, max, origin);
	geom.add_transform(0,999,0,0,0,0,0);
	geom.set_outer_cell(999);
	geom.update();
	geom.check();
	init_test.N = 10000;
	init_test.stack_size_multiplier = 1;
	init_test.set_image_type("cell");
	EXPECT_NO_THROW(init_test.init(geom,compute_device,accel_type));
}

/**
 * \brief set trace type test
 *
 * \details creates an OptiX stuff object and a wgeometry object. sets wgeometry parameters, initializes 
 * the OptiX object, and sets the trace type.
 */
TEST_F(optixStuffTest, SetTraceType)
{
	optix_stuff set_trace_type_test;
        wgeometry geom;
        compute_device = 0;
        accel_type = "Sbvh";
        material = 1;
        is_fissile = 1;
        n_isotopes = 1;
        density = 19.816;
        isotopes.push_back(94239);
        mat_fracs.push_back(1);
        geom.add_material(material,is_fissile,n_isotopes,density,isotopes,mat_fracs);
        type = 3;
        min[0] = -5.1; min[1] = -5.1; min[2] = -5.1;
        max[0] =  5.1; max[1] =  5.1; max[2] =  5.1;
        origin[0] = 0; origin[1] = 0; origin[2] = 0;
        geom.add_primitive(type, material, min, max, origin);
        geom.add_transform(0,999,0,0,0,0,0);
        geom.set_outer_cell(999);
        geom.update();
        geom.check();
        set_trace_type_test.N = 10000;
        set_trace_type_test.stack_size_multiplier = 1;
        set_trace_type_test.init(geom,compute_device,accel_type);
	set_trace_type_test.set_trace_type(1);
}
/**
 * \brief default trace test
 *
 * \details creates an optix stuff object and a wgeometry object. sets up the wgeometry object and gets
 * its extreme coordinate values. gets random values (angles, surface distances, positions) for every 
 * history. copies this data to the GPU. tests that no exception is thrown from a call to the trace 
 * function.
 */
TEST_F(optixStuffTest, DefaultTrace)
{
	optix_stuff default_trace_test;
	default_trace_test.N=10000;
	default_trace_test.stack_size_multiplier=1;
	wgeometry geom;
        compute_device = 0;
        accel_type = "Sbvh";
        material = 1;
        is_fissile = 1;
        n_isotopes = 1;
        density = 19.816;
        isotopes.push_back(94239);
        mat_fracs.push_back(1);
        geom.add_material(material,is_fissile,n_isotopes,density,isotopes,mat_fracs);
        type = 0; //non-zero type gets RT exceptions
        min[0] = -5.1; min[1] = -5.1; min[2] = -5.1;
        max[0] =  5.1; max[1] =  5.1; max[2] =  5.1;
        origin[0] = 0; origin[1] = 0; origin[2] = 0;
        geom.add_primitive(type, material, min, max, origin);
        geom.add_transform(0,999,0,0,0,0,0);
        geom.set_outer_cell(999);
        geom.update();
        geom.check();
        default_trace_test.init(geom,compute_device,accel_type);
        source_point * positions_local = new source_point[N];
        x_min = default_trace_test.outer_cell_dims[0];
        y_min = default_trace_test.outer_cell_dims[1];
        z_min = default_trace_test.outer_cell_dims[2];
        x_max = default_trace_test.outer_cell_dims[3];
        y_max = default_trace_test.outer_cell_dims[4];
        z_max = default_trace_test.outer_cell_dims[5];
	for(index=0;index<N;index++){
		mu 				 = 2.0*default_trace_test.get_rand()-1.0;
		theta				 = 2.0*pi*default_trace_test.get_rand();
		positions_local[index].surf_dist = 500000;   
		positions_local[index].x         = 0.9*((x_max - x_min) * default_trace_test.get_rand() + x_min); 
		positions_local[index].y         = 0.9*((y_max - y_min) * default_trace_test.get_rand() + y_min); 
		positions_local[index].z         = 0.9*((z_max - z_min) * default_trace_test.get_rand() + z_min);
		positions_local[index].xhat      = sqrtf(1-mu*mu) * cosf(theta);
		positions_local[index].yhat      = sqrtf(1-mu*mu) * sinf(theta);
		positions_local[index].zhat      = mu;
	}
        copy_to_device((void*)default_trace_test.positions_ptr,positions_local,N*sizeof(source_point));
	EXPECT_NO_THROW(default_trace_test.trace());
}

/**
 * \brief typed trace test
 *
 * \details creates an optix stuff object and a wgeometry object. sets up the wgeometry object and gets
 * its extreme coordinate values. gets random values (angles, surface distances, positions) for every 
 * history. copies this data to the GPU. tests that no exception is thrown from a call to the trace 
 * function where the trace type is passed into the function.
 */
TEST_F(optixStuffTest, TypeTrace)
{
	optix_stuff type_trace_test(10000,1);
	wgeometry geom;
        compute_device = 0;
        accel_type = "Sbvh";
        material = 1;
        is_fissile = 1;
        n_isotopes = 1;
        density = 19.816;
        isotopes.push_back(94239);
        mat_fracs.push_back(1);
        geom.add_material(material,is_fissile,n_isotopes,density,isotopes,mat_fracs);
        type = 0; //non-zero type gets RT exceptions
        min[0] = -5.1; min[1] = -5.1; min[2] = -5.1;
        max[0] =  5.1; max[1] =  5.1; max[2] =  5.1;
        origin[0] = 0; origin[1] = 0; origin[2] = 0;
        geom.add_primitive(type, material, min, max, origin);
        geom.add_transform(0,999,0,0,0,0,0);
        geom.set_outer_cell(999);
        geom.update();
        geom.check();
        type_trace_test.init(geom,compute_device,accel_type);
        source_point * positions_local = new source_point[N];
        x_min = type_trace_test.outer_cell_dims[0];
        y_min = type_trace_test.outer_cell_dims[1];
        z_min = type_trace_test.outer_cell_dims[2];
        x_max = type_trace_test.outer_cell_dims[3];
        y_max = type_trace_test.outer_cell_dims[4];
        z_max = type_trace_test.outer_cell_dims[5];
	for(index=0;index<N;index++){
		mu 				 = 2.0*type_trace_test.get_rand()-1.0;
		theta				 = 2.0*pi*type_trace_test.get_rand();
		positions_local[index].surf_dist = 500000;   
		positions_local[index].x         = 0.9*((x_max - x_min) * type_trace_test.get_rand() + x_min); 
		positions_local[index].y         = 0.9*((y_max - y_min) * type_trace_test.get_rand() + y_min); 
		positions_local[index].z         = 0.9*((z_max - z_min) * type_trace_test.get_rand() + z_min);
		positions_local[index].xhat      = sqrtf(1-mu*mu) * cosf(theta);
		positions_local[index].yhat      = sqrtf(1-mu*mu) * sinf(theta);
		positions_local[index].zhat      = mu;
	}
        copy_to_device((void*)type_trace_test.positions_ptr,positions_local,N*sizeof(source_point));
	EXPECT_NO_THROW(type_trace_test.trace(3)); //type 2 causes RT_EXCEPTION_BUFFER_INDEX_OUT_OF_BOUNDS
}

/**
 * \brief history-numbered typed trace test
 *
 * \details creates an optix stuff object and a wgeometry object. sets up the wgeometry object and gets
 * its extreme coordinate values. gets random values (angles, surface distances, positions) for every 
 * history. copies this data to the GPU. tests that no exception is thrown from a call to the trace 
 * function where the number of histories and trace type are passed into the function.
 */
TEST_F(optixStuffTest, NumTypeTrace)
{
	optix_stuff num_type_trace_test(10000,1);
	wgeometry geom;
        compute_device = 0;
        accel_type = "Sbvh";
        material = 1;
        is_fissile = 1;
        n_isotopes = 1;
        density = 19.816;
        isotopes.push_back(94239);
        mat_fracs.push_back(1);
        geom.add_material(material,is_fissile,n_isotopes,density,isotopes,mat_fracs);
        type = 0; //non-zero type gets RT exceptions
        min[0] = -5.1; min[1] = -5.1; min[2] = -5.1;
        max[0] =  5.1; max[1] =  5.1; max[2] =  5.1;
        origin[0] = 0; origin[1] = 0; origin[2] = 0;
        geom.add_primitive(type, material, min, max, origin);
        geom.add_transform(0,999,0,0,0,0,0);
        geom.set_outer_cell(999);
        geom.update();
        geom.check();
        num_type_trace_test.init(geom,compute_device,accel_type);
        source_point * positions_local = new source_point[N];
        x_min = num_type_trace_test.outer_cell_dims[0];
        y_min = num_type_trace_test.outer_cell_dims[1];
        z_min = num_type_trace_test.outer_cell_dims[2];
        x_max = num_type_trace_test.outer_cell_dims[3];
        y_max = num_type_trace_test.outer_cell_dims[4];
        z_max = num_type_trace_test.outer_cell_dims[5];
	for(index=0;index<N;index++){
		mu 				 = 2.0*num_type_trace_test.get_rand()-1.0;
		theta				 = 2.0*pi*num_type_trace_test.get_rand();
		positions_local[index].surf_dist = 500000;   
		positions_local[index].x         = 0.9*((x_max - x_min) * num_type_trace_test.get_rand() + x_min); 
		positions_local[index].y         = 0.9*((y_max - y_min) * num_type_trace_test.get_rand() + y_min); 
		positions_local[index].z         = 0.9*((z_max - z_min) * num_type_trace_test.get_rand() + z_min);
		positions_local[index].xhat      = sqrtf(1-mu*mu) * cosf(theta);
		positions_local[index].yhat      = sqrtf(1-mu*mu) * sinf(theta);
		positions_local[index].zhat      = mu;
	}
        copy_to_device((void*)num_type_trace_test.positions_ptr,positions_local,N*sizeof(source_point));
	EXPECT_NO_THROW(num_type_trace_test.trace(3,5000)); //type 2 causes RT_EXCEPTION_BUFFER_INDEX_OUT_OF_BOUNDS
}

/**
 * \brief get random number test
 *
 * \details gets a random number and tests that it is between zero and one (inclusive ).
 */
TEST_F(optixStuffTest, GetRand)
{
	optix_stuff get_rand_test;
	rand = get_rand_test.get_rand();
	EXPECT_LE(0, rand);
	EXPECT_GE(1, rand);
}

/**
 * \brief set image type test
 *
 * \details tests that an exception is not thrown when either of the permissible image types are set.
 */
TEST_F(optixStuffTest, SetImageType)
{
	optix_stuff image_type_test;
	EXPECT_NO_THROW(image_type_test.set_image_type("cell"));
	EXPECT_NO_THROW(image_type_test.set_image_type("material"));
}

/**
 * \brief get outer cell test
 *
 * \details creates a wgeometry object and an optix stuff object. initializes the optix object with the 
 * geometry and tests that the 'get outer cell' function returns the set value.
 */
TEST_F(optixStuffTest, GetOuterCell)
{
        optix_stuff outer_cell_test;
        wgeometry geom;
        compute_device = 0;
        accel_type = "Sbvh";
        material = 1;
        is_fissile = 1;
        n_isotopes = 1;
        density = 19.816;
        isotopes.push_back(94239);
        mat_fracs.push_back(1);
        geom.add_material(material,is_fissile,n_isotopes,density,isotopes,mat_fracs);
        type = 3;
        min[0] = -5.1; min[1] = -5.1; min[2] = -5.1;
        max[0] =  5.1; max[1] =  5.1; max[2] =  5.1;
        origin[0] = 0; origin[1] = 0; origin[2] = 0;
        geom.add_primitive(type, material, min, max, origin);
        geom.add_transform(0,999,0,0,0,0,0);
        geom.set_outer_cell(999);
        geom.update();
        geom.check();
        outer_cell_test.N = 10000;
        outer_cell_test.stack_size_multiplier = 1;
	outer_cell_test.init(geom,0,"Sbvh");
	outer_cell = outer_cell_test.get_outer_cell();
	EXPECT_EQ(999u, outer_cell);
}
