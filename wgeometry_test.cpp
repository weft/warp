#include <limits.h>
#include "wgeometry_test.h"

/**
 * \brief wgeometry construction test
 *
 * \details creates a wgeometry object. tests that the number of boxes, cylinders, hexagons, 
 * primitives, transforms, materials, and isotopes all default to zero. tests that the outer cell and 
 * fissile flag both default to zero.
 */
TEST_F(wgeometryTest, Construction)
{
	wgeometry default_geom;
	n_box = default_geom.n_box;
	n_cyl = default_geom.n_cyl;
        n_hex = default_geom.n_hex;
        n_primitives = default_geom.n_primitives;
        n_transforms = default_geom.n_transforms;
        outer_cell = default_geom.outer_cell;
        n_materials = default_geom.n_materials;
        n_isotopes = default_geom.n_isotopes;
        fissile_flag = default_geom.fissile_flag;
	EXPECT_EQ(0u,n_box);
	EXPECT_EQ(0u,n_cyl);
	EXPECT_EQ(0u,n_hex);
	EXPECT_EQ(0u,n_primitives);
	EXPECT_EQ(0u,n_transforms);
	EXPECT_EQ(0u,outer_cell);
	EXPECT_EQ(0u,n_materials);
	EXPECT_EQ(0u,n_isotopes);
	EXPECT_EQ(0u,fissile_flag);
}

/**
 * \brief add default primtive test
 *
 * \details creates a wgeometry object and adds a default primitive to it. tests that the number of 
 * primtives is incremented and that the primitive is added to the primitive vector.
 */
TEST_F(wgeometryTest, AddDefaultPrimitive)
{
	wgeometry default_prim_test_geom;
	n_primitives = default_prim_test_geom.n_primitives;
	prim_vec_len = default_prim_test_geom.primitives.size();
	EXPECT_EQ(0u,n_primitives);
	EXPECT_EQ(0,prim_vec_len);
	default_prim_test_geom.add_primitive();
	n_primitives = default_prim_test_geom.n_primitives;
	prim_vec_len = default_prim_test_geom.primitives.size();
	EXPECT_EQ(1u,n_primitives);
	EXPECT_EQ(1, prim_vec_len);
}

/**
 * \brief add valued primitive test
 *
 * \details creates a wgeometry object and adds a valued primitive to it. tests that the number of 
 * primitives is incremented, that the primitive is added to the primitive vector, and that the type, 
 * material, and coordinate values (extrema and location) of the primitive are equal to the 
 * user-specified values.
 */
TEST_F(wgeometryTest, AddValuedPrimitive)
{
	wgeometry valued_prim_test_geom;
	min[0] = -1; min[1] = -1; min[2] = -1;
	max[0] =  1; max[1] =  1; max[2] =  1;
	loc[0] =  0; loc[0] =  0; loc[2] =  0; 
	n_primitives = valued_prim_test_geom.n_primitives;
        prim_vec_len = valued_prim_test_geom.primitives.size();
        EXPECT_EQ(0u,n_primitives);
        EXPECT_EQ(0,prim_vec_len);
	valued_prim_test_geom.add_primitive(0,0,min,max,loc);
	n_primitives = valued_prim_test_geom.n_primitives;
        prim_vec_len = valued_prim_test_geom.primitives.size();
	prim_type = valued_prim_test_geom.primitives[0].type;
	material = valued_prim_test_geom.primitives[0].material;
	for(int i = 0; i < 3; i++)
	{
		min[i] = valued_prim_test_geom.primitives[0].min[i];
		max[i] = valued_prim_test_geom.primitives[0].max[i];
		loc[i] = valued_prim_test_geom.primitives[0].location[i];
	}
        EXPECT_EQ(1u,n_primitives);
        EXPECT_EQ(1,prim_vec_len);
	EXPECT_EQ(0u,prim_type);
	EXPECT_EQ(0u,material);
	for(int i = 0; i < 3; i++)
	{
		EXPECT_FLOAT_EQ(-1, min[i]);
		EXPECT_FLOAT_EQ(1, max[i]);
		EXPECT_FLOAT_EQ(0, loc[i]);
	}
}

/**
 * \brief update test
 *
 * \details creates a wgeometry object, adds several primitives to it, and adds transforms to those 
 * primitives. calls the wgeometry update function. tests that the numbers of boxes, cylinders, 
 * hexagons, spheres, transforms, and isotopes are all updated appropriately.
 */
TEST_F(wgeometryTest, Update)
{
	wgeometry update_test_geom; 
	min[0] = -1; min[1] = -1; min[2] = -1;
	max[0] =  1; max[1] =  1; max[2] =  1;
	loc[0] =  0; loc[0] =  0; loc[2] =  0; 
	update_test_geom.add_primitive(0,0,min,max,loc);
	update_test_geom.add_primitive(1,0,min,max,loc);
	update_test_geom.add_primitive(2,0,min,max,loc);
	update_test_geom.add_primitive(3,0,min,max,loc);
	update_test_geom.add_transform(0);
	update_test_geom.add_transform(1);
	update_test_geom.add_transform(2);
	update_test_geom.add_transform(3);
	update_test_geom.update();
	n_box = update_test_geom.n_box;
	n_cyl = update_test_geom.n_cyl;
	n_hex = update_test_geom.n_hex;
	n_sph = update_test_geom.n_sph;
	n_transforms = update_test_geom.n_transforms;
	n_isotopes = update_test_geom.n_isotopes;
	EXPECT_EQ(1u, n_box);
	EXPECT_EQ(1u, n_cyl);
	EXPECT_EQ(1u, n_hex);
	EXPECT_EQ(1u, n_sph);
	EXPECT_EQ(4u, n_transforms);
	EXPECT_EQ(0u, n_isotopes);
}

/**
 * \brief get primitve count test
 *
 * \details creates a wgoemetry object, adds five primitives to it, and calls the 'get primitive count' 
 * function. tests that the function returns the correct number of primitives.
 */
TEST_F(wgeometryTest, GetPrimitiveCount)
{
	wgeometry prim_count_test_geom;
	for(int i = 0; i < 5; i++)
	{
		prim_count_test_geom.add_primitive();
	}
	n_primitives = prim_count_test_geom.get_primitive_count();
	EXPECT_EQ(5u, n_primitives);
}

/**
 * \brief get transform count test
 *
 * \details creates a wgeometry object, adds a default primitive to it, adds five transforms to the 
 * primitive, and calls the update and 'get transform count' functions. tests that the wgeometry is 
 * updated and that the get function returns the correct number of transforms.
 */
TEST_F(wgeometryTest, GetTransformCount)
{
	wgeometry tform_count_test_geom;
	tform_count_test_geom.add_primitive();
	for(int i = 0; i < 5; i++)
	{
		tform_count_test_geom.primitives[0].add_transform();
	}
	tform_count_test_geom.update();
	n_transforms = tform_count_test_geom.get_transform_count();
	EXPECT_EQ(5u, n_transforms);
}

/**
 * \brief get outer cell test
 *
 * \details creates a wgeometry object and sets the outer cell. calls the 'get outer cell' function and 
 * tests that it returns the set outer cell value. adds a primitive and a transform to the wgeometry, 
 * setting a new outer cell value. tests that the 'get outer cell' function returns this new set value.
 */
TEST_F(wgeometryTest, OuterCellTest)
{
	wgeometry outer_cell_test_geom;
	outer_cell_test_geom.set_outer_cell(0,1);
	outer_cell = outer_cell_test_geom.get_outer_cell();
	EXPECT_EQ(0u, outer_cell);
	outer_cell_test_geom.add_primitive();
	outer_cell_test_geom.add_transform(0,1,0,0,0,3.14159,3.14159);
	outer_cell_test_geom.set_outer_cell(1,1);
	outer_cell = outer_cell_test_geom.get_outer_cell();
	EXPECT_EQ(1u, outer_cell);
}

/**
 * \brief get minimum cell test
 *
 * \details creates a wgeometry object and tests that the default minimum cell is -1u. adds a primitive 
 * to the wgeometry object and adds three transforms to that primitive, each with a different minimum 
 * cell. calls the 'get minimum cell' function again and tests that it returns the lowest minimum cell 
 * number added to the wgeometry object.
 */
TEST_F(wgeometryTest, MinCellTest)
{
	wgeometry min_cell_test_geom;
	mincell = min_cell_test_geom.get_minimum_cell();
	EXPECT_EQ(-1u, mincell);
	min_cell_test_geom.add_primitive();
	min_cell_test_geom.add_transform(0,0,0,0,0,0,0);
	min_cell_test_geom.add_transform(0,1,0,0,0,0,0);
	min_cell_test_geom.add_transform(0,999,0,0,0,0,0);
	min_cell_test_geom.update();
	mincell = min_cell_test_geom.get_minimum_cell();
	EXPECT_EQ(0u, mincell);
}

/**
 * \brief get maximum cell test
 *
 * \details creates a wgeometry object and tests that the default maximum cell is 0. adds a primitive 
 * to the wgeometry object and adds three transforms to that primitive, each with a different maximum 
 * cell. calls the 'get maximum cell' function again and tests that it returns the highest maximum cell 
 * number added to the wgeometry object.
 */
TEST_F(wgeometryTest, MaxCellTest)
{
	wgeometry max_cell_test_geom;
	maxcell = max_cell_test_geom.get_maximum_cell();
	EXPECT_EQ(0u, maxcell);
	max_cell_test_geom.add_primitive();
	max_cell_test_geom.add_transform(0,0,0,0,0,0,0);
	max_cell_test_geom.add_transform(0,1,0,0,0,0,0);
	max_cell_test_geom.add_transform(0,999,0,0,0,0,0);
	maxcell = max_cell_test_geom.get_maximum_cell();
	EXPECT_EQ(999u, maxcell);
}

/**
 * \brief add material test
 *
 * \details creates a wgeometry object and adds a material to it. tests that the number of isotopes, 
 * material number, material index, density, and fissile flag are equal to the user-specified values. 
 * tests that the number of materials is incremented and that the material is added to the material 
 * vector.
 */
TEST_F(wgeometryTest, AddMaterial)
{
	wgeometry add_mat_test_geom;
	mat_vec_len = add_mat_test_geom.materials.size();
	n_materials = add_mat_test_geom.n_materials;
	EXPECT_EQ(0, mat_vec_len);
	EXPECT_EQ(0u, n_materials);
	n_isotopes = 2;
	isotopes.push_back(92235);
	isotopes.push_back(92238);
	mat_fracs.push_back(0.1);
	mat_fracs.push_back(0.9);
	density = 15;
	add_mat_test_geom.add_material(0,1,n_isotopes,density,isotopes,mat_fracs);
	n_isotopes = add_mat_test_geom.materials[0].num_isotopes;
	mat_num = add_mat_test_geom.materials[0].matnum;
	mat_index = add_mat_test_geom.materials[0].id;
	density = add_mat_test_geom.materials[0].density;
	is_fissile = add_mat_test_geom.materials[0].is_fissile;
	mat_vec_len = add_mat_test_geom.materials.size();
	n_materials = add_mat_test_geom.n_materials;
	EXPECT_EQ(2u, n_isotopes);
	EXPECT_EQ(0u, mat_num);
	EXPECT_EQ(0, mat_index);
	EXPECT_FLOAT_EQ(15, density);
	EXPECT_EQ(1u, is_fissile);
	EXPECT_EQ(1, mat_vec_len);
	EXPECT_EQ(1u, n_materials);
}

/**
 * \brief check test
 *
 * \details creates a wgeometry object and adds a material, a primitive, and a transform to it. calls 
 * the update function and the check function. tests that the check function returns 0, indicating no 
 * errors.
 */
TEST_F(wgeometryTest, Check)
{
	wgeometry check_test_geom;
	n_isotopes = 2;
	isotopes.push_back(92235);
	isotopes.push_back(92238);
	mat_fracs.push_back(0.1);
	mat_fracs.push_back(0.9);
	density = 15;
	check_test_geom.add_material(0,1,n_isotopes,density,isotopes,mat_fracs);
	prim_type = 0;
	material = 0;
	min[0] = -1; min[1] = -1; min[2] = -1;
	max[0] =  1; max[1] =  1; max[2] =  1;
	loc[0] =  0; loc[0] =  0; loc[2] =  0; 
	prim_id = check_test_geom.add_primitive(prim_type,material,min,max,loc);
	check_test_geom.add_transform(prim_id,999,0,0,0,0,0);
	check_test_geom.set_outer_cell(999,1);
	check_test_geom.update();
	check = check_test_geom.check();
	EXPECT_EQ(0, check);
}

/**
 * \brief get outer cell dimensions test
 *
 * \details creates a wgeometry object, adds a primitive and three transforms to it, and sets the outer 
 * cell. calls the 'get outer cell dims' function and tests that it returns the type of the primitive.
 */
TEST_F(wgeometryTest, GetOuterCellDims)
{
	float dims [6];
	dims[0] = -2; dims[1] = -2; dims[2] = -2;
	dims[3] =  2; dims[4] =  2; dims[5] =  2;
	wgeometry outer_cell_dim_test_geom;
	outer_cell_dim_test_geom.add_primitive();
	outer_cell_dim_test_geom.add_transform(0,0,0,0,0,3.14159,3.14159);
	outer_cell_dim_test_geom.add_transform(0,1,-1,-1,-1,3.14159,3.14159);
	outer_cell_dim_test_geom.add_transform(0,999,1,1,1,3.14159,3.14159);
	outer_cell_dim_test_geom.set_outer_cell(999,1);
	prim_type = outer_cell_dim_test_geom.get_outer_cell_dims(dims);
	EXPECT_EQ(0u,prim_type);
}

/**
 * \brief get material count test
 *
 * \details creates a wgeometry object and adds five materials to it. calls the 'get material count; 
 * function and tests that it returns the correct number of materials.
 */
TEST_F(wgeometryTest, GetMaterialCount)
{
	wgeometry mat_count_test_geom;
	n_isotopes = 1;
	isotopes.push_back(92235);
	mat_fracs.push_back(1.0);
	density = 19;
	is_fissile = 1;
	for(int i = 0; i < 5; i++)
	{
		mat_count_test_geom.add_material(i,is_fissile,n_isotopes,density,isotopes,mat_fracs);
	}
	n_materials = mat_count_test_geom.get_material_count();
	EXPECT_EQ(5u, n_materials);
}

/**
 * \brief check fissile test
 *
 * \details creates a wgeometry object and tests that it defaults to non-fissile. adds a fissile material 
 * to the wgeometry and calls the update, check, and 'is fissile' functions. tests that the wgeometry is 
 * updated to contain a fissile material.
 */
TEST_F(wgeometryTest, CheckFissile)
{
	wgeometry fissile_test_geom;
	is_fissile = fissile_test_geom.fissile_flag;
	EXPECT_EQ(0u, is_fissile);
	n_isotopes = 1;
	isotopes.push_back(92235);
	mat_fracs.push_back(1.0);
	density = 19;
	is_fissile = 1;
	fissile_test_geom.add_material(0,is_fissile,n_isotopes,density,isotopes,mat_fracs);
	fissile_test_geom.add_primitive();
	fissile_test_geom.add_transform(0,0,0,0,0,0,0);
	fissile_test_geom.update();
	fissile_test_geom.check();
	is_fissile = fissile_test_geom.fissile_flag;
	EXPECT_EQ(1u, is_fissile);
}

/**
 * \brief add default transform test
 *
 * \details creates a wgeometry object and adds a default primitive and a default transform to it. tests 
 * that the cell number, material number, transform coordinates, and transform angles default to zero. 
 * tests that the number of transforms is incremented and that the transform is added to the transform 
 * vector.
 */
TEST_F(wgeometryTest, AddDefaultTransform)
{
	wgeometry default_tform_test_geom;
	default_tform_test_geom.add_primitive();
	n_transforms = default_tform_test_geom.n_transforms;
	transform_vec_len = default_tform_test_geom.primitives[0].transforms.size();
	EXPECT_EQ(0u, n_transforms);
	EXPECT_EQ(0, transform_vec_len);
	default_tform_test_geom.add_transform(0);
	default_tform_test_geom.update();
	cellnum = default_tform_test_geom.primitives[0].transforms[0].cellnum;
	material = default_tform_test_geom.primitives[0].transforms[0].cellmat;
	dx = default_tform_test_geom.primitives[0].transforms[0].dx;
	dy = default_tform_test_geom.primitives[0].transforms[0].dy;
	dz = default_tform_test_geom.primitives[0].transforms[0].dz;
	theta = default_tform_test_geom.primitives[0].transforms[0].theta;
	phi = default_tform_test_geom.primitives[0].transforms[0].phi;
	n_transforms = default_tform_test_geom.n_transforms;
	transform_vec_len = default_tform_test_geom.primitives[0].transforms.size();
	EXPECT_EQ(0u, cellnum);
	EXPECT_EQ(0u, material);
	EXPECT_FLOAT_EQ(0, dx);
	EXPECT_FLOAT_EQ(0, dy);
	EXPECT_FLOAT_EQ(0, dz);
	EXPECT_FLOAT_EQ(0, theta);
	EXPECT_FLOAT_EQ(0, phi);
	EXPECT_EQ(1u, n_transforms);
	EXPECT_EQ(1, transform_vec_len);
}

/**
 * \brief add cell-numbered transform test
 *
 * \details creates a wgeometry object and adds a default primitive and a cell-numbered transform to it. 
 * tests that the cell number, transform coordinates, and transform angles are all equal to the 
 * user-specified values. tests that the number of transforms is incremented and that the transform is 
 * added to the transform vector. tests that the material number defaults to zero.
 */
TEST_F(wgeometryTest, NumValAddTransform)
{
	wgeometry num_val_tform_test_geom;
	num_val_tform_test_geom.add_primitive();
	n_transforms = num_val_tform_test_geom.n_transforms;
	transform_vec_len = num_val_tform_test_geom.primitives[0].transforms.size();
	EXPECT_EQ(0u, n_transforms);
	EXPECT_EQ(0, transform_vec_len);
	num_val_tform_test_geom.add_transform(0,1,2,2,2,3.14159,3.14159);
	num_val_tform_test_geom.update();
	cellnum = num_val_tform_test_geom.primitives[0].transforms[0].cellnum;
	material = num_val_tform_test_geom.primitives[0].transforms[0].cellmat;
	dx = num_val_tform_test_geom.primitives[0].transforms[0].dx;
	dy = num_val_tform_test_geom.primitives[0].transforms[0].dy;
	dz = num_val_tform_test_geom.primitives[0].transforms[0].dz;
	theta = num_val_tform_test_geom.primitives[0].transforms[0].theta;
	phi = num_val_tform_test_geom.primitives[0].transforms[0].phi;
	n_transforms = num_val_tform_test_geom.n_transforms;
	transform_vec_len = num_val_tform_test_geom.primitives[0].transforms.size();
	EXPECT_EQ(1u, cellnum);
	EXPECT_EQ(0u, material);
	EXPECT_FLOAT_EQ(2, dx);
	EXPECT_FLOAT_EQ(2, dy);
	EXPECT_FLOAT_EQ(2, dz);
	EXPECT_FLOAT_EQ(3.14159, theta);
	EXPECT_FLOAT_EQ(3.14159, phi);
	EXPECT_EQ(1u, n_transforms);
	EXPECT_EQ(1, transform_vec_len);
}

/**
 * \brief add cell- and material-numbered transform test
 *
 * \details creates a wgeometry object and adds a default primitive and a cell-numbered transform to it. 
 * tests that the cell number, material, transform coordinates, and transform angles are all equal to the 
 * user-specified values. tests that the number of transforms is incremented and that the transform is 
 * added to the transform vector.
 */
TEST_F(wgeometryTest, NumMatAddTransform)
{
	wgeometry num_mat_tform_test_geom;
	num_mat_tform_test_geom.add_primitive();
	n_transforms = num_mat_tform_test_geom.n_transforms;
	transform_vec_len = num_mat_tform_test_geom.primitives[0].transforms.size();
	EXPECT_EQ(0u, n_transforms);
	EXPECT_EQ(0, transform_vec_len);
	num_mat_tform_test_geom.add_transform(0,1,1,2,2,2,3.14159,3.14159);
	num_mat_tform_test_geom.update();
	cellnum = num_mat_tform_test_geom.primitives[0].transforms[0].cellnum;
	material = num_mat_tform_test_geom.primitives[0].transforms[0].cellmat;
	dx = num_mat_tform_test_geom.primitives[0].transforms[0].dx;
	dy = num_mat_tform_test_geom.primitives[0].transforms[0].dy;
	dz = num_mat_tform_test_geom.primitives[0].transforms[0].dz;
	theta = num_mat_tform_test_geom.primitives[0].transforms[0].theta;
	phi = num_mat_tform_test_geom.primitives[0].transforms[0].phi;
	n_transforms = num_mat_tform_test_geom.n_transforms;
	transform_vec_len = num_mat_tform_test_geom.primitives[0].transforms.size();
	EXPECT_EQ(1u, cellnum);
	EXPECT_EQ(1u, material);
	EXPECT_FLOAT_EQ(2, dx);
	EXPECT_FLOAT_EQ(2, dy);
	EXPECT_FLOAT_EQ(2, dz);
	EXPECT_FLOAT_EQ(3.14159, theta);
	EXPECT_FLOAT_EQ(3.14159, phi);
	EXPECT_EQ(1u, n_transforms);
	EXPECT_EQ(1, transform_vec_len);
}

/**
 * \brief make hexagonal array test
 *
 * \details creates a wgeometry object and makes a hexagonal array. tests that the number of transforms 
 * is incremented, that the transform is added to the transform vector, and that the cell number, 
 * material number, transform coordinates, and transform angles are all equal to the user-specified 
 * values.
 */
TEST_F(wgeometryTest, MakeHexArray)
{
	wgeometry hex_array_test_geom;
	hex_array_test_geom.add_primitive();
	n_transforms = hex_array_test_geom.n_transforms;
	transform_vec_len = hex_array_test_geom.primitives[0].transforms.size();
	EXPECT_EQ(0u, n_transforms);
	EXPECT_EQ(0, transform_vec_len);
	hex_array_test_geom.make_hex_array(0,1,0,0,0,0);
	hex_array_test_geom.update();
	n_transforms = hex_array_test_geom.n_transforms;
	transform_vec_len = hex_array_test_geom.primitives[0].transforms.size();
	cellnum = hex_array_test_geom.primitives[0].transforms[0].cellnum;
	material = hex_array_test_geom.primitives[0].transforms[0].cellmat;
	dx = hex_array_test_geom.primitives[0].transforms[0].dx;
	dy = hex_array_test_geom.primitives[0].transforms[0].dy;
	dz = hex_array_test_geom.primitives[0].transforms[0].dz;
	theta = hex_array_test_geom.primitives[0].transforms[0].theta;
	phi = hex_array_test_geom.primitives[0].transforms[0].phi;
	EXPECT_EQ(1u, n_transforms);
	EXPECT_EQ(1, transform_vec_len);
	EXPECT_EQ(0u, cellnum);
	EXPECT_EQ(0u, material);
	EXPECT_FLOAT_EQ(0, dx);
	EXPECT_FLOAT_EQ(0, dy);
	EXPECT_FLOAT_EQ(0, dz);
	EXPECT_FLOAT_EQ(0, theta);
	EXPECT_FLOAT_EQ(0, phi);
}
