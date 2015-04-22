#include <limits.h>
#include "primitive_test.h"

/**
 * \brief default primitive construction test
 *
 * \details tests that a default primitive is constructed with no transforms, a box type,
 * an ID number of 0, the material numbered 0, and all coordinates (minima, maxima, and 
 * location) set to zero. tests that the number of primitives is incremented.
 */
TEST_F(primitiveTest, DefaultConstructionTest)
{
	primitive default_construct_test_prim;
	for(int i = 0; i < 3; i++)
	{
		min[i] = default_construct_test_prim.min[i];
		max[i] = default_construct_test_prim.max[i];
		loc[i] = default_construct_test_prim.location[i];
	}
	num_prim = default_construct_test_prim.num_primitives;
	type = default_construct_test_prim.type;
	prim_id = default_construct_test_prim.primitive_id;
	n_transforms = default_construct_test_prim.n_transforms;
	material = default_construct_test_prim.material;
	for(int i = 0; i < 3; i++)
	{
		EXPECT_EQ(0, min[i]);
		EXPECT_EQ(0, max[i]);
		EXPECT_EQ(0, loc[i]);
	}
	EXPECT_EQ(1, num_prim);
        EXPECT_EQ(0, type);
	EXPECT_EQ(0, prim_id);
	EXPECT_EQ(0, n_transforms);
        EXPECT_EQ(0u, material);
}

/**
 * \brief valued primitive construction test
 *
 * \details tests that a primitive is constructed with user-specified type, ID, 
 * material, coordinate extrema, and location. tests that no transforms are created and that the number of primitives is incremented. 
 */
TEST_F(primitiveTest, ValuedConstructionTest)
{
	std::vector<float> min (3);
	std::vector<float> max (3);
	std::vector<float> loc (3);
	min[0] = -1; min[1] = -1; min[2] = -1;
	max[0] =  1; max[1] =  1; max[2] =  1;
	loc[0] =  0; loc[0] =  0; loc[2] =  0; 
	primitive valued_construct_test_prim = primitive(0,0,min,max,loc);
	for(int i = 0; i < 3; i++)
	{
		min[i] = valued_construct_test_prim.min[i];
		max[i] = valued_construct_test_prim.max[i];
		loc[i] = valued_construct_test_prim.location[i];
	}
	num_prim = valued_construct_test_prim.num_primitives;
	type = valued_construct_test_prim.type;
	prim_id = valued_construct_test_prim.primitive_id;
	n_transforms = valued_construct_test_prim.n_transforms;
	material = valued_construct_test_prim.material;
	transform_vec_len = valued_construct_test_prim.transforms.size();
	for(int i = 0; i < 3; i++)
	{
		EXPECT_EQ(-1, min[i]);
		EXPECT_EQ(1, max[i]);
		EXPECT_EQ(0, loc[i]);
	}
	EXPECT_EQ(1, num_prim);
        EXPECT_EQ(0, type);
	EXPECT_EQ(0, prim_id);
	EXPECT_EQ(0, n_transforms);
        EXPECT_EQ(0u, material);
	EXPECT_EQ(0, transform_vec_len);
}

/**
 * \brief default add transform test
 *
 * \details creates a default primitive and adds a default transform to it. tests 
 * that the transform's coordinates, angles, cell number, and material number are 
 * all zero. tests that the number of transforms is incremented and that the transform 
 * is added to the transform vector.
 */
TEST_F(primitiveTest, DefaultAddTransformTest)
{
	primitive default_add_tform_test_prim;
	n_transforms = default_add_tform_test_prim.n_transforms;
	transform_vec_len = default_add_tform_test_prim.transforms.size();
	EXPECT_EQ(0, n_transforms);
	EXPECT_EQ(0, transform_vec_len);
	default_add_tform_test_prim.add_transform();
	cellnum = default_add_tform_test_prim.transforms[0].cellnum;
	material = default_add_tform_test_prim.transforms[0].cellmat;
	dx = default_add_tform_test_prim.transforms[0].dx;
	dy = default_add_tform_test_prim.transforms[0].dy;
	dz = default_add_tform_test_prim.transforms[0].dz;
	theta = default_add_tform_test_prim.transforms[0].theta;
	phi = default_add_tform_test_prim.transforms[0].phi;
	n_transforms = default_add_tform_test_prim.n_transforms;
	transform_vec_len = default_add_tform_test_prim.transforms.size();
	EXPECT_EQ(0u, cellnum);
	EXPECT_EQ(0u, material);
	EXPECT_FLOAT_EQ(0, dx);
	EXPECT_FLOAT_EQ(0, dy);
	EXPECT_FLOAT_EQ(0, dz);
	EXPECT_FLOAT_EQ(0, theta);
	EXPECT_FLOAT_EQ(0, phi);
	EXPECT_EQ(1, n_transforms);
	EXPECT_EQ(1, transform_vec_len);
}

/**
 * \brief add cell-numbered transform test
 *
 * \details creates a default primitive and adds a cell-numbered transform to 
 * it. tests that the cell number, coordinates, and angles are all equal to the 
 * user-specified values. tests that the transform material is the default zero 
 * value. tests that the number of transforms is incremented and that the transform 
 * is added to the transform vector.
 */
TEST_F(primitiveTest, NumValAddTransformTest)
{
        primitive num_val_add_tform_test_prim;
        n_transforms = num_val_add_tform_test_prim.n_transforms;
        transform_vec_len = num_val_add_tform_test_prim.transforms.size();
        EXPECT_EQ(0, n_transforms);
        EXPECT_EQ(0, transform_vec_len);
        num_val_add_tform_test_prim.add_transform(1,2,2,2,3.14159,3.14159);
        cellnum = num_val_add_tform_test_prim.transforms[0].cellnum;
        material = num_val_add_tform_test_prim.transforms[0].cellmat;
        dx = num_val_add_tform_test_prim.transforms[0].dx;
        dy = num_val_add_tform_test_prim.transforms[0].dy;
        dz = num_val_add_tform_test_prim.transforms[0].dz;
        theta = num_val_add_tform_test_prim.transforms[0].theta;
        phi = num_val_add_tform_test_prim.transforms[0].phi;
        n_transforms = num_val_add_tform_test_prim.n_transforms;
        transform_vec_len = num_val_add_tform_test_prim.transforms.size();
        EXPECT_EQ(1u, cellnum);
        EXPECT_EQ(0u, material);
        EXPECT_FLOAT_EQ(2, dx);
        EXPECT_FLOAT_EQ(2, dy);
        EXPECT_FLOAT_EQ(2, dz);
        EXPECT_FLOAT_EQ(3.14159, theta);
        EXPECT_FLOAT_EQ(3.14159, phi);
        EXPECT_EQ(1, n_transforms);
        EXPECT_EQ(1, transform_vec_len);
}

/**
 * \brief add cell- and material-numbered transform test
 *
 * \details creates a default primitive and adds a cell- and material-numbered 
 * transform to it. tests that the cell and material numbers, coordinates, and 
 * angles are equal to the user-specified values. tests that the number of transforms 
 * is incremented and that the transform is added to the transform vector.
 */ 
TEST_F(primitiveTest, NumMatAddTransformTest)
{
        primitive num_mat_val_add_tform_test_prim;
        n_transforms = num_mat_val_add_tform_test_prim.n_transforms;
        transform_vec_len = num_mat_val_add_tform_test_prim.transforms.size();
        EXPECT_EQ(0, n_transforms);
        EXPECT_EQ(0, transform_vec_len);
        num_mat_val_add_tform_test_prim.add_transform(1,1,2,2,2,3.14159,3.14159);
        cellnum = num_mat_val_add_tform_test_prim.transforms[0].cellnum;
        material = num_mat_val_add_tform_test_prim.transforms[0].cellmat;
        dx = num_mat_val_add_tform_test_prim.transforms[0].dx;
        dy = num_mat_val_add_tform_test_prim.transforms[0].dy;
        dz = num_mat_val_add_tform_test_prim.transforms[0].dz;
        theta = num_mat_val_add_tform_test_prim.transforms[0].theta;
        phi = num_mat_val_add_tform_test_prim.transforms[0].phi;
        n_transforms = num_mat_val_add_tform_test_prim.n_transforms;
        transform_vec_len = num_mat_val_add_tform_test_prim.transforms.size();
        EXPECT_EQ(1u, cellnum);
        EXPECT_EQ(1u, material);
        EXPECT_FLOAT_EQ(2, dx);
        EXPECT_FLOAT_EQ(2, dy);
        EXPECT_FLOAT_EQ(2, dz);
        EXPECT_FLOAT_EQ(3.14159, theta);
        EXPECT_FLOAT_EQ(3.14159, phi);
        EXPECT_EQ(1, n_transforms);
        EXPECT_EQ(1, transform_vec_len);
}

/**
 * \brief make hex array test
 *
 * \details creates a default primitve and adds a hexagonal array transform to 
 * it. tests that the cell and material numbers, coordinates, and angles are 
 * all equal to the user-specified values. tests that the number of transforms 
 * is incremented and that the transform is added to the transform vector.
 */
TEST_F(primitiveTest, HexArrayTest)
{
	primitive hex_array_test_prim;
	n_transforms = hex_array_test_prim.n_transforms;
	transform_vec_len = hex_array_test_prim.transforms.size();
	EXPECT_EQ(0, n_transforms);
	EXPECT_EQ(0, transform_vec_len);
	hex_array_test_prim.make_hex_array(1,0,0,0,0);
	n_transforms = hex_array_test_prim.n_transforms;
	transform_vec_len = hex_array_test_prim.transforms.size();
	cellnum = hex_array_test_prim.transforms[0].cellnum;
	material = hex_array_test_prim.transforms[0].cellmat;
	dx = hex_array_test_prim.transforms[0].dx;
	dy = hex_array_test_prim.transforms[0].dy;
	dz = hex_array_test_prim.transforms[0].dz;
	theta = hex_array_test_prim.transforms[0].theta;
	phi = hex_array_test_prim.transforms[0].phi;
	EXPECT_EQ(1, n_transforms);
	EXPECT_EQ(1, transform_vec_len);
	EXPECT_EQ(0u, cellnum);
	EXPECT_EQ(0u, material);
	EXPECT_FLOAT_EQ(0, dx);
	EXPECT_FLOAT_EQ(0, dy);
	EXPECT_FLOAT_EQ(0, dz);
	EXPECT_FLOAT_EQ(0, theta);
	EXPECT_FLOAT_EQ(0, phi);
}

/**
 * \brief everything test
 *
 * \details executes all primitive functions at least once, testing all of the 
 * parameters affected by the call(s) of the functions.
 */
TEST_F(primitiveTest, EverythingTest)
{
	primitive all_test_prim;
	for(int i = 0; i < 3; i++)
        {
                min[i] = all_test_prim.min[i];
                max[i] = all_test_prim.max[i];
                loc[i] = all_test_prim.location[i];
        }
        num_prim = all_test_prim.num_primitives;
        type = all_test_prim.type;
        prim_id = all_test_prim.primitive_id;
        n_transforms = all_test_prim.n_transforms;
        material = all_test_prim.material;
        for(int i = 0; i < 3; i++)
        {
                EXPECT_EQ(0, min[i]);
                EXPECT_EQ(0, max[i]);
                EXPECT_EQ(0, loc[i]);
        }
        EXPECT_EQ(1, num_prim);
        EXPECT_EQ(0, type);
        EXPECT_EQ(0, prim_id);
        EXPECT_EQ(0, n_transforms);
        EXPECT_EQ(0u, material);
	std::vector<float> min (3);
        std::vector<float> max (3);
        std::vector<float> loc (3);
        min[0] = -1; min[1] = -1; min[2] = -1;
        max[0] =  1; max[1] =  1; max[2] =  1;
        loc[0] =  0; loc[0] =  0; loc[2] =  0;
        all_test_prim = primitive(0,0,min,max,loc);
        for(int i = 0; i < 3; i++)
        {
                min[i] = all_test_prim.min[i];
                max[i] = all_test_prim.max[i];
                loc[i] = all_test_prim.location[i];
        }
        num_prim = all_test_prim.num_primitives;
        type = all_test_prim.type;
        prim_id = all_test_prim.primitive_id;
        n_transforms = all_test_prim.n_transforms;
        material = all_test_prim.material;
        for(int i = 0; i < 3; i++)
        {
                EXPECT_EQ(-1, min[i]);
                EXPECT_EQ(1, max[i]);
                EXPECT_EQ(0, loc[i]);
        }
        EXPECT_EQ(1, num_prim);
        EXPECT_EQ(0, type);
        EXPECT_EQ(1, prim_id);
        EXPECT_EQ(0, n_transforms);
        EXPECT_EQ(0u, material);
        all_test_prim.add_transform(1,2,2,2,3.14159,3.14159);
        cellnum = all_test_prim.transforms[0].cellnum;
        material = all_test_prim.transforms[0].cellmat;
        dx = all_test_prim.transforms[0].dx;
        dy = all_test_prim.transforms[0].dy;
        dz = all_test_prim.transforms[0].dz;
        theta = all_test_prim.transforms[0].theta;
        phi = all_test_prim.transforms[0].phi;
        n_transforms = all_test_prim.n_transforms;
        transform_vec_len = all_test_prim.transforms.size();
        EXPECT_EQ(1u, cellnum);
        EXPECT_EQ(0u, material);
        EXPECT_FLOAT_EQ(2, dx);
        EXPECT_FLOAT_EQ(2, dy);
        EXPECT_FLOAT_EQ(2, dz);
        EXPECT_FLOAT_EQ(3.14159, theta);
        EXPECT_FLOAT_EQ(3.14159, phi);
        EXPECT_EQ(1, n_transforms);
        EXPECT_EQ(1, transform_vec_len);
        all_test_prim.add_transform(1,1,2,2,2,3.14159,3.14159);
        cellnum = all_test_prim.transforms[1].cellnum;
        material = all_test_prim.transforms[1].cellmat;
        dx = all_test_prim.transforms[1].dx;
        dy = all_test_prim.transforms[1].dy;
        dz = all_test_prim.transforms[1].dz;
        theta = all_test_prim.transforms[1].theta;
        phi = all_test_prim.transforms[1].phi;
        n_transforms = all_test_prim.n_transforms;
        transform_vec_len = all_test_prim.transforms.size();
        EXPECT_EQ(1u, cellnum);
        EXPECT_EQ(1u, material);
        EXPECT_FLOAT_EQ(2, dx);
        EXPECT_FLOAT_EQ(2, dy);
        EXPECT_FLOAT_EQ(2, dz);
        EXPECT_FLOAT_EQ(3.14159, theta);
        EXPECT_FLOAT_EQ(3.14159, phi);
        EXPECT_EQ(2, n_transforms);
        EXPECT_EQ(2, transform_vec_len);
        all_test_prim.make_hex_array(1,0,0,0,0);
        n_transforms = all_test_prim.n_transforms;
        transform_vec_len = all_test_prim.transforms.size();
        cellnum = all_test_prim.transforms[2].cellnum;
        material = all_test_prim.transforms[2].cellmat;
        dx = all_test_prim.transforms[2].dx;
        dy = all_test_prim.transforms[2].dy;
        dz = all_test_prim.transforms[2].dz;
        theta = all_test_prim.transforms[2].theta;
        phi = all_test_prim.transforms[2].phi;
        EXPECT_EQ(3, n_transforms);
        EXPECT_EQ(3, transform_vec_len);
        EXPECT_EQ(0u, cellnum);
        EXPECT_EQ(0u, material);
        EXPECT_FLOAT_EQ(0, dx);
        EXPECT_FLOAT_EQ(0, dy);
        EXPECT_FLOAT_EQ(0, dz);
        EXPECT_FLOAT_EQ(0, theta);
        EXPECT_FLOAT_EQ(0, phi);
        all_test_prim.add_transform();
        cellnum = all_test_prim.transforms[3].cellnum;
        material = all_test_prim.transforms[3].cellmat;
        dx = all_test_prim.transforms[3].dx;
        dy = all_test_prim.transforms[3].dy;
        dz = all_test_prim.transforms[3].dz;
        theta = all_test_prim.transforms[3].theta;
        phi = all_test_prim.transforms[3].phi;
        n_transforms = all_test_prim.n_transforms;
        transform_vec_len = all_test_prim.transforms.size();
        EXPECT_EQ(0u, cellnum);
        EXPECT_EQ(0u, material);
        EXPECT_FLOAT_EQ(0, dx);
        EXPECT_FLOAT_EQ(0, dy);
        EXPECT_FLOAT_EQ(0, dz);
        EXPECT_FLOAT_EQ(0, theta);
        EXPECT_FLOAT_EQ(0, phi);
        EXPECT_EQ(4, n_transforms);
        EXPECT_EQ(4, transform_vec_len);
}
