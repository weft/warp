#include <limits.h>
#include "whistory_test.h"
#include <vector>

/**
 * \brief construction test
 *
 * \details constructs a history object with a set number of histories and a wgoemetry object.
 */
TEST_F(whistoryTest, Construction)
{
	N = 10000;
	wgeometry geom;
	whistory construct_hist(N, geom);
}

/**
 * \brief set run type test
 *
 * \details creates a wgeometry object that contains a fissile material and constructs a whistory with 
 * that wgeometry. sets the run type of the whistory.
 */
TEST_F(whistoryTest, SetRunType)
{
        N = 10000;
        wgeometry geom;
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
        geom.set_outer_cell(999,1);
        geom.update();
        geom.check();
        whistory hist(N, geom);
	hist.set_run_type("criticality");
}

/**
 * \brief set tally cell test
 *
 * \details creates a wgeometry object that contains a fissile material and constructs a whistory with 
 * that wgeometry. sets the tally cell of the whistory.
 */
TEST_F(whistoryTest, SetTallyCell)
{
        N = 10000;
        wgeometry geom;
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
        geom.set_outer_cell(999,1);
        geom.update();
        geom.check();
        whistory hist(N, geom);
	hist.set_tally_cell(999);
}

/**
 * \brief set run parameters test
 *
 * \details creates a wgeometry object that contains a fissile material and constructs a whistory with 
 * that wgeometry. sets the run parameters (number of cycles to run and cycles to skip) of the whistory.
 */
TEST_F(whistoryTest, SetRunParam)
{
        N = 10000;
        wgeometry geom;
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
        geom.set_outer_cell(999,1);
        geom.update();
        geom.check();
        whistory hist(N, geom);
	hist.set_run_param(40,20);
}

/**
 * \brief set filename test
 *
 * \details creates a wgeometry object that contains a fissile material and constructs a whistory with 
 * that wgeometry. sets the filename for data to be written to.
 */
TEST_F(whistoryTest, SetFilename)
{
        N = 10000;
        wgeometry geom;
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
        geom.set_outer_cell(999,1);
        geom.update();
        geom.check();
        whistory hist(N, geom);
	EXPECT_NO_THROW(hist.set_filename("whistory_test_log.txt"));
}

/**
 * \brief set device test
 *
 * \details creates a wgeometry object that contains a fissile material and constructs a whistory with 
 * that wgeometry. tests that no error is thrown when the device is set.
 */
TEST_F(whistoryTest, SetDevice)
{
        N = 10000;
        wgeometry geom;
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
        geom.set_outer_cell(999,1);
        geom.update();
        geom.check();
        whistory hist(N, geom);
	EXPECT_NO_THROW(hist.set_device(0));
}

/**
 * \brief initialization test
 *
 * \details creates a wgeometry object that contains a fissile material and constructs a whistory with 
 * that wgeometry. initializes the whistory.
 */
TEST_F(whistoryTest, Init)
{
	N = 10000;
        wgeometry geom;
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
	geom.set_outer_cell(999,1);
	geom.update();
	geom.check();
        whistory init_hist(N, geom);
	init_hist.init();
}

