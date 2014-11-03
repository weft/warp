#ifndef OPTIX_STUFF_TEST_H
#define OPTIX_STUFF_TEST_H
#include "gtest/gtest.h"
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
 * \brief OptiX stuff test
 */
class optixStuffTest : public ::testing::Test {
 public:
   unsigned N; /**< number of histories */
   unsigned mult; /**< stack size multiplier */
   unsigned compute_device; /**< which GPU to use (always 0) */
   unsigned mincell; /**< minimum cell */
   unsigned maxcell; /**< maximum cell */
   unsigned outer_cell; /**< outermost cell (usually used for tallying) */
   unsigned n_materials; /**< number of materials */
   unsigned material; /**< material number */
   unsigned is_fissile; /**< indicates whether or not a material is fissile */
   unsigned n_isotopes; /**< number of isotopes */
   unsigned index; /**< incrementation variable */
   float density; /**< material density */
   float time; /**< runtime */
   float mu; /**< direction angle */
   float theta; /**< rotation angle */
   float pi; /**< constant value of 3.14159 */
   float x_min; /**< minimum x-coordinate */
   float y_min; /**< minimum y-coordinate */
   float z_min; /**< minimum z-coordinate */
   float x_max; /**< maximum x-coordinate */
   float y_max; /**< maximum y-coordinate */
   float z_max; /**< maximum z-coordinate */
   float rand; /**< random number between zero and one */
   int type; /**< primitive type */
   std::vector<float> min; /**< coordinate minima vector */
   std::vector<float> max; /**< coordinate maxima vector */
   std::vector<float> origin; /**< coordinate origin vector */
   std::vector<float> mat_fracs; /**< material fraction(s) vector */
   std::vector<unsigned> isotopes; /**< isotopes vector */
   std::string accel_type; /**< acceleration type */
   /**
    * @fn SetUp
    *
    * resizes the coordinate vectors for Cartesian (x,y,z) coordinate input, sets the vlaue of pi
    */ 
   virtual void SetUp()
   {
	min.resize(3);
	max.resize(3);
	origin.resize(3);
	pi = 3.14159;
   }
   virtual void TearDown() {}
};
#endif
