#ifndef WGEOMETRY_TEST_H
#define WGEOMETRY_TEST_H
#include "gtest/gtest.h"
#include "wgeometry.h"

/**
 * \brief wgeometry test class
 */
class wgeometryTest : public ::testing::Test {
 public:
   unsigned n_box; /**< number of boxes */
   unsigned n_cyl; /**< number of cylinders */
   unsigned n_hex; /**< number of hexagons */
   unsigned n_sph; /**< number of spheres */
   unsigned n_primitives; /**< number of primitives */
   unsigned n_transforms; /**< number of transforms */
   unsigned outer_cell; /**< outermost cell (usually used for tallying) */
   unsigned n_materials; /**< number of materials */
   unsigned n_isotopes; /**< number of isotopes */
   unsigned fissile_flag; /**< indicates whether or not a material is fissile */
   int transform_vec_len; /**< length of transform vector */
   int prim_vec_len; /**< length of primitive vector */
   unsigned prim_type; /**< primitive type: 0 = box, 1 = cylinder, 2 = hexagon, 3 = sphere */
   int prim_id; /**< primitive ID number */
   unsigned material; /**< material number */
   unsigned mincell; /**< minimum (usually innermost) cell */
   unsigned maxcell; /**< maximum (usually outermost) cell */
   unsigned cellnum; /**< cell number */
   unsigned mat_num; /**< material number */
   int mat_index; /**< material index */
   int mat_vec_len; /**< length of material vector */
   unsigned is_fissile; /**< indicates whether or not a material is fissile */
   int check; /**< value returned by the check function */
   float density; /**< material density */
   float dx; /**< x-coordinate of transform */ 
   float dy; /**< y-coordinate of transform */
   float dz; /**< z-coordinate of transform */
   float theta; /**< azimuthal angle of transform */
   float phi; /**< polar angle of transform */
   std::vector<float> min; /**< coordinate minima vector */
   std::vector<float> max; /**< coordinate maxima vector */
   std::vector<float> loc; /**< coordinate location vector */
   std::vector<float> mat_fracs; /**< material fraction(s) vector */
   std::vector<unsigned> isotopes; /**< isotopes vector */
   /**
    * @fn SetUp
    *
    * resizes the coordinate vectors for Cartesian (x,y,z) coordinate input
    */ 
   virtual void SetUp()
   {
	min.resize(3);
	max.resize(3);
	loc.resize(3);
   }
   virtual void TearDown() {}
}; 
#endif
