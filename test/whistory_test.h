#ifndef WHISTORY_TEST_H
#define WHISTORY_TEST_H
#include "gtest/gtest.h"
#include "warp.h"
#include "optix_stuff.h"

/**
 * \brief whistory test class
 */
class whistoryTest : public ::testing::Test {
 public:
   unsigned N; /**< number of histories */
   unsigned material; /**< material number */
   unsigned is_fissile; /**< indicates whether or not a material is fissile */
   unsigned n_isotopes; /**< number of isotopes */
   unsigned type; /**< primitive type: 0 = box, 1 = cylinder, 2 = hexagon */
   float density; /**< material density */
   std::vector<float> min; /**< coordinate minima vector */
   std::vector<float> max; /**< coordinate maxima vector */
   std::vector<float> origin; /**< coordinate origin vector */
   std::vector<float> mat_fracs; /**< material fraction(s) vector */
   std::vector<unsigned> isotopes; /**< isotopes vector */
   /**
    * @fn SetUp
    *
    * resizes the coordinate vectors for Cartesian (x,y,z) coordinate input
    */
   virtual void SetUp() {
	min.resize(3);
	max.resize(3);
	origin.resize(3);
   }
   virtual void TearDown() {}
}; 
#endif
