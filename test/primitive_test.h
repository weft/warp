#ifndef PRIMITIVE_TEST_H
#define PRIMITIVE_TEST_H
#include "gtest/gtest.h"
#include "primitive.h"

/**
 * \brief primitive test class
 */
class primitiveTest : public ::testing::Test {
 public:
 float min[3]; /**< coordinate minima array */
 float max[3]; /**< coordinate maxima array */
 float loc[3]; /**< coordinate location array */
 int num_prim; /**< number of primitives */
 int type; /**< primitive type: 0 = box, 1 = cylinder, 2 = hexagon */
 int prim_id; /**< primitive ID number */
 int n_transforms; /**< number of transforms */
 int transform_vec_len; /**< length of transform vector */
 unsigned material; /**< material number */
 unsigned cellnum; /**< cell number */
 float dx; /**< transform x-coordinate */
 float dy; /**< transform y-coordinate */
 float dz; /**< transform z-coordinate */
 float theta; /**< transform azimuthal angle */
 float phi;  /**< transform polar angle */
 virtual void SetUp() {}
 virtual void TearDown() {}
};
#endif
