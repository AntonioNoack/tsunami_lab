/**
 * @author Antonio Noack
 * @section DESCRIPTION
 * Tests the discontinuity setup.
 **/
#include <catch2/catch.hpp>
#include "Discontinuity1d.h"
#include "../constants.h"

#define t_real tsunami_lab::t_real

TEST_CASE( "Tests for 1d discontinuity.", "[Discontinuity1d-setup]" ) {

  t_real l_hLeft = 2, l_hRight = 17, l_huLeft = 1, l_huRight = 5, l_split = 10;
  t_real l_epsilon = 1e-3, l_delta = 10;
  
  tsunami_lab::setups::Discontinuity1d l_setup(l_hLeft, l_hRight, l_huLeft, l_huRight, l_split);  
  
  REQUIRE(l_setup.getHeight(l_split-l_epsilon, 0) == l_hLeft);
  REQUIRE(l_setup.getHeight(l_split-l_delta, 0) == l_hLeft);
  
  REQUIRE(l_setup.getHeight(l_split+l_epsilon, 0) == l_hRight);
  REQUIRE(l_setup.getHeight(l_split+l_delta, 0) == l_hRight);
  
  REQUIRE(l_setup.getMomentumX(l_split-l_epsilon, 0) == l_huLeft);
  REQUIRE(l_setup.getMomentumX(l_split-l_delta, 0) == l_huLeft);
  REQUIRE(l_setup.getMomentumY(l_split-l_epsilon, 0) == 0);
  REQUIRE(l_setup.getMomentumY(l_split-l_delta, 0) == 0);
  
  REQUIRE(l_setup.getMomentumX(l_split+l_epsilon, 0) == l_huRight);
  REQUIRE(l_setup.getMomentumX(l_split+l_delta, 0) == l_huRight);
  REQUIRE(l_setup.getMomentumY(l_split+l_epsilon, 0) == 0);
  REQUIRE(l_setup.getMomentumY(l_split+l_delta, 0) == 0);
  
}
