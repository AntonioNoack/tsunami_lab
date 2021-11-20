/**
 * @author Antonio Noack
 * @section DESCRIPTION
 * Tests the discontinuity setup.
 **/
#include <catch2/catch.hpp>
#include "TsunamiEvent1d.h"
#include "../constants.h"

#define t_real tsunami_lab::t_real

TEST_CASE( "Tests for 1d tsunami event.", "[TsunamiEvent1d]" ) {
  
  t_real l_displacement = 10;
  t_real l_scale = 1;
  
  t_real l_data[12] = { -10, -1, 10, 20, 30, 40, -30, -30, -30, -30, -30, -30 };
  
  // scale = 1, offset = 0, displacement from 6 to 11, in the second half
  tsunami_lab::setups::TsunamiEvent1d l_setup( l_data, 12, l_scale, 6, 11, l_displacement );
  
  REQUIRE(l_setup.getBathymetry( 0, 0) == -20);// coast
  REQUIRE(l_setup.getBathymetry( 1, 0) == -20);// coast
  REQUIRE(l_setup.getBathymetry( 2, 0) == +20);// mountain
  REQUIRE(l_setup.getBathymetry( 3, 0) == +20);
  REQUIRE(l_setup.getBathymetry( 4, 0) == +30);
  REQUIRE(l_setup.getBathymetry( 5, 0) == +40);// end of mountain
  REQUIRE(l_setup.getBathymetry( 6, 0) == -30);// start of pool
  REQUIRE(l_setup.getBathymetry( 7, 0) == -30);
  REQUIRE(l_setup.getBathymetry( 8, 0) == -30);
  REQUIRE(l_setup.getBathymetry( 9, 0) == -30);
  REQUIRE(l_setup.getBathymetry(10, 0) == -30);
  REQUIRE(l_setup.getBathymetry(11, 0) == -30);// end of pool
  
  // check displacement inside the pool
  REQUIRE(l_setup.getDisplacement( 6, 0) == Approx(0));
  REQUIRE(l_setup.getDisplacement( 7, 0) == Approx(+l_displacement * std::sin(  72 * M_PI / 180 )));
  REQUIRE(l_setup.getDisplacement( 8, 0) == Approx(+l_displacement * std::sin( 144 * M_PI / 180 )));
  REQUIRE(l_setup.getDisplacement( 9, 0) == Approx(-l_displacement * std::sin( 144 * M_PI / 180 )));
  REQUIRE(l_setup.getDisplacement(10, 0) == Approx(-l_displacement * std::sin(  72 * M_PI / 180 )));
  REQUIRE(l_setup.getDisplacement(11, 0) == Approx(0));
  
  REQUIRE(l_setup.getBathymetry( 3.5, 0) == +25); // interpolation test
  
  REQUIRE(l_setup.getHeight( 0, 0) == +20);// coast
  REQUIRE(l_setup.getHeight( 1, 0) == +20);// coast
  REQUIRE(l_setup.getHeight( 2, 0) ==   0);// mountain
  REQUIRE(l_setup.getHeight( 3, 0) ==   0);
  REQUIRE(l_setup.getHeight( 4, 0) ==   0);
  REQUIRE(l_setup.getHeight( 5, 0) ==   0);// end of mountain
  REQUIRE(l_setup.getHeight( 6, 0) == +30);// start of pool
  REQUIRE(l_setup.getHeight( 7, 0) == +30);
  REQUIRE(l_setup.getHeight( 8, 0) == +30);
  REQUIRE(l_setup.getHeight( 9, 0) == +30);
  REQUIRE(l_setup.getHeight(10, 0) == +30);
  REQUIRE(l_setup.getHeight(11, 0) == +30);// end of pool
  
}
