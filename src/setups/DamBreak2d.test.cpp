/**
 * @author Antonio Noack
 *
 * @section DESCRIPTION
 * Tests the 2d, circular dam break setup.
 **/
#include <catch2/catch.hpp>
#include "DamBreak2d.h"

TEST_CASE( "Test the two-dimensional dam break setup.", "[DamBreak2d]" ) {
  
  tsunami_lab::setups::DamBreak2d l_setup( 10, 20, 0, 0, 3 );

  // inside
  REQUIRE( l_setup.getHeight(    2, 0 ) == 10 );

  REQUIRE( l_setup.getMomentumX( 2, 0 ) == 0 );

  REQUIRE( l_setup.getMomentumY( 2, 0 ) == 0 );

  REQUIRE( l_setup.getHeight(    2, 1 ) == 10 );

  REQUIRE( l_setup.getMomentumX( 2, 1 ) == 0 );

  REQUIRE( l_setup.getMomentumY( 2, 1 ) == 0 );

  // outside
  REQUIRE( l_setup.getHeight(    4, 0 ) == 20 );

  REQUIRE( l_setup.getMomentumX( 4, 0 ) == 0 );

  REQUIRE( l_setup.getMomentumY( 4, 0 ) == 0 );

  REQUIRE( l_setup.getHeight(    2.5, 2.5 ) == 20 );

  REQUIRE( l_setup.getMomentumX( 2.5, 2.5 ) == 0 );

  REQUIRE( l_setup.getMomentumY( 2.5, 2.5 ) == 0 );  
}
