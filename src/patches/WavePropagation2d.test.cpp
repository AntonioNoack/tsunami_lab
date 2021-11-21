/**
 * @author Antonio Noack
 * 
 * @section DESCRIPTION
 * Unit tests for the two-dimensional wave propagation patch. Re-purposed from 1d currently.
 **/
#include <catch2/catch.hpp>
#include <fstream> // to read csv files
#include <sstream> // to split csv
#include <iostream> // debug output
#include <algorithm> // std::max
#include <cmath> // std::abs

#define private public

#include "WavePropagation2d.h"
#include "../constants.h"
#include "../setups/Discontinuity1d.h"
#include "../solvers/FWave.h" // for gravity constant

#define t_real tsunami_lab::t_real
#define t_idx tsunami_lab::t_idx

TEST_CASE( "Test the 2d wave propagation solver.", "[WaveProp2d]" ) {
  /*
   * Test case:
   *
   *   Single dam break problem between cell 49 and 50.
   *     left | right
   *       10 | 8
   *        0 | 0
   *
   *   Elsewhere steady state.
   *
   * The net-updates at the respective edge are given as
   * (see derivation in Roe solver):
   *    left          | right
   *      9.394671362 | -9.394671362
   *    -88.25985     | -88.25985
   */

  // construct solver and setup a dambreak problem
  tsunami_lab::patches::WavePropagation2d m_waveProp( 100, 1 );

  for( std::size_t l_ce = 0; l_ce < 50; l_ce++ ) {
    m_waveProp.setHeight( l_ce, 0, 10 );
    m_waveProp.setMomentumX( l_ce, 0, 0 );
  }
  for( std::size_t l_ce = 50; l_ce < 100; l_ce++ ) {
    m_waveProp.setHeight( l_ce, 0, 8 );
    m_waveProp.setMomentumX( l_ce, 0, 0 );
  }

  // set outflow boundary condition
  m_waveProp.setGhostOutflow();

  // perform a time step
  m_waveProp.timeStep( 0.1 );

  // steady state
  for( std::size_t l_ce = 0; l_ce < 49; l_ce++ ) {
    REQUIRE( m_waveProp.getHeight()[l_ce]    == Approx(10) );
    REQUIRE( m_waveProp.getMomentumX()[l_ce] == Approx(0) );
  }

  // dam-break
  // margin added for support of different gravity constants (9.81 vs ~9.8066)
  REQUIRE( m_waveProp.getHeight()[49]    == Approx(10 - 0.1 * 9.394671362).margin(0.001) );
  REQUIRE( m_waveProp.getMomentumX()[49] == Approx( 0 + 0.1 * 88.25985).margin(0.01) );

  REQUIRE( m_waveProp.getHeight()[50]    == Approx(8 + 0.1 * 9.394671362).margin(0.001) );
  REQUIRE( m_waveProp.getMomentumX()[50] == Approx(0 + 0.1 * 88.25985).margin(0.01) );

  // steady state
  for( std::size_t l_ce = 51; l_ce < 100; l_ce++ ) {
    REQUIRE( m_waveProp.getHeight()[l_ce]    == Approx(8) );
    REQUIRE( m_waveProp.getMomentumX()[l_ce] == Approx(0) );
  }
}

TEST_CASE( "Shock-Shock by reflective boundary condition in pseudo-2d.", "[WaveProp2d][Shock-Shock-By-Boundary]" ) {
  /*
   * Test case:
   *   A flow from left to right, but on the right side is a wall
   */

  // construct solver and setup problem
  tsunami_lab::patches::WavePropagation2d m_waveProp( 10, 1 );

  for( std::size_t l_ce = 0; l_ce < 10; l_ce++ ) {
    bool l_isWater = l_ce < 9;
    m_waveProp.setHeight( l_ce, 0, l_isWater ? 10 : 0 );
    m_waveProp.setMomentumX( l_ce, 0, 5 );
    m_waveProp.setBathymetry( l_ce, 0, l_isWater ? -20 : +20 );// floor / wall
  }

  // set outflow boundary condition
  m_waveProp.setGhostOutflow();

  // perform a series of time steps
  for(int i=0;i<50;i++) m_waveProp.timeStep( m_waveProp.computeMaxTimestep(1.0) );
  
  // then everything should be higher than 10, because the shock-shock wave increases the height
  // also the momentum should be much lower, close to zero
  for( std::size_t l_ce = 0; l_ce < 9; l_ce++ ){
    REQUIRE( m_waveProp.getHeight()[l_ce] > 10 );
    REQUIRE( m_waveProp.getMomentumX()[l_ce] < 1 );
  }

}
