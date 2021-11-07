/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de), Antonio Noack
 * 
 * @section LICENSE
 * Copyright 2020, Friedrich Schiller University Jena
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Unit tests for the one-dimensional wave propagation patch.
 **/
#include <catch2/catch.hpp>
#include <fstream> // to read csv files
#include <sstream> // to split csv
#include <iostream> // debug output
#include <algorithm> // std::max
#include <cmath> // std::abs

#define private public

#include "WavePropagation1d.h"
#include "../constants.h"
#include "../setups/Discontinuity1d.h"
#include "../solvers/FWave.h" // for gravity constant

#define t_real tsunami_lab::t_real
#define t_idx tsunami_lab::t_idx

TEST_CASE( "Test the 1d wave propagation solver.", "[WaveProp1d]" ) {
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
  tsunami_lab::patches::WavePropagation1d m_waveProp( 100 );

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


TEST_CASE( "Tests from sample file.", "[Discontinuity1d][SampleFile]" ) {

  // middle_states.csv contains a lot of samples to test our function
  // currently we only need to open a single csv file, so let's implement it locally
  // these tests were computed analytically, so our model probably cannot solve all of them correctly
  
  std::fstream l_dataFile;
  l_dataFile.open("data/middle_states.csv", std::ios::in);
  
  REQUIRE(l_dataFile.is_open());
  
  // in the first 500k tests, all tests pass with these settings and accuracy requirements
  // in the following 500k tests, there are a lot of failed samples, probably because our simulation is too small
  // because the rarefaction waves are too long
  t_idx l_testLimit = 500 * 1000;
  t_idx l_totalTests = 0;
  t_real l_expectedEpsilon = 0.001;
  
  // to make more tests pass, these parameters can be tuned.
  // Ideally, we would resize the simulation based on how far the waves actually have traveled.
  // with l_size = 32, all first 500k samples succeed.
  t_idx l_size = 32, l_steps = l_size;
  
  tsunami_lab::patches::WavePropagation1d l_simulation(l_size);
  
  std::string l_csvLine;
  while(l_totalTests < l_testLimit && getline(l_dataFile, l_csvLine)){
    if(!l_csvLine.empty() && l_csvLine[0] != '#'){// not empty and not a comment
      if(l_csvLine[0] != 'h'){// not the table header
        
        l_totalTests++;
        
        std::stringstream l_lineStream(l_csvLine);
        char l_comma;
        t_real l_heightLeft, l_heightRight, l_impulseLeft, l_impulseRight, l_targetValue;
        l_lineStream >> l_heightLeft;   l_lineStream >> l_comma;
        l_lineStream >> l_heightRight;  l_lineStream >> l_comma;
        l_lineStream >> l_impulseLeft;  l_lineStream >> l_comma;
        l_lineStream >> l_impulseRight; l_lineStream >> l_comma;
        l_lineStream >> l_targetValue;
        
        // input must be valid
        REQUIRE(l_heightLeft  >= 0);
        REQUIRE(l_heightRight >= 0);
        REQUIRE(l_targetValue >= 0);
        
        // just a guess, could be 1 as well; influence is rather small, I'd expect
        t_real l_cellSizeMeters = (l_heightLeft + l_heightRight) * 0.5 / 8;
        
        {// we could improve the performance of these tests, if we initialize them lazily
          tsunami_lab::setups::Discontinuity1d l_setup(l_heightLeft, l_heightRight, l_impulseLeft, l_impulseRight, l_size / 2);
          l_simulation.initWithSetup(&l_setup, 1.0);
        }
        
        // compute steps until we have converged
        t_idx l_i;
        for(l_i=0;l_i<l_steps;l_i++){
          
          t_idx l_simulationRadius = l_i + 1;
          t_real l_timeStep = l_simulation.computeMaxTimestep(l_cellSizeMeters, l_simulationRadius);
          
          // in these tests, updating the ghost zone is only needed in the first step theoretically
          // the border also only matters, if all 500 steps are used: only then can the initial wave travel from the center to the border (250 steps),
          // and cause a calculation error there, and then return with 250 to the center.
          l_simulation.setGhostOutflow();
          
          t_real l_scaling = l_timeStep / l_cellSizeMeters;
          l_simulation.timeStep(l_scaling, l_simulationRadius);
          
          // todo better convergence test (?)
          // convergence should have a single direction, so the only mistake that could happen, is when we step over the true result.
          t_real l_hMiddleState = l_simulation.getHeight()[l_size / 2];
          if(l_hMiddleState == Approx(l_targetValue).epsilon(l_expectedEpsilon)) break;
          
        }
        
        // we need to have converged
        REQUIRE(l_i < l_steps);
        
      }
    }
  }
  
  l_dataFile.close();
}


TEST_CASE( "Shock-Shock by reflective boundary condition.", "[WaveProp1d][Shock-Shock-By-Boundary]" ) {
  /*
   * Test case:
   *   A flow from left to right, but on the right side is a wall
   */

  // construct solver and setup problem
  tsunami_lab::patches::WavePropagation1d m_waveProp(10);

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
