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
 * One-dimensional wave propagation patch.
 **/
#include <algorithm> // std::max
#include <cmath> // std::sqrt
#include <iostream>
#include <cassert> // assert
#include "WavePropagation1d.h"
#include "../setups/Setup.h"
#include "../solvers/FWave.h"
#include "../solvers/Roe.h"

tsunami_lab::patches::WavePropagation1d::WavePropagation1d( t_idx i_nCells ) {
  m_nCells = i_nCells;

  // allocate memory including a single ghost cell on each side
  t_idx l_cellCount = m_nCells + 2;
  for( unsigned short l_st = 0; l_st < 2; l_st++ ) {
    m_h [l_st] = new t_real[l_cellCount];
    m_hu[l_st] = new t_real[l_cellCount];
  }
  
  m_bathymetry = new t_real[l_cellCount];

  // init to zero
  for( unsigned short l_st = 0; l_st < 2; l_st++ ) {
    for( t_idx l_ce = 0; l_ce < l_cellCount; l_ce++ ) {
      m_h[l_st][l_ce] = 0;
      m_hu[l_st][l_ce] = 0;
    }
  }
  for( t_idx l_ce = 0; l_ce < l_cellCount; l_ce++ ) {
    m_bathymetry[l_ce] = -20;
  }
}

tsunami_lab::patches::WavePropagation1d::WavePropagation1d( t_idx i_nCells, tsunami_lab::setups::Setup* i_setup, t_real i_scale ) {
  m_nCells = i_nCells;

  // allocate memory including a single ghost cell on each side
  t_idx l_cellCount = m_nCells + 2;
  for( unsigned short l_st = 0; l_st < 2; l_st++ ) {
    m_h[l_st]  = new t_real[l_cellCount];
    m_hu[l_st] = new t_real[l_cellCount];
  }
  
  m_bathymetry = new t_real[l_cellCount];

  initWithSetup(i_setup, i_scale);
}

void tsunami_lab::patches::WavePropagation1d::initWithSetup( tsunami_lab::setups::Setup* i_setup, t_real i_scale ) {
  for( unsigned short l_st = 0; l_st < 2; l_st++ ) {
    t_real* l_h  = m_h [l_st];
    t_real* l_hu = m_hu[l_st];
    for( t_idx l_ce = 0; l_ce < m_nCells + 2; l_ce++ ) {
      t_real l_x = (l_ce - (t_real) 0.5) * i_scale;
      l_h [l_ce] = i_setup->getHeight(l_x, 0);
      l_hu[l_ce] = i_setup->getMomentumX(l_x, 0);
    }
  }
  for( t_idx l_ce = 0; l_ce < m_nCells + 2; l_ce++ ) {
    t_real l_x = (l_ce - (t_real) 0.5) * i_scale;
    m_bathymetry[l_ce] = i_setup->getBathymetry(l_x, 0) + i_setup->getDisplacement(l_x, 0);
  }
}

tsunami_lab::patches::WavePropagation1d::~WavePropagation1d() {
  for( unsigned short l_st = 0; l_st < 2; l_st++ ) {
    delete[] m_h[l_st];
    delete[] m_hu[l_st];
  }
  delete[] m_bathymetry;
}

void tsunami_lab::patches::WavePropagation1d::timeStep( t_real i_scaling ) {
  timeStep(i_scaling, m_nCells);
}

void tsunami_lab::patches::WavePropagation1d::timeStep( t_real i_scaling, t_idx i_updateRadius ) {
  
  // pointers to old and new data
  t_real* l_hOld  = m_h[m_step];
  t_real* l_huOld = m_hu[m_step];

  m_step = !m_step;
  t_real* l_hNew  = m_h[m_step];
  t_real* l_huNew = m_hu[m_step];
  
  t_real* l_b = m_bathymetry;
  
  t_idx l_middleIndex = m_nCells/2+1;
  t_idx l_startIndex = (t_idx) std::max((int64_t) (l_middleIndex - i_updateRadius), (int64_t) 0);
  t_idx l_endIndex   = std::min(l_middleIndex + i_updateRadius, (t_idx) m_nCells+1);

  // init new cell quantities
  for( t_idx l_ce = l_startIndex; l_ce <= l_endIndex; l_ce++ ) {
    l_hNew [l_ce] = l_hOld [l_ce];
    l_huNew[l_ce] = l_huOld[l_ce];
  }
  
  const bool l_useFWaveSolver = m_useFWaveSolver;

  // iterate over edges and update with Riemann solutions
  // ja, das ist nicht mit SIMD einfach parallelisierbar,
  // aber bisher spielt die Performance auch noch keine große Rolle
  for( t_idx l_ceL = l_startIndex; l_ceL < l_endIndex; l_ceL++ ) {
      
    // determine right cell-id
    t_idx l_ceR = l_ceL+1;

    // compute net-updates
    t_real l_netUpdatesL[2];
    t_real l_netUpdatesR[2];
    
    t_real l_hL  = l_hOld[l_ceL];
    t_real l_hR  = l_hOld[l_ceR];
    t_real l_huL = l_huOld[l_ceL];
    t_real l_huR = l_huOld[l_ceR];
    t_real l_bL  = l_b[l_ceL], l_bL0 = l_bL;
    t_real l_bR  = l_b[l_ceR], l_bR0 = l_bR;
    
    if(l_bL > 0 || l_bR > 0){
      // one cell is dry -> reflecting boundary condition
      if(l_bR > 0){
        // right cell is dry
        l_hR = l_hL;
        l_bR = l_bL;
        l_huR = -l_huL;
      } else {
        // left cell is dry
        l_hL = l_hR;
        l_bL = l_bR;
        l_huL = -l_huR;
      }
    }
    
    if(l_useFWaveSolver){
      solvers::FWave::netUpdates( l_hL, l_hR, l_huL, l_huR, l_bL, l_bR,
                                  l_netUpdatesL, l_netUpdatesR );
    } else {
      solvers::Roe::netUpdates( l_hL, l_hR, l_huL, l_huR,
                                l_netUpdatesL, l_netUpdatesR );
    }
    
    // update the cells' quantities
    if(l_bL0 <= 0){
      l_hNew [l_ceL] -= i_scaling * l_netUpdatesL[0];
      l_huNew[l_ceL] -= i_scaling * l_netUpdatesL[1];
    } else l_huNew[l_ceL] = l_hNew[l_ceL] = 0;
    
    if(l_bR0 <= 0){
      l_hNew [l_ceR] -= i_scaling * l_netUpdatesR[0];
      l_huNew[l_ceR] -= i_scaling * l_netUpdatesR[1];
    } else l_huNew[l_ceR] = l_hNew[l_ceR] = 0;
  }
}

tsunami_lab::t_real tsunami_lab::patches::WavePropagation1d::computeMaxTimestep( t_real i_cellSizeMeters, t_idx i_updateRadius ){
  
  t_real* l_h = m_h[m_step];
  t_real* l_hu = m_hu[m_step];
  
  t_idx l_middleIndex = m_nCells/2;
  t_idx l_startIndex = (t_idx) std::max((int64_t) (l_middleIndex - i_updateRadius), (int64_t) 0) + 1;
  t_idx l_endIndex   = std::min(l_middleIndex + i_updateRadius, (t_idx) m_nCells) + 1;
  
  t_real l_maxVelocity = 0;
  t_real l_gravity = tsunami_lab::solvers::FWave::m_gravity;
  
  for(t_idx l_x=l_startIndex;l_x<l_endIndex;l_x++){
    t_real l_height = l_h[l_x], l_height0 = l_height;
    if(l_height > 0){
      // worst case consideration for height; alternatively, we could look at the worst case velocity
      l_height = std::max(l_height, l_h[l_x-1]);
      l_height = std::max(l_height, l_h[l_x+1]);
      t_real l_impulse = l_hu[l_x];
      t_real l_velocity = std::abs(l_impulse) / l_height0;
      t_real l_expectedVelocity = l_velocity + std::sqrt(l_gravity * l_height);
      l_maxVelocity = std::max(l_maxVelocity, l_expectedVelocity);
    }
  }
  
  return 0.5 * i_cellSizeMeters / l_maxVelocity;
  
}

tsunami_lab::t_real tsunami_lab::patches::WavePropagation1d::computeMaxTimestep( t_real i_cellSizeMeters ){
  return computeMaxTimestep( i_cellSizeMeters, m_nCells );
}

void tsunami_lab::patches::WavePropagation1d::setGhostOutflow() {

  t_real* l_b = m_bathymetry;
  t_real* l_h = m_h[m_step];
  t_real* l_hu = m_hu[m_step];

  // set left boundary
  l_b [0] = l_b [1];
  l_h [0] = l_h [1];
  l_hu[0] = l_hu[1];

  // set right boundary
  l_b [m_nCells+1] = l_b [m_nCells];
  l_h [m_nCells+1] = l_h [m_nCells];
  l_hu[m_nCells+1] = l_hu[m_nCells];
}
