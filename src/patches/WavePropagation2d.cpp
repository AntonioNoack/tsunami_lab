/**
 * @author Antonio Noack
 * 
 * @section DESCRIPTION
 * One-dimensional wave propagation patch.
 **/
#include <algorithm> // std::max
#include <cmath> // std::sqrt
#include <iostream>
#include <cassert> // assert
#include "WavePropagation2d.h"
#include "../setups/Setup.h"
#include "../solvers/FWave.h"
#include "../solvers/Roe.h"

tsunami_lab::patches::WavePropagation2d::WavePropagation2d( t_idx i_nCellsX, t_idx i_nCellsY ) {

  m_nCellsX = i_nCellsX;
  m_nCellsY = i_nCellsY;
  
  m_nCells = (m_nCellsX+2) * (m_nCellsY+2);

  // allocate memory including a single ghost cell on each side
  t_idx l_cellCount = m_nCells;
  for( unsigned short l_st = 0; l_st < 2; l_st++ ) {
    m_h [l_st] = new t_real[l_cellCount];
    m_hu[l_st] = new t_real[l_cellCount];
	m_hv[l_st] = new t_real[l_cellCount];
  }
  
  m_bathymetry = new t_real[l_cellCount];

  // init to zero
  for( unsigned short l_st = 0; l_st < 2; l_st++ ) {
    for( t_idx l_ce = 0; l_ce < l_cellCount; l_ce++ ) {
      m_h [l_st][l_ce] = 0;
      m_hu[l_st][l_ce] = 0;
	  m_hv[l_st][l_ce] = 0;
    }
  }
  for( t_idx l_ce = 0; l_ce < l_cellCount; l_ce++ ) {
    m_bathymetry[l_ce] = -20;
  }
}

tsunami_lab::patches::WavePropagation2d::WavePropagation2d( t_idx i_nCellsX, t_idx i_nCellsY, tsunami_lab::setups::Setup* i_setup, t_real i_scaleX, t_real i_scaleY ) {
  
  m_nCellsX = i_nCellsX;
  m_nCellsY = i_nCellsY;
  
  m_nCells = (m_nCellsX+2) * (m_nCellsY+2);

  // allocate memory including a single ghost cell on all sides
  for( unsigned short l_st = 0; l_st < 2; l_st++ ) {
    m_h [l_st] = new t_real[m_nCells];
    m_hu[l_st] = new t_real[m_nCells];
	m_hv[l_st] = new t_real[m_nCells];
  }
  
  m_bathymetry = new t_real[m_nCells];

  initWithSetup( i_setup, i_scaleX, i_scaleY );
}

void tsunami_lab::patches::WavePropagation2d::initWithSetup( tsunami_lab::setups::Setup* i_setup, t_real i_scaleX, t_real i_scaleY ) {
  for( unsigned short l_st = 0; l_st < 2; l_st++ ) {
    t_real* l_h  = m_h [l_st];
    t_real* l_hu = m_hu[l_st];
    t_real* l_hv = m_hv[l_st];
    for( t_idx l_iy = 0, l_i = 0; l_iy < m_nCellsY + 2; l_iy++ ) {
      t_real l_y = (l_iy - 1) * i_scaleX;
      for( t_idx l_ix = 0; l_ix < m_nCellsX + 2; l_ix++, l_i++ ) {
        t_real l_x = (l_ix - 1) * i_scaleY;
        l_h [l_i] = i_setup->getHeight(    l_x, l_y );
        l_hu[l_i] = i_setup->getMomentumX( l_x, l_y );
        l_hv[l_i] = i_setup->getMomentumY( l_x, l_y );
      }
    }
  }
  
  for( t_idx l_iy = 0, l_i = 0; l_iy < m_nCellsY + 2; l_iy++ ) {
    t_real l_y = (l_iy - 1) * i_scaleX;
    for( t_idx l_ix = 0; l_ix < m_nCellsX + 2; l_ix++, l_i++ ) {
      t_real l_x = (l_ix - 1) * i_scaleY;
      m_bathymetry[l_i] = i_setup->getBathymetry( l_x, l_y );
    }
  }
}

tsunami_lab::patches::WavePropagation2d::~WavePropagation2d() {
  for( unsigned short l_st = 0; l_st < 2; l_st++ ) {
    delete[] m_h [l_st];
    delete[] m_hu[l_st];
    delete[] m_hv[l_st];
  }
  delete[] m_bathymetry;
}

void tsunami_lab::patches::WavePropagation2d::internalUpdate( t_real i_scaling, t_idx l_ceL, t_idx l_ceR, 
    t_real* l_hOld, t_real* l_huOld, t_real* l_hNew, t_real* l_huNew ) {
  
  t_real* l_b = m_bathymetry;
  
  // compute net-updates
  t_real l_netUpdatesL[2];
  t_real l_netUpdatesR[2];
  
  t_real l_hL  = l_hOld [l_ceL];
  t_real l_hR  = l_hOld [l_ceR];
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
  
  if(m_useFWaveSolver){
    solvers::FWave::netUpdates( l_hL, l_hR, l_huL, l_huR, l_bL, l_bR, l_netUpdatesL, l_netUpdatesR );
  } else {
    solvers::Roe::netUpdates( l_hL, l_hR, l_huL, l_huR, l_netUpdatesL, l_netUpdatesR );
  }
  
  // update the cells' quantities
  if(l_bL0 < 0){
    l_hNew [l_ceL] -= i_scaling * l_netUpdatesL[0];
    l_huNew[l_ceL] -= i_scaling * l_netUpdatesL[1];
  } else l_huNew[l_ceL] = l_hNew[l_ceL] = 0;
  
  if(l_bR0 < 0){
    l_hNew [l_ceR] -= i_scaling * l_netUpdatesR[0];
    l_huNew[l_ceR] -= i_scaling * l_netUpdatesR[1];
  } else l_huNew[l_ceR] = l_hNew[l_ceR] = 0;
  
}

void tsunami_lab::patches::WavePropagation2d::timeStep( t_real i_scaling ) {
  
  // pointers to old and new data
  t_real* l_hOld  = m_h [0];
  t_real* l_huOld = m_hu[m_step];
  t_real* l_hvOld = m_hv[m_step];
  
  m_step = !m_step;
  
  t_real* l_hNew  = m_h [1];
  t_real* l_huNew = m_hu[m_step];
  t_real* l_hvNew = m_hv[m_step];

  // init new cell quantities
  for( t_idx l_ce = 0; l_ce < m_nCells; l_ce++ ) {
    l_hNew [l_ce] = l_hOld [l_ce];
    l_huNew[l_ce] = l_huOld[l_ce];
    l_hvNew[l_ce] = l_hvOld[l_ce];
  }
  
    //////////////////////////////
   // half step in x direction //
  //////////////////////////////
  
  t_idx l_stride = getStride();

  // iterate over edges and update with Riemann solutions
  for( t_idx l_iy = 0; l_iy < m_nCellsY + 2; l_iy++ ) {
    for( t_idx l_ix = 0; l_ix < m_nCellsX + 2 - 1; l_ix++ ) {
      // determine cell-ids left and right
      t_idx l_ceL = l_ix + l_iy * l_stride;
      t_idx l_ceR = l_ceL + 1;
      internalUpdate( i_scaling, l_ceL, l_ceR, l_hOld, l_huOld, l_hNew, l_huNew );
    }
  }
  
  // pointers to old and new data
  l_hOld  = l_hNew;
  l_hNew  = m_h [0];

  // init new cell quantities
  for( t_idx l_ce = 0; l_ce < m_nCells; l_ce++ ) {
    l_hNew [l_ce] = l_hOld [l_ce];
  }
  
    //////////////////////////////
   // half step in y direction //
  //////////////////////////////
  
  // iterate over edges and update with Riemann solutions
  for( t_idx l_iy = 0; l_iy < m_nCellsY + 2 - 1; l_iy++ ) {
    for( t_idx l_ix = 0; l_ix < m_nCellsX + 2; l_ix++ ) {
      // determine cell-ids above and below
      t_idx l_ceL = l_ix + l_iy * l_stride;
      t_idx l_ceR = l_ceL + l_stride;
      internalUpdate( i_scaling, l_ceL, l_ceR, l_hOld, l_huOld, l_hNew, l_huNew );
    }
  }
  
}

tsunami_lab::t_real tsunami_lab::patches::WavePropagation2d::computeMaxTimestep( t_real i_cellSizeMeters ){
  
  t_real* l_h  = m_h [0];
  t_real* l_hu = m_hu[0];
  t_real* l_hv = m_hv[0];
  
  t_real l_maxVelocity = 0;
  t_real l_gravity = tsunami_lab::solvers::FWave::m_gravity;
  
  t_idx  l_stride = getStride();
  
  for( t_idx l_iy = 1; l_iy <= m_nCellsY; l_iy++){
    t_idx l_i = l_iy * l_stride + 1;// +1, because we start iterating at l_ix = 1
    for( t_idx l_ix = 1; l_ix <= m_nCellsX; l_ix++, l_i++){
      t_real l_height = l_h[l_i], l_height0 = l_height;
      // worst case consideration for height; alternatively, we could look at the worst case velocity
      l_height = std::max(l_height, l_h[l_i-1]);// left
      l_height = std::max(l_height, l_h[l_i+1]);// right
      l_height = std::max(l_height, l_h[l_i-l_stride]);// top
      l_height = std::max(l_height, l_h[l_i+l_stride]);// bottom
      t_real l_impulse = std::max(std::abs(l_hu[l_i]), std::abs(l_hv[l_i]));
      t_real l_velocity = l_impulse / l_height0;
      t_real l_expectedVelocity = l_velocity + std::sqrt(l_gravity * l_height);
      l_maxVelocity = std::max(l_maxVelocity, l_expectedVelocity);
    }
  }
  
  // 0.45 instead of 0.50, because after the first half step,
  // h/hv may have changed and be incorrect. So be save, and do a smaller step
  
  // in the future, we could cancel the computation and restart, if this step guess failed.
  return 0.45 * i_cellSizeMeters / l_maxVelocity;
  
}

void tsunami_lab::patches::WavePropagation2d::setGhostOutflow() {
  
  t_real* l_b  = m_bathymetry;
  t_real* l_h  = m_h [0];
  t_real* l_hu = m_hu[0];
  t_real* l_hv = m_hv[0];
  
  t_idx l_stride = getStride();
  
  // set left boundary
  for(t_idx l_y = 0; l_y < m_nCellsY+2; l_y++){
	t_idx l_i0 = l_y * l_stride;
	t_idx l_i1 = l_i0 + 1;
    l_b [l_i0] = l_b [l_i1];
    l_h [l_i0] = l_h [l_i1];
    l_hu[l_i0] = l_hu[l_i1];
    l_hv[l_i0] = l_hv[l_i1];
  }
  
  // set right boundary
  for(t_idx l_y = 0; l_y < m_nCellsY+2; l_y++){
	t_idx l_i0 = m_nCellsX + 1 + l_y * l_stride;
	t_idx l_i1 = l_i0 - 1;
    l_b [l_i0] = l_b [l_i1];
    l_h [l_i0] = l_h [l_i1];
    l_hu[l_i0] = l_hu[l_i1];
    l_hv[l_i0] = l_hv[l_i1];
  }
  
  // set top boundary
  for(t_idx l_x = 0; l_x < m_nCellsX+2; l_x++){
	t_idx l_i0 = l_x;
	t_idx l_i1 = l_i0 + l_stride;
    l_b [l_i0] = l_b [l_i1];
    l_h [l_i0] = l_h [l_i1];
    l_hu[l_i0] = l_hu[l_i1];
    l_hv[l_i0] = l_hv[l_i1];
  }
  
  // set bottom boundary
  for(t_idx l_x = 0; l_x < m_nCellsX+2; l_x++){
	t_idx l_i0 = l_x + (m_nCellsY + 1) * l_stride;
	t_idx l_i1 = l_i0 - l_stride;
    l_b [l_i0] = l_b [l_i1];
    l_h [l_i0] = l_h [l_i1];
    l_hu[l_i0] = l_hu[l_i1];
    l_hv[l_i0] = l_hv[l_i1];
  }
  
}
