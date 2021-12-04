/**
 * @author Antonio Noack
 * @section DESCRIPTION
 * Checkpointed event.
 **/ 
#include "CheckPoint.h"

// std::max, std::min
#include <algorithm>

tsunami_lab::t_idx tsunami_lab::setups::CheckPoint::getIndex( t_real i_x, t_real i_y ) const {

  // calculate index in buffer, and clamp the coordinates to the defined range
  // the ghost cells are copied anyways
  t_real l_indexX = std::min(std::max((t_idx) i_x, (t_idx) 0), m_sizeX-1);
  t_real l_indexY = std::min(std::max((t_idx) i_y, (t_idx) 0), m_sizeY-1);
  
  return l_indexX + l_indexY * m_stride;
  
}

tsunami_lab::t_real tsunami_lab::setups::CheckPoint::getHeight( t_real i_x, t_real i_y ) const {
  return m_height[getIndex(i_x, i_y)];
}

tsunami_lab::t_real tsunami_lab::setups::CheckPoint::getBathymetry( t_real i_x, t_real i_y ) const {
  return m_bathymetry[getIndex(i_x, i_y)];
}

tsunami_lab::t_real tsunami_lab::setups::CheckPoint::getDisplacement( t_real, t_real ) const {
  return 0;
}

tsunami_lab::t_real tsunami_lab::setups::CheckPoint::getMomentumX( t_real i_x, t_real i_y ) const {
  return m_momentumX[getIndex(i_x, i_y)];
}

tsunami_lab::t_real tsunami_lab::setups::CheckPoint::getMomentumY( t_real i_x, t_real i_y ) const {
  return m_momentumY[getIndex(i_x, i_y)];
}
