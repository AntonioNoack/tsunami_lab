/**
 * @author Antonio Noack
 * @section DESCRIPTION
 * One-dimensional tsunami event.
 **/ 
#include "TsunamiEvent2d.h"

// std::max, std::min
#include <algorithm>
#include <iostream>

// for use of M_PI
#define _USE_MATH_DEFINES
#include <math.h>

tsunami_lab::t_real tsunami_lab::setups::TsunamiEvent2d::getInterpolatedValue( t_real i_x, t_real i_y, t_idx i_sizeX, t_idx i_sizeY, t_idx i_stride, t_real* i_data ) const {

  // currently assumes m_sizeX, m_sizeY >= 2
  
  // calculate index in buffer
  t_real l_indexX = std::min(std::max(i_x + m_offset, (t_real) 0), (t_real) i_sizeX);
  t_real l_indexY = std::min(std::max(i_y + m_offset, (t_real) 0), (t_real) i_sizeY);
  
  // find the integer index
  t_idx l_index0X = (t_idx) l_indexX;
  t_idx l_index0Y = (t_idx) l_indexY;
  
  // right/bottom side check
  l_index0X = std::min(l_index0X, i_sizeX-2);
  l_index0Y = std::min(l_index0Y, i_sizeY-2);
  
  // find interpolation factor
  t_real l_fractX1 = l_indexX - l_index0X;
  t_real l_fractY1 = l_indexY - l_index0Y;
  
  t_real l_fractX2 = 1-l_fractX1;
  t_real l_fractY2 = 1-l_fractY1;
  
  t_idx l_index0 = l_index0X + l_index0Y * i_stride;
  t_idx l_index1 = l_index0 + i_stride;
  
  // interpolate on both axes
  t_real l_x0 = i_data[l_index0] * l_fractX2 + i_data[l_index0+1] * l_fractX1;
  t_real l_x1 = i_data[l_index1] * l_fractX2 + i_data[l_index1+1] * l_fractX1;
  
  return l_x0 * l_fractY2 + l_x1 * l_fractY1;
  
}

tsunami_lab::t_real tsunami_lab::setups::TsunamiEvent2d::getHeight( t_real i_x, t_real i_y ) const {
  t_real l_bathIn = getInterpolatedValue(i_x * m_scaleBath, i_y * m_scaleBath, m_sizeXBath, m_sizeYBath, m_strideBath, m_bathymetry);
  t_real l_delta  = m_shoreCliffHeight;
  return l_bathIn < 0 ? std::max(-l_bathIn, l_delta) : 0;
}

tsunami_lab::t_real tsunami_lab::setups::TsunamiEvent2d::getBathymetry( t_real i_x, t_real i_y ) const {
  // if we are on the shore, clamp the value either on top of the cliff or in the cliff-deep water
  t_real l_bathIn = getInterpolatedValue(i_x * m_scaleBath, i_y * m_scaleBath, m_sizeXBath, m_sizeYBath, m_strideBath, m_bathymetry);
  t_real l_delta  = m_shoreCliffHeight;
  return std::abs(l_bathIn) < l_delta ?
    l_bathIn < 0 ? -l_delta : l_delta :
    l_bathIn;
}

tsunami_lab::t_real tsunami_lab::setups::TsunamiEvent2d::getDisplacement( t_real i_x, t_real i_y ) const {
  return getInterpolatedValue( i_x * m_scaleDisp, i_y * m_scaleDisp, m_sizeXDisp, m_sizeYDisp, m_strideDisp, m_displacement );
}

tsunami_lab::t_real tsunami_lab::setups::TsunamiEvent2d::getMomentumX( t_real, t_real ) const {
  return 0;
}

tsunami_lab::t_real tsunami_lab::setups::TsunamiEvent2d::getMomentumY( t_real, t_real ) const {
  return 0;
}
