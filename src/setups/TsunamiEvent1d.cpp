/**
 * @author Antonio Noack
 * @section DESCRIPTION
 * One-dimensional tsunami event.
 **/ 
#include "TsunamiEvent1d.h"

// std::max, std::min
#include <algorithm>

// for use of M_PI
#define _USE_MATH_DEFINES
#include <math.h>

tsunami_lab::t_real tsunami_lab::setups::TsunamiEvent1d::getRawBathymetry( t_real i_x ) const {

  // calculate index in buffer
  t_real l_index = i_x * m_scale + m_offset;
  
  // left side check
  if(l_index <= 0) return m_bathymetry[0];
  
  // find the integer index
  t_idx  l_index0 = (t_idx) l_index;
  
  // right side check
  if(l_index0 >= m_size-1) return m_bathymetry[m_size - 1];

  // find interpolation factor
  t_real l_fract = l_index - l_index0;
  return m_bathymetry[l_index0] * (1-l_fract) + m_bathymetry[l_index0+1] * l_fract;// interpolate
  
}

tsunami_lab::t_real tsunami_lab::setups::TsunamiEvent1d::getHeight( t_real i_x, t_real ) const {
  t_real l_bathIn = getRawBathymetry(i_x);
  t_real l_delta  = m_shoreCliffHeight;
  return l_bathIn < 0 ? std::max(-l_bathIn, l_delta) : 0;
}

tsunami_lab::t_real tsunami_lab::setups::TsunamiEvent1d::getBathymetry( t_real i_x, t_real ) const {
  // if we are on the shore, clamp the value either on top of the cliff or in the cliff-deep water
  t_real l_bathIn = getRawBathymetry(i_x);
  t_real l_delta  = m_shoreCliffHeight;
  t_real l_valueBeforeDisplacement = std::abs(l_bathIn) < l_delta ?
    l_bathIn < 0 ? -l_delta : l_delta :
    l_bathIn;
  // an earth quake happened, so add its displacement
  return l_valueBeforeDisplacement + getDisplacement( i_x, 0 );
}

tsunami_lab::t_real tsunami_lab::setups::TsunamiEvent1d::getDisplacement( t_real i_x, t_real ) const {
  // in the task, a 0 was forgotten
  t_real l_left = m_displacementStart + m_offset;
  t_real l_right = m_displacementEnd + m_offset;
  t_real l_waveHeight = m_displacement;// in the original formula, +PI is added, which is just changing the sign
  t_real l_normalizedPosition = (i_x-l_left)/(l_right-l_left);
  return i_x > l_left && i_x < l_right ? l_waveHeight * std::sin(l_normalizedPosition * (2 * M_PI)) : 0;
}

tsunami_lab::t_real tsunami_lab::setups::TsunamiEvent1d::getMomentumX( t_real, t_real ) const {
  return 0;
}

tsunami_lab::t_real tsunami_lab::setups::TsunamiEvent1d::getMomentumY( t_real, t_real ) const {
  return 0;
}