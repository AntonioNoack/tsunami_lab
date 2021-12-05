/**
 * @author Antonio Noack
 * @section DESCRIPTION
 * An artificial tsunami in a swimming pool for testing; it is kind of 2d.
 * Notice: The tests for this setup are indirectly included in the test for the TsunamiEvent2d.
 **/ 
#include "ArtificialTsunami2d.h"

// for use of M_PI
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>

tsunami_lab::t_real tsunami_lab::setups::ArtificialTsunami2d::getHeight( t_real, t_real ) const {
  return +100;
}

tsunami_lab::t_real tsunami_lab::setups::ArtificialTsunami2d::getBathymetry( t_real, t_real ) const {
  return -100;
}

tsunami_lab::t_real tsunami_lab::setups::ArtificialTsunami2d::getDisplacement( t_real i_x, t_real i_y ) const {
  t_real l_x = (i_x - 500)/500 + m_offsetX;// i_x always starts at zero, not at -500, then convert them to the [-1,1], as used by the formula
  t_real l_y = (i_y - 500)/500 + m_offsetY;
  if(l_x < -1 || l_x > +1 || l_y < -1 || l_y > +1) return 0;
  t_real l_f = -std::sin(l_x * M_PI);// sin(x+pi) = -sin(x)
  t_real l_g = 1 - l_y * l_y;
  return 5 * l_f * l_g;
}

tsunami_lab::t_real tsunami_lab::setups::ArtificialTsunami2d::getMomentumX( t_real, t_real ) const {
  return 0;
}

tsunami_lab::t_real tsunami_lab::setups::ArtificialTsunami2d::getMomentumY( t_real, t_real ) const {
  return 0;
}
