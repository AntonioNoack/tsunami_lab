/**
 * @author Antonio Noack
 * 
 * @section DESCRIPTION
 * One-dimensional subcritical flow example.
 **/
#include "SubcriticalFlow1d.h"

tsunami_lab::t_real tsunami_lab::setups::SubcriticalFlow1d::getHeight( t_real i_x, t_real ) const {
  return -getBathymetry( i_x, 0 );
}

tsunami_lab::t_real tsunami_lab::setups::SubcriticalFlow1d::getBathymetry( t_real i_x, t_real ) const {
  return i_x >= 8 && i_x <= 12 ? -1.8 - 0.05 * (i_x-10) * (i_x-10) : -2;
}

tsunami_lab::t_real tsunami_lab::setups::SubcriticalFlow1d::getMomentumX( t_real, t_real ) const {
  return 4.42;
}

tsunami_lab::t_real tsunami_lab::setups::SubcriticalFlow1d::getMomentumY( t_real, t_real ) const {
  return 0;
}