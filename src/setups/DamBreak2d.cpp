/**
 * @author Antonio Noack
 *
 * @section DESCRIPTION
 * Two-dimensional dam break problem.
 **/
#include "DamBreak2d.h"

tsunami_lab::setups::DamBreak2d::DamBreak2d( t_real i_heightInner,
                                             t_real i_heightOuter,
                                             t_real i_locationX,
                                             t_real i_locationY,
                                             t_real i_radius ) {
  m_heightInner = i_heightInner;
  m_heightOuter = i_heightOuter;
  m_locationX = i_locationX;
  m_locationY = i_locationY;
  m_radius = i_radius;
}

tsunami_lab::t_real tsunami_lab::setups::DamBreak2d::getHeight( t_real i_x, t_real i_y ) const {
  t_real l_deltaX = i_x - m_locationX;
  t_real l_deltaY = i_y - m_locationY;
  return l_deltaX * l_deltaX + l_deltaY * l_deltaY < m_radius * m_radius ? m_heightInner : m_heightOuter;
}

tsunami_lab::t_real tsunami_lab::setups::DamBreak2d::getBathymetry( t_real, t_real ) const {
  return -20;
}

tsunami_lab::t_real tsunami_lab::setups::DamBreak2d::getMomentumX( t_real, t_real ) const {
  return 0;
}

tsunami_lab::t_real tsunami_lab::setups::DamBreak2d::getMomentumY( t_real, t_real ) const {
  return 0;
}