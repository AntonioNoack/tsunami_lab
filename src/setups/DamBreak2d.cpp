/**
 * @author Antonio Noack
 *
 * @section DESCRIPTION
 * Two-dimensional dam break problem.
 **/
#include "DamBreak2d.h"
#include <iostream>
#include <cmath>

tsunami_lab::t_real tsunami_lab::setups::DamBreak2d::getHeight( t_real i_x, t_real i_y ) const {
  t_real l_deltaX = i_x - m_locationX;
  t_real l_deltaY = i_y - m_locationY;
  t_real l_f = std::min(std::max(std::sqrt(l_deltaX * l_deltaX + l_deltaY * l_deltaY) - m_radius + (t_real) 0.5, (t_real) 0.0), (t_real) 1.0);
  t_real l_fluidHeight = m_heightInner * (1-l_f) + m_heightOuter * l_f;// linear interpolation between inner and outer area
  // reduce height, if there is an obstacle
  bool l_isInObstacle = i_x >= m_obstacle.x0 && i_x < m_obstacle.x1 && i_y >= m_obstacle.y0 && i_y < m_obstacle.y1;
  if( l_isInObstacle ){// the surface height shall be constant
    l_fluidHeight -= (m_obstacle.bathymetryOverride - m_bathymetry);
  }
  return std::max(l_fluidHeight, (t_real) 0.0);
}

tsunami_lab::t_real tsunami_lab::setups::DamBreak2d::getBathymetry( t_real i_x, t_real i_y ) const {
  bool l_isInObstacle = i_x >= m_obstacle.x0 && i_x < m_obstacle.x1 && i_y >= m_obstacle.y0 && i_y < m_obstacle.y1;
  return l_isInObstacle ? m_obstacle.bathymetryOverride : m_bathymetry;
}

tsunami_lab::t_real tsunami_lab::setups::DamBreak2d::getMomentumX( t_real, t_real ) const {
  return 0;
}

tsunami_lab::t_real tsunami_lab::setups::DamBreak2d::getMomentumY( t_real, t_real ) const {
  return 0;
}