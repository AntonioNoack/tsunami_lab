/**
 * @author Antonio Noack
 *
 * @section DESCRIPTION
 * Two-dimensional dam break problem.
 **/
#include "DamBreak2d.h"
#include <iostream>

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
  t_real l_fluidHeight = l_deltaX * l_deltaX + l_deltaY * l_deltaY < m_radius * m_radius ? m_heightInner : m_heightOuter;
  // reduce height, if there is an obstacle
  bool l_isInObstacle = i_x >= m_obstacle.x0 && i_x < m_obstacle.x1 && i_y >= m_obstacle.y0 && i_y < m_obstacle.y1;
  if( l_isInObstacle ){
    t_real l_defaultBathymetry = getBathymetry(m_obstacle.x0 - 1, 0);
    l_fluidHeight -= (m_obstacle.bathymetryOverride - l_defaultBathymetry);
  }
  return std::max(l_fluidHeight, (t_real) 0);
}

tsunami_lab::t_real tsunami_lab::setups::DamBreak2d::getBathymetry( t_real i_x, t_real i_y ) const {
  bool l_isInObstacle = i_x >= m_obstacle.x0 && i_x < m_obstacle.x1 && i_y >= m_obstacle.y0 && i_y < m_obstacle.y1;
  return l_isInObstacle ? m_obstacle.bathymetryOverride : -20;
}

tsunami_lab::t_real tsunami_lab::setups::DamBreak2d::getMomentumX( t_real, t_real ) const {
  return 0;
}

tsunami_lab::t_real tsunami_lab::setups::DamBreak2d::getMomentumY( t_real, t_real ) const {
  return 0;
}