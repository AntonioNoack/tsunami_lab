/**
 * @author Antonio Noack
 * @section DESCRIPTION
 * One-dimensional single discontinuity problem.
 **/
#include "Discontinuity1d.h"

tsunami_lab::t_real tsunami_lab::setups::Discontinuity1d::getHeight( t_real i_x, t_real ) const {
  return i_x < m_locationSplit ? m_heightLeft : m_heightRight;
}

tsunami_lab::t_real tsunami_lab::setups::Discontinuity1d::getBathymetry( t_real i_x, t_real ) const {
  return i_x < m_locationSplit ? m_bathymetryLeft : m_bathymetryRight;
}

tsunami_lab::t_real tsunami_lab::setups::Discontinuity1d::getMomentumX( t_real i_x, t_real ) const {
  return i_x < m_locationSplit ? m_impulseLeft : m_impulseRight;
}

tsunami_lab::t_real tsunami_lab::setups::Discontinuity1d::getMomentumY( t_real, t_real ) const {
  return 0;
}