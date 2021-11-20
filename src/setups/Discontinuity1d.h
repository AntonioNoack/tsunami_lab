/**
 * @author Antonio Noack
 * @section DESCRIPTION
 * One-dimensional single discontinuity problem.
 **/
#ifndef TSUNAMI_LAB_SETUPS_DISCONTINUITY_1D_H
#define TSUNAMI_LAB_SETUPS_DISCONTINUITY_1D_H

#include "Setup.h"

namespace tsunami_lab {
  namespace setups {
    class Discontinuity1d;
  }
}

/**
 * 1d dam break setup.
 **/
class tsunami_lab::setups::Discontinuity1d: public Setup {
  private:
    
    //! fluid height
    t_real m_heightLeft = 0, m_heightRight = 0;
    
    //! initial impulse
    t_real m_impulseLeft = 0, m_impulseRight = 0;
    
    //! initial bathymetry
    t_real m_bathymetryLeft = -20, m_bathymetryRight = -20;
    
    //! location of the split
    t_real m_locationSplit = 0;
  
  public:
    /**
     * Constructor.
     *
     * @param i_heightLeft water height on the left side of the discontinuity.
     * @param i_heightRight water height on the right side of the discontinuity.
     * @param i_impulseLeft water momentum in x direction on the left side of the discontinuity.
     * @param i_impulseRight water momentum in x direction on the right side of the discontinuity.
     * @param i_locationSplit location (x-coordinate) of the discontinuity.
     **/
    explicit Discontinuity1d( 
                    t_real i_heightLeft,  t_real i_heightRight,
                    t_real i_impulseLeft, t_real i_impulseRight,
                    t_real i_locationSplit ) : 
                    m_heightLeft(i_heightLeft), m_heightRight(i_heightRight),
                    m_impulseLeft(i_impulseLeft), m_impulseRight(i_impulseRight),
                    m_locationSplit(i_locationSplit) {}
    /**
     * Constructor.
     *
     * @param i_heightLeft water height on the left side of the discontinuity.
     * @param i_heightRight water height on the right side of the discontinuity.
     * @param i_impulseLeft water momentum in x direction on the left side of the discontinuity.
     * @param i_impulseRight water momentum in x direction on the right side of the discontinuity.
     * @param i_bathymetryLeft bathymetry on the left side.
     * @param i_bathymetryRight bathymetry on the right side.
     * @param i_locationSplit location (x-coordinate) of the discontinuity.
     **/
    explicit Discontinuity1d( 
                    t_real i_heightLeft,     t_real i_heightRight,
                    t_real i_impulseLeft,    t_real i_impulseRight,
                    t_real i_bathymetryLeft, t_real i_bathymetryRight,
                    t_real i_locationSplit ) : 
                    m_heightLeft(i_heightLeft), m_heightRight(i_heightRight),
                    m_impulseLeft(i_impulseLeft), m_impulseRight(i_impulseRight),
                    m_bathymetryLeft(i_bathymetryLeft), m_bathymetryRight(i_bathymetryRight),
                    m_locationSplit(i_locationSplit) {}
    
    /**
     * Gets the water height at a given point.
     *
     * @param i_x x-coordinate of the queried point.
     * @return height at the given point.
     **/
    t_real getHeight( t_real i_x, t_real ) const;

    /**
     * Gets the water depth at a given point.
     * Positive values mean above sea level, negative values mean below sea level.
     *
     * @param i_x x-coordinate of the queried point.
     * @return water depth at the given point.
     **/
    t_real getBathymetry( t_real i_x, t_real ) const;
    
    /**
     * Gets the momentum in x-direction.
     *
     * @param i_x x-coordinate of the queried point.
     * @return momentum in x-direction.
     **/
    t_real getMomentumX( t_real i_x, t_real ) const;

    /**
     * Gets the momentum in y-direction.
     *
     * @return momentum in y-direction.
     **/
    t_real getMomentumY( t_real, t_real ) const;

    /**
     * Gets the earth quake displacement.
     * Positive values mean the ground moved upwards.
     *
     * @return displacement in meters at the given point.
     **/
    t_real getDisplacement( t_real, t_real ) const {
      return 0;
    }
    
};

#endif // TSUNAMI_LAB_SETUPS_DISCONTINUITY_1D_H