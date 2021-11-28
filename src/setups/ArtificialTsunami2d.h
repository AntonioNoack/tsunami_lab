/**
 * @author Antonio Noack
 * @section DESCRIPTION
 * An artificial tsunami in a swimming pool for testing in 2d
 **/
#ifndef TSUNAMI_LAB_SETUPS_ARTIFICIAL_TSUNAMI_2D_H
#define TSUNAMI_LAB_SETUPS_ARTIFICIAL_TSUNAMI_2D_H

#include "Setup.h"

// assert
#include <cassert>

namespace tsunami_lab {
  namespace setups {
    class ArtificialTsunami2d;
  }
}

/**
 * 2d dam break setup.
 **/
class tsunami_lab::setups::ArtificialTsunami2d: public Setup {
  private:
    // offset of the displacement
    t_real m_offsetX, m_offsetY;
  public:
    /**
     * Constructor.
     **/
    explicit ArtificialTsunami2d(){}
    
    /**
     * Gets the water height at a given point.
     *
     * @return height at the given point.
     **/
    t_real getHeight( t_real, t_real ) const;

    /**
     * Gets the water depth at a given point.
     * Positive values mean above sea level, negative values mean below sea level.
     *
     * @return water depth at the given point.
     **/
    t_real getBathymetry( t_real, t_real ) const;

    /**
     * Gets the earth quake displacement.
     * Positive values mean the ground moved upwards.
     *
     * @param i_x x-coordinate of the queried point.
     * @param i_y y-coordinate of the queried point.
     * @return displacement in meters at the given point.
     **/
    t_real getDisplacement( t_real i_x, t_real i_y ) const;
    
    /**
     * Gets the momentum in x-direction.
     *
     * @return momentum in x-direction.
     **/
    t_real getMomentumX( t_real, t_real ) const;

    /**
     * Gets the momentum in y-direction.
     *
     * @return momentum in y-direction.
     **/
    t_real getMomentumY( t_real, t_real ) const;
    
};

#endif // TSUNAMI_LAB_SETUPS_ARTIFICIAL_TSUNAMI_2D_H