/**
 * @author Antonio Noack
 *
 * @section DESCRIPTION
 * One-dimensional subcritical flow.
 **/
#ifndef TSUNAMI_LAB_SETUPS_SUBCRITICAL_FLOW_1D_H
#define TSUNAMI_LAB_SETUPS_SUBCRITICAL_FLOW_1D_H

#include "Setup.h"

namespace tsunami_lab {
  namespace setups {
    class SubcriticalFlow1d;
  }
}

/**
 * 1d dam break setup.
 **/
class tsunami_lab::setups::SubcriticalFlow1d: public Setup {
  public:
    /**
     * Constructor.
     **/
    SubcriticalFlow1d() {}

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

#endif