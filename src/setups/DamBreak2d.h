/**
 * @author Antonio Noack
 *
 * @section DESCRIPTION
 * Two-dimensional dam break problem.
 **/
#ifndef TSUNAMI_LAB_SETUPS_DAM_BREAK_2D_H
#define TSUNAMI_LAB_SETUPS_DAM_BREAK_2D_H

#include "Setup.h"

namespace tsunami_lab {
  namespace setups {
    class DamBreak2d;
  }
}

/**
 * 1d dam break setup.
 **/
class tsunami_lab::setups::DamBreak2d: public Setup {
  private:
    
    //! height on the inside 
    t_real m_heightInner = 0;
    
    //! height on the outside
    t_real m_heightOuter = 0;

    //! x coordinate of the center of the dam
    t_real m_locationX = 0;

    //! x coordinate of the center of the dam
    t_real m_locationY = 0;

    //! radius of the dam
    t_real m_radius = 0;

  public:
    /**
     * Constructor.
     *
     * @param i_heightInner water height on the inside of the dam.
     * @param i_heightOuter water height on the outside of the dam.
     * @param i_locationX location (x-coordinate) of the dam center.
     * @param i_locationY location (y-coordinate) of the dam center.
     * @param i_radius radius of the dam.
     **/
    DamBreak2d( t_real i_heightInner,
                t_real i_heightOuter,
                t_real i_locationX,
                t_real i_locationY,
                t_real i_radius );

    /**
     * Gets the water height at a given point.
     *
     * @param i_x x-coordinate of the queried point.
     * @param i_y y-coordinate of the queried point.
     * @return height at the given point.
     **/
    t_real getHeight( t_real i_x, t_real i_y ) const;

    /**
     * Gets the water depth at a given point.
     * Positive values mean above sea level, negative values mean below sea level.
     *
     * @return water depth at the given point.
     **/
    t_real getBathymetry( t_real, t_real ) const;

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