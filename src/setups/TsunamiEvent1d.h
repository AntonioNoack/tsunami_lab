/**
 * @author Antonio Noack
 * @section DESCRIPTION
 * One-dimensional single discontinuity problem.
 **/
#ifndef TSUNAMI_LAB_SETUPS_DISCONTINUITY_1D_H
#define TSUNAMI_LAB_SETUPS_DISCONTINUITY_1D_H

#include "Setup.h"

// assert
#include <cassert>

namespace tsunami_lab {
  namespace setups {
    class TsunamiEvent1d;
  }
}

/**
 * 1d dam break setup.
 **/
class tsunami_lab::setups::TsunamiEvent1d: public Setup {
  private:
    
    t_real m_scale = 1, m_offset = 0, m_shoreCliffHeight = 20;
    t_idx  m_size;
    t_real m_displacementStart, m_displacementEnd;
    t_real m_displacement = 0;
    
    t_real* m_bathymetry;
    
    /**
     * Get the internal, interpolated bathymetry at that point.
     *
     * @param i_x x-coordinate of the queried point.
     * @return raw bathymetry without coast-fix.
     **/
    t_real getRawBathymetry( t_real i_x ) const;
  
  public:
    /**
     * Constructor.
     *
     * @param i_bathymetry bathymetry data, a vector of length i_size, must not be null.
     * @param i_size size of the bathymetry data array, must be > 0.
     * @parma i_scale how much the scene is scaled up, typically <= 1.
     * @param i_displacementStart index of the start of the displacement sine wave.
     * @param i_displacementEnd index of the end of the displacement sine wave.
     * @param i_displacement height of displacement.
     **/
    explicit TsunamiEvent1d( t_real* i_bathymetry,
                             t_idx   i_size,
                             t_real  i_scale,
                             t_real  i_displacementStart,
                             t_real  i_displacementEnd,
                             t_real  i_displacement ) : 
      m_scale(i_scale), m_offset(0), m_size(i_size),
      m_displacementStart(i_displacementStart),
      m_displacementEnd(i_displacementEnd),
      m_displacement(i_displacement),
      m_bathymetry(i_bathymetry) {
      
      assert(i_size > 0);
      assert(i_bathymetry != nullptr);
      
    }
    
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
     * Gets the earth quake displacement.
     * Positive values mean the ground moved upwards.
     *
     * @param i_x x-coordinate of the queried point.
     * @return displacement in meters at the given point.
     **/
    t_real getDisplacement( t_real i_x, t_real ) const;
    
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

#endif // TSUNAMI_LAB_SETUPS_DISCONTINUITY_1D_H