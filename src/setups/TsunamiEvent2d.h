/**
 * @author Antonio Noack
 * @section DESCRIPTION
 * Two-dimensional tsunami event problem.
 **/
#ifndef TSUNAMI_LAB_SETUPS_TSUNAMI_EVENT_2D_H
#define TSUNAMI_LAB_SETUPS_TSUNAMI_EVENT_2D_H

#include "Setup.h"

// assert
#include <cassert>

namespace tsunami_lab {
  namespace setups {
    class TsunamiEvent2d;
  }
}

/**
 * 2d dam break setup.
 **/
class tsunami_lab::setups::TsunamiEvent2d: public Setup {
  private:
    
    t_real m_offset = 0, m_shoreCliffHeight = 20;
    t_idx  m_scaleBath = 1, m_sizeXBath, m_sizeYBath, m_strideBath;
    t_idx  m_scaleDisp = 1, m_sizeXDisp, m_sizeYDisp, m_strideDisp;
     
    t_real* m_bathymetry;
    t_real* m_displacement;
    
    /**
     * Interpolate a value in a 2d grid.
     *
     * @param i_x x-coordinate of the queried point times scale.
     * @param i_y y-coordinate of the queried point times scale.
     * @param i_stride stride of the data, typically size.x.
	 * @param i_sizeX data size in x dimension.
	 * @parma i_sizeY data size in y direction.
     * @param i_data float data with size of m_sizeX x m_sizeY and stride m_stride.
     * @return linearly interpolated value at the queried point.
     **/
    t_real getInterpolatedValue( t_real i_x, t_real i_y, t_idx i_sizeX, t_idx i_sizeY, t_idx i_stride, t_real* i_data ) const;
    
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
    explicit TsunamiEvent2d( t_real* i_bathymetry,
                             t_idx   i_sizeXBath,
                             t_idx   i_sizeYBath,
                             t_idx   i_strideBath,
                             t_real  i_scaleBath,
                             t_real* i_displacement,
                             t_idx   i_sizeXDisp,
                             t_idx   i_sizeYDisp,
                             t_idx   i_strideDisp,
							 t_idx   i_scaleDisp) : 
      m_offset(0.5),
      m_scaleBath(i_scaleBath), m_sizeXBath(i_sizeXBath), m_sizeYBath(i_sizeYBath), m_strideBath(i_strideBath),
      m_scaleDisp(i_scaleDisp), m_sizeXDisp(i_sizeXDisp), m_sizeYDisp(i_sizeYDisp), m_strideDisp(i_strideDisp),
      m_bathymetry(i_bathymetry),
      m_displacement(i_displacement) {
      
      assert(i_sizeXBath > 0 && i_sizeYBath > 0);
      assert(i_sizeXDisp > 0 && i_sizeYDisp > 0);
      assert(i_bathymetry != nullptr);
      assert(i_displacement != nullptr);
      
    }
    
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
     * @param i_x x-coordinate of the queried point.
     * @param i_y y-coordinate of the queried point.
     * @return water depth at the given point.
     **/
    t_real getBathymetry( t_real i_x, t_real i_y ) const;

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

#endif // TSUNAMI_LAB_SETUPS_TSUNAMI_EVENT_2D_H