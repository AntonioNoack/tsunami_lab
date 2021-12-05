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
 * 2d tsunami setup.
 **/
class tsunami_lab::setups::TsunamiEvent2d: public Setup {
  private:
    
    t_real m_offset = 0, m_shoreCliffHeight = 20;
    t_real m_scaleBath = 1; t_idx m_sizeXBath, m_sizeYBath, m_strideBath;
    t_real m_scaleDisp = 1; t_idx m_sizeXDisp, m_sizeYDisp, m_strideDisp;
     
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
     * @param i_data float data with size of i_sizeX x i_sizeY and stride i_stride.
     * @return linearly interpolated value at the queried point.
     **/
    t_real getInterpolatedValue( t_real i_x, t_real i_y, t_idx i_sizeX, t_idx i_sizeY, t_idx i_stride, t_real* i_data ) const;
    
  public:
    /**
     * Constructor.
     *
     * @param i_bathymetry bathymetry data, a vector of length i_size, must not be null.
     * @param i_sizeXBath, i_sizeYBath size of the bathymetry data.
	 * @param i_strideBath bathymetry data stride.
	 * @param i_scaleBath bathymetry data scale.
     * @param i_displacement displacement data, a vector of length i_size, must not be null.
     * @param i_sizeXDisp, i_sizeYDisp size of the displacement data.
	 * @param i_strideDisp displacement data stride.
	 * @param i_scaleDisp displacement data scale.
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
							 t_real  i_scaleDisp) : 
      m_offset(0.5),
      m_scaleBath(i_scaleBath), m_sizeXBath(i_sizeXBath), m_sizeYBath(i_sizeYBath), m_strideBath(i_strideBath),
      m_scaleDisp(i_scaleDisp), m_sizeXDisp(i_sizeXDisp), m_sizeYDisp(i_sizeYDisp), m_strideDisp(i_strideDisp),
      m_bathymetry(i_bathymetry),
      m_displacement(i_displacement) {
      
      assert(i_sizeXBath > 1 && i_sizeYBath > 1 && i_scaleBath > 0);
      assert(i_sizeXDisp > 1 && i_sizeYDisp > 1 && i_scaleDisp > 0);
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