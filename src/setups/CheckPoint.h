/**
 * @author Antonio Noack
 * @section DESCRIPTION
 * Checkpoint setup for pausing computations.
 **/
#ifndef TSUNAMI_LAB_SETUPS_CHECKPOINT_H
#define TSUNAMI_LAB_SETUPS_CHECKPOINT_H

#include "Setup.h"

// assert
#include <cassert>

namespace tsunami_lab {
  namespace setups {
    class CheckPoint;
  }
}

/**
 * checkpoint setup.
 **/
class tsunami_lab::setups::CheckPoint: public Setup {
  private:
    
    t_idx m_sizeX, m_sizeY, m_stride;
    
    t_real* m_height;
    t_real* m_bathymetry;
    t_real* m_momentumX;
    t_real* m_momentumY;
    
    /**
     * Get the index in the 2d grid; Clamped for Ghost-Cells.
     *
     * @param i_x x-coordinate of the queried point.
     * @param i_y y-coordinate of the queried point.
     * @return index for the queried point.
     **/
    t_idx getIndex( t_real i_x, t_real i_y ) const;
    
  public:
    /**
     * Constructor.
     *
     * todo
     **/
    explicit CheckPoint( t_real* i_height,
                         t_real* i_bathymetry,
                         t_real* i_momentumX,
                         t_real* i_momentumY,
                         t_idx   i_sizeX,
                         t_idx   i_sizeY,
                         t_idx   i_stride ) : 
      m_sizeX(i_sizeX), m_sizeY(i_sizeY), m_stride(i_stride),
      m_height(i_height),
      m_bathymetry(i_bathymetry),
      m_momentumX(i_momentumX),
      m_momentumY(i_momentumY) {
      
      assert(i_sizeX > 0 && i_sizeY > 0 && i_stride > 0);
      assert(i_height != nullptr);
      assert(i_bathymetry != nullptr);
      assert(i_momentumX != nullptr);
      assert(i_momentumY != nullptr);
      
    }
    
    ~CheckPoint() {
      delete[] m_height;
      delete[] m_bathymetry;
      delete[] m_momentumX;
      delete[] m_momentumY;
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
     * Gets the momentum in x-direction.
     *
     * @param i_x x-coordinate of the queried point.
     * @param i_y y-coordinate of the queried point.
     * @return momentum in x-direction.
     **/
    t_real getMomentumX( t_real i_x, t_real i_y ) const;

    /**
     * Gets the momentum in y-direction.
     *
     * @param i_x x-coordinate of the queried point.
     * @param i_y y-coordinate of the queried point.
     * @return momentum in y-direction.
     **/
    t_real getMomentumY( t_real i_x, t_real i_y ) const;

    /**
     * Gets the earth quake displacement.
     * Positive values mean the ground moved upwards.
     * Not required for this setup.
     *
     * @return displacement in meters at the given point.
     **/
    t_real getDisplacement( t_real, t_real ) const;
    
};

#endif // TSUNAMI_LAB_SETUPS_CHECKPOINT_H