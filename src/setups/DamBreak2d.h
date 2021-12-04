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
    struct ObstacleRect {
      //! the obstacle has the bounds [x0,x1[, [y0,y1[ in the grid.
      t_real x0, x1, y0, y1;
      t_real bathymetryOverride;
    };
    class DamBreak2d;
  }
}

/**
 * 2d dam break setup.
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

    //! bathymetry of the whole region except for the obstacle; caution: values > 0 make the water stop flowing!
    t_real m_bathymetry = 0;
    
    //! obstacle, by default empty
    ObstacleRect m_obstacle = { 0, 0, 0, 0, 0 };

  public:
    /**
     * Constructor.
     *
     * @param i_heightInner water height on the inside of the dam.
     * @param i_heightOuter water height on the outside of the dam.
     * @param i_locationX location (x-coordinate) of the dam center.
     * @param i_locationY location (y-coordinate) of the dam center.
     * @param i_radius radius of the dam.
     * @param i_bathymetry bathymetry in the whole region except for the obstacle.
     **/
    DamBreak2d( t_real i_heightInner,
                t_real i_heightOuter,
                t_real i_locationX,
                t_real i_locationY,
                t_real i_radius,
				t_real i_bathymetry ): m_heightInner(i_heightInner), m_heightOuter(i_heightOuter),
				m_locationX(i_locationX), m_locationY(i_locationY), m_radius(i_radius), m_bathymetry(i_bathymetry) {}
    
    /**
     * Defines the (optional) obstacle.
     *
     * @param i_x0 start (inclusive) of obstacle in x-dimension.
     * @param i_x1 end (exclusive) of obstacle in x-dimension.
     * @param i_y0 start (inclusive) of obstacle in y-dimension.
     * @param i_y1 end (exclusive) of obstacle in y-dimension.
     * @param i_bathymetryOverride bathymetry in this rectangle.
     **/
    void setObstacle( t_real i_x0,
                      t_real i_x1,
                      t_real i_y0,
                      t_real i_y1,
                      t_real i_bathymetryOverride ) {
      m_obstacle.x0 = i_x0;
      m_obstacle.x1 = i_x1;
      m_obstacle.y0 = i_y0;
      m_obstacle.y1 = i_y1;
      m_obstacle.bathymetryOverride = i_bathymetryOverride;
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
     * @return momentum in x-direction.
     **/
    t_real getMomentumX( t_real, t_real ) const;

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

#endif