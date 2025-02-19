/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 * 
 * @section LICENSE
 * Copyright 2020, Friedrich Schiller University Jena
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Simulation setup.
 **/
#ifndef TSUNAMI_LAB_SETUPS_SETUP_H
#define TSUNAMI_LAB_SETUPS_SETUP_H

#include "../constants.h"

namespace tsunami_lab {
  namespace setups {
    class Setup;
  }
}

/**
 * Base setup.
 **/
class tsunami_lab::setups::Setup {
  private:
    //! scale for setup; set on init of WavePropagation
    t_real m_initScaleX = 1, m_initScaleY = 1;
  
  public:
    /**
     * Virtual destructor for base class.
     **/
    virtual ~Setup(){};

    /**
     * Gets the water height at a given point.
     *
     * @param i_x x-coordinate of the queried point.
     * @param i_y y-coordinate of the queried point.
     * @return water height at the given point.
     **/
    virtual t_real getHeight( t_real i_x,
                              t_real i_y ) const = 0;

    /**
     * Gets the water depth at a given point.
     * Positive values mean above sea level, negative values mean below sea level.
     *
     * @param i_x x-coordinate of the queried point.
     * @param i_y y-coordinate of the queried point.
     * @return water depth at the given point.
     **/
    virtual t_real getBathymetry( t_real i_x,
                                  t_real i_y ) const = 0;

    /**
     * Gets the earth quake displacement.
     * Positive values mean the ground moved upwards.
     *
     * @param i_x x-coordinate of the queried point.
     * @param i_y y-coordinate of the queried point.
     * @return displacement in meters at the given point.
     **/
    virtual t_real getDisplacement( t_real i_x, t_real i_y ) const = 0;

    /**
     * Gets the momentum in x-direction.
     *
     * @param i_x x-coordinate of the queried point.
     * @param i_y y-coordinate of the queried point.
     * @return momentum in x-direction.
     **/
    virtual t_real getMomentumX( t_real i_x,
                                 t_real i_y ) const = 0;

    /**
     * Gets the momentum in y-direction.
     *
     * @param i_x x-coordinate of the queried point.
     * @param i_y y-coordinate of the queried point.
     * @return momentum in y-direction.
     **/
    virtual t_real getMomentumY( t_real i_x,
                                 t_real i_y ) const = 0;
    
    /**
     * Gets the init scales for x and y.
     *
     * @param o_initScaleX scale on x axis.
     * @param o_initScaleY scale on y axis.
     **/
    void getInitScale(t_real &o_initScaleX, t_real &o_initScaleY){
      o_initScaleX = m_initScaleX;
      o_initScaleY = m_initScaleY;
    }
    
    /**
     * Sets the init scales for x and y.
     *
     * @param i_initScaleX scale on x axis.
     * @param i_initScaleY scale on y axis.
     **/
    void setInitScale(t_real i_initScaleX, t_real i_initScaleY){
      m_initScaleX = i_initScaleX;
      m_initScaleY = i_initScaleY;
    }
};

#endif