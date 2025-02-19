/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de), Antonio Noack
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
 * One-dimensional wave propagation patch.
 **/
#ifndef TSUNAMI_LAB_PATCHES_WAVE_PROPAGATION_1D
#define TSUNAMI_LAB_PATCHES_WAVE_PROPAGATION_1D

#include "WavePropagation.h"
#include "../setups/Setup.h"

namespace tsunami_lab {
  namespace patches {
    class WavePropagation1d;
  }
}

class tsunami_lab::patches::WavePropagation1d: public WavePropagation {
  private:
    //! current step which indicates the active values in the arrays below
    unsigned short m_step = 0;
    
    //! number of cells discretizing the computational domain
    t_idx m_nCells = 0;
    
    //! water heights for the current and next time step for all cells
    t_real * m_h[2] = { nullptr, nullptr };
    
    //! momenta for the current and next time step for all cells
    t_real * m_hu[2] = { nullptr, nullptr };
    
    //! bathymetry in meters for all cells
    t_real * m_bathymetry = nullptr;
    
    //! if true, use FWave, else use Roe solver
    bool m_useFWaveSolver = true;
    
  public:
    /**
     * Constructs the 1d wave propagation solver.
     *
     * @param i_nCells number of cells.
     **/
    WavePropagation1d( t_idx i_nCells );
    
    /**
     * Constructs the 1d wave propagation solver and applies the setup.
     *
     * @param i_nCells number of cells.
     * @param i_setup setup for cell initialization.
     * @param i_scale scale for the scene; e.g. you can multiply the number of cells by x, and set the scale to 1/x, and your setup will still work.
     **/
    WavePropagation1d( t_idx i_nCells, tsunami_lab::setups::Setup* i_setup, t_real i_scale );

    /**
     * Destructor which frees all allocated memory.
     **/
    ~WavePropagation1d();
    
    /**
     * Initializes the internal state with a setup.
     *
     * @param i_setup setup for cell initialization.
     * @param i_scale scale for the scene; e.g. you can multiply the number of cells by x, and set the scale to 1/x, and your setup will still work.
     **/
    void initWithSetup( tsunami_lab::setups::Setup* i_setup, t_real i_scale );
    
    /**
     * Computes the maximum time step that is allowed without breaking the CFL condition.
     **/
    t_real computeMaxTimestep( t_real i_cellSizeMeters );
    
    /**
     * Computes the maximum time step that is allowed without breaking the CFL condition within a certain radius from the center.
     **/
    t_real computeMaxTimestep( t_real i_cellSizeMeters, t_idx i_updateRadius );
    
    /**
     * Performs a time step.
     *
     * @param i_scaling scaling of the time step (dt / dx).
     **/
    void timeStep( t_real i_scaling );
    
    /**
     * Performs a time step on the inner region.
     *
     * @param i_scaling scaling of the time step (dt / dx).
     * @param i_size radius of what is actually updated.
     **/
    void timeStep( t_real i_scaling, t_idx i_updateRadius );
    
    /**
     * Sets the values of the ghost cells according to outflow boundary conditions.
     **/
    void setGhostOutflow();
    
    /**
     * Gets the stride in y-direction. x-direction is stride-1.
     *
     * @return stride in y-direction.
     **/
    t_idx getStride(){
      return m_nCells+2;
    }
    
    /**
     * Gets cells' water heights.
     *
     * @return water heights.
     */
    t_real const * getHeight(){
      return m_h[m_step]+1;
    }
    
    /**
     * Gets the cells' momenta in x-direction.
     *
     * @return momenta in x-direction.
     **/
    t_real const * getMomentumX(){
      return m_hu[m_step]+1;
    }
    
    /**
     * Dummy function which returns a nullptr.
     **/
    t_real const * getMomentumY(){
      return nullptr;
    }
    
    /**
     * Gets the cells' bathymetry.
     *
     * @return bathymetry.
     **/
    t_real const * getBathymetry(){
      return m_bathymetry+1;
    }
    
    /**
     * Sets the bathymetry of the cell to the given value.
     *
     * @param i_ix id of the cell in x-direction.
     * @param i_b bathymetry.
     **/
    void setBathymetry( t_idx  i_ix,
                        t_idx,
                        t_real i_b ) {
      m_bathymetry[i_ix+1] = i_b;
    }
    
    /**
     * Sets the height of the cell to the given value.
     *
     * @param i_ix id of the cell in x-direction.
     * @param i_h water height.
     **/
    void setHeight( t_idx  i_ix,
                    t_idx,
                    t_real i_h ) {
      m_h[m_step][i_ix+1] = i_h;
    }
    
    /**
     * Sets the momentum in x-direction to the given value.
     *
     * @param i_ix id of the cell in x-direction.
     * @param i_hu momentum in x-direction.
     **/
    void setMomentumX( t_idx  i_ix,
                       t_idx,
                       t_real i_hu ) {
      m_hu[m_step][i_ix+1] = i_hu;
    }
    
    /**
     * Dummy function since there is no y-momentum in the 1d solver.
     **/
    void setMomentumY( t_idx,
                       t_idx,
                       t_real ) {};
};

#endif