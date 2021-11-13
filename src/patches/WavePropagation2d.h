/**
 * @author Antonio Noack
 * 
 * @section DESCRIPTION
 * Two-dimensional wave propagation patch.
 **/
#ifndef TSUNAMI_LAB_PATCHES_WAVE_PROPAGATION_2D
#define TSUNAMI_LAB_PATCHES_WAVE_PROPAGATION_2D

#include "WavePropagation.h"
#include "../setups/Setup.h"

namespace tsunami_lab {
  namespace patches {
    class WavePropagation2d;
  }
}

class tsunami_lab::patches::WavePropagation2d: public WavePropagation {
  private:
    
	//! current index for momentum buffers
	t_idx m_step = 0;
	
    //! number of cells discretizing the computational domain
    t_idx m_nCells = 0;
	
	//! number of cells on the x and y axis
	t_idx m_nCellsX = 0, m_nCellsY;
    
    //! water heights for the current and next time step for all cells
	//! updated twice per step
    t_real * m_h[2] = { nullptr, nullptr };
    
    //! momenta in x direction for the current and next time step for all cells
	//! updated once per step
    t_real * m_hu[2] = { nullptr, nullptr };
    
    //! momenta in y direction for the current and next time step for all cells
	//! updated once per step
    t_real * m_hv[2] = { nullptr, nullptr };
    
    //! bathymetry in meters for all cells
	//! currently only updated at the start of simulation
    t_real * m_bathymetry = nullptr;
    
    //! if true, use FWave, else use Roe solver
    bool m_useFWaveSolver = true;
    
  public:
    /**
     * Constructs the 2d wave propagation solver.
     *
     * @param i_nCellsX number of cells on the x axis.
	 * @param i_nCellsY number of cells on the y axis.
     **/
    WavePropagation2d( t_idx i_nCellsX, t_idx i_nCellsY );
    
    /**
     * Constructs the 1d wave propagation solver and applies the setup.
     *
     * @param i_nCellsX number of cells on the x axis.
	 * @param i_nCellsY number of cells on the y axis.
     * @param i_setup setup for cell initialization.
     * @param i_scaleX scale for the scene in x direction; e.g. you can multiply the number of cells by x, and set the scale to 1/x, and your setup will still work.
     * @param i_scaleY scale for the scene in y direction
     **/
    WavePropagation2d( t_idx i_nCellsX, t_idx i_nCellsY, tsunami_lab::setups::Setup* i_setup, t_real i_scaleX, t_real i_scaleY );

    /**
     * Destructor which frees all allocated memory.
     **/
    ~WavePropagation2d();
    
    /**
     * Initializes the internal state with a setup.
     *
     * @param i_setup setup for cell initialization.
     * @param i_scaleX scale for the scene in x direction; e.g. you can multiply the number of cells by x, and set the scale to 1/x, and your setup will still work.
     * @param i_scaleY scale for the scene in y direction
     **/
    void initWithSetup( tsunami_lab::setups::Setup* i_setup, t_real i_scaleX, t_real i_scaleY );
    
    /**
     * Computes the maximum time step that is allowed without breaking the CFL condition.
     **/
    t_real computeMaxTimestep( t_real i_cellSizeMeters );
	
	/**
	 * Internal function, which computes and applies the update from two neighbor cells to each other
	 
	 todo params
     **/
    void internalUpdate( t_real i_scaling, t_idx l_ceL, t_idx l_ceR, 
        t_real* l_hOld, t_real* l_huOld, t_real* l_hNew, t_real* l_huNew );
    
    /**
     * Performs a time step.
     *
     * @param i_scaling scaling of the time step (dt / dx).
     **/
    void timeStep( t_real i_scaling );
    
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
      return m_nCellsX+2;
    }
    
    /**
     * Gets cells' water heights.
     *
     * @return water heights.
     */
    t_real const * getHeight(){
      return m_h[0]+1;
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
     * Gets the cells' momenta in y-direction.
     *
     * @return momenta in y-direction.
     **/
    t_real const * getMomentumY(){
      return m_hv[m_step]+1;
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
     * @param i_iy id of the cell in y-direction.
     * @param i_b bathymetry.
     **/
    void setBathymetry( t_idx  i_ix,
                        t_idx  i_iy,
                        t_real i_b ) {
      m_bathymetry[(i_ix+1) + (i_iy+1) * (m_nCellsX+2)] = i_b;
    }
    
    /**
     * Sets the height of the cell to the given value.
     *
     * @param i_ix id of the cell in x-direction.
     * @param i_iy id of the cell in y-direction.
     * @param i_h water height.
     **/
    void setHeight( t_idx  i_ix,
                    t_idx  i_iy,
                    t_real i_h ) {
      m_h[0][(i_ix+1) + (i_iy+1) * (m_nCellsX+2)] = i_h;
    }
    
    /**
     * Sets the momentum in x-direction to the given value.
     *
     * @param i_ix id of the cell in x-direction.
     * @param i_iy id of the cell in y-direction.
     * @param i_hu momentum in x-direction.
     **/
    void setMomentumX( t_idx  i_ix,
                       t_idx  i_iy,
                       t_real i_hu ) {
      m_hu[m_step][(i_ix+1) + (i_iy+1) * (m_nCellsX+2)] = i_hu;
    }
    
    /**
     * Sets the momentum in y-direction to the given value.
     *
     * @param i_ix id of the cell in x-direction.
     * @param i_iy id of the cell in y-direction.
     * @param i_hv momentum in y-direction.
     **/
    void setMomentumY( t_idx  i_ix,
                       t_idx  i_iy,
                       t_real i_hv) {
      m_hv[m_step][(i_ix+1) + (i_iy+1) * (m_nCellsX+2)] = i_hv;
	};
};

#endif