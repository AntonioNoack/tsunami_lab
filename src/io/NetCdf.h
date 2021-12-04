/**
 * @author Antonio Noack
 * 
 * @section DESCRIPTION
 * IO-routines for writing a snapshot/checkpoint as a NetCDF file.
 **/
#ifndef TSUNAMI_LAB_IO_NETCDF_H
#define TSUNAMI_LAB_IO_NETCDF_H

#include <string>
#include <vector>

#include "../constants.h"
#include "../io/Station.h"
#include "../patches/WavePropagation.h"
#include "../setups/Setup.h"
#include "../setups/CheckPoint.h"

namespace tsunami_lab {
  namespace io {
    class NetCDF;
  }
}

class tsunami_lab::io::NetCDF {
  private:
    static void downsample( t_idx         i_sizeInX,
                            t_idx         i_sizeInY,
                            t_idx         i_strideIn,
                            t_real const* i_dataIn,
                            t_idx         i_sizeOutX,
                            t_idx         i_sizeOutY,
                            t_idx         i_strideOut,
                            t_real*       i_dataOut,
                            t_idx         i_step );
  public:
    /**
     * Writes the data as a NetCDF file.
     *
     * @param i_cellSizeMeters cell size in x- and y-direction.
     * @param i_nx number of cells in x-direction.
     * @param i_ny number of cells in y-direction.
     * @param i_gridOffsetX x-coordinate in meters of first cell, excluding the ghost cells.
     * @param i_gridOffsetY y-coordinate in meters of first cell, excluding the ghost cells.
     * @param i_step only every step-th cell is written.
     * @param i_stride stride of the data arrays in y-direction (x is assumed to be stride-1).
     * @param i_h water height of the cells.
     * @param i_hu momentum in x-direction of the cells; optional: use nullptr if not required.
     * @param i_hv momentum in y-direction of the cells; optional: use nullptr if not required.
     * @param i_b bathymetry data of the cells.
     * @param i_setup setup for displacement data.
     * @param i_time time in seconds of the current frame.
     * @param i_frameIndex n-1 for n-th frame, 0 for first frame.
     * @param i_deflateLevel compression level, 0 = large/fastest, 9 = compact/slowest, 5 is recommended
     * @param i_fileName target file, where the data is written to.
     * @return 0 if the function was successful, -1 or error code else.
     **/
    static int appendTimeframe( t_real                       i_cellSizeMeters,
                                t_idx                        i_nx,
                                t_idx                        i_ny,
                                t_real                       i_gridOffsetX,
                                t_real                       i_gridOffsetY,
                                t_idx                        i_step,
                                t_idx                        i_stride,
                                t_real               const * i_h,
                                t_real               const * i_hu,
                                t_real               const * i_hv,
                                t_real               const * i_b,
                                tsunami_lab::setups::Setup * i_setup,
                                t_real                       i_time,
                                t_idx                        i_frameIndex,
                                int                          i_deflateLevel,
                                std::string                  i_fileName );
    
    /**
     * Reads a checkpint from a file.
     *
     * @param i_fileName input file name.
     * @param o_nx field size on x axis.
     * @param o_ny field size on y axis.
     * @param o_cellSizeMeters cell size in meters.
     * @param o_cflFactor cfl factor for 2d simulations, 0.5 else.
     * @param o_simulationTime last computed simulation time in seconds.
     * @param o_timeStepIndex last time step index.
     * @param o_stations stations with previously collected values.
     * @return a setup with height, bathymetry and impulse data.
     **/
    static tsunami_lab::setups::Setup* loadCheckpoint( std::string                             i_fileName,
                                                       t_idx                                 & o_nx,
                                                       t_idx                                 & o_ny,
                                                       t_real                                & o_cellSizeMeters,
                                                       t_real                                & o_cflFactor,
                                                       double                                & o_simulationTime,
                                                       t_idx                                 & o_timeStepIndex,
                                                       std::vector<tsunami_lab::io::Station> & o_stations );
    
    /**
     * Writes the current simulation data as a NetCDF file / checkpoint.
     *
     * @param i_fileName output file name.
     * @param i_nx field size x axis.
     * @param i_ny field size y axis.
     * @param i_cellSizeMeters cell size in meters.
     * @param i_cflFactor cfl factor for 2d simulations.
     * @param i_simulationTime current simulation time in seconds.
     * @param i_timeStepIndex n for the n-th time step.
     * @param i_stations stations.
     * @param i_waveProp wave propagation instance.
     * @return 0 if successful, error code else
     **/
    static int storeCheckpoint( std::string                             i_fileName,
                                t_idx                                   i_nx,
                                t_idx                                   i_ny,
                                t_real                                  i_cellSizeMeters,
                                t_real                                  i_cflFactor,
                                double                                  i_simulationTime,
                                t_idx                                   i_timeStepIndex,
                                std::vector<tsunami_lab::io::Station> & i_stations,
                                tsunami_lab::patches::WavePropagation * i_waveProp );
    
    /**
     * Reads a 2d array of data from a NetCDF file.
     *
     * @param i_fileName file name for the file to be loaded.
     * @param i_variableName name of the variable to be loaded; default value for GMT data: z.
     * @param o_sizeX size of the loaded data in x-direction.
     * @param o_sizeY size of the loaded data in y-direction.
     * @param o_cellSizeMeters size of a single cell in meters; square cells are assumed.
     * @param o_gridOffsetX x-coordinate in meters of first cell, excluding the ghost cells.
     * @param o_gridOffsetY y-coordinate in meters of first cell, excluding the ghost cells.
     * @param o_data data vector, where the data is written to.
     * @return 0 if the function was successful, -1 or error code else.
     **/
    static int load2dArray( std::string           i_fileName,
                            std::string           i_variableName,
                            t_idx               & o_sizeX,
                            t_idx               & o_sizeY,
                            t_real              & o_cellSizeMeters,
                            t_real              & o_gridOffsetX,
                            t_real              & o_gridOffsetY,
                            std::vector<t_real> & o_data );
};

#endif // TSUNAMI_LAB_IO_NETCDF_H