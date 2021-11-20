/**
 * @author Antonio Noack
 * 
 * @section DESCRIPTION
 * IO-routines for writing a snapshot as a NetCDF file.
 **/
#ifndef TSUNAMI_LAB_IO_NETCDF_H
#define TSUNAMI_LAB_IO_NETCDF_H

#include <cstring>
#include <iostream>
#include <vector>
#include <optional>

#include "../constants.h"
#include "../setups/Setup.h"

namespace tsunami_lab {
  namespace io {
    class NetCDF;
  }
}

class tsunami_lab::io::NetCDF {
  public:
    /**
     * Writes the data as NetCDF to the given stream.
     *
     * @param i_cellSizeMeters cell size in x- and y-direction.
     * @param i_nx number of cells in x-direction.
     * @param i_ny number of cells in y-direction.
     * @param i_step only every step-th cell is written.
     * @param i_stride stride of the data arrays in y-direction (x is assumed to be stride-1).
     * @param i_h water height of the cells.
     * @param i_hu momentum in x-direction of the cells; optional: use nullptr if not required.
     * @param i_hv momentum in y-direction of the cells; optional: use nullptr if not required.
     * @param i_b bathymetry data of the cells.
     * @param i_setup setup for displacement data.
     * @param i_time time in seconds of the current frame.
     * @param i_frameIndex n-1 for n-th frame, 0 for first frame.
     * @param i_fileName target file, where the data is written to.
	 * @return 0 if the function was successful, -1 or error code else.
     **/
    static int appendTimeframe( t_real                       i_dxy,
                                t_idx                        i_nx,
                                t_idx                        i_ny,
                                t_idx                        i_step,
                                t_idx                        i_stride,
                                t_real               const * i_h,
                                t_real               const * i_hu,
                                t_real               const * i_hv,
                                t_real               const * i_b,
								tsunami_lab::setups::Setup * i_setup,
								t_real                       i_time,
								t_idx                        i_frameIndex,
								std::string                  i_fileName );
    /**
     * Writes the data as NetCDF to the given stream.
     *
     * @param i_fileName file name for the file to be loaded.
     * @param i_variableName name of the variable to be loaded; default value for GMT data: z.
     * @param o_sizeX size of the loaded data in x-direction.
     * @param o_sizeY size of the loaded data in y-direction.
     * @param o_cellSizeMeters size of a single cell in meters; square cells are assumed.
     * @param o_data data vector, where the data is written to.
	 * @return 0 if the function was successful, -1 or error code else.
     **/
    static int load2dArray( std::string           i_fileName,
	                        std::string           i_variableName,
							t_idx               & o_sizeX,
							t_idx               & o_sizeY,
							t_real              & o_cellSizeMeters,
							std::vector<t_real> & o_data );
};

#endif // TSUNAMI_LAB_IO_NETCDF_H