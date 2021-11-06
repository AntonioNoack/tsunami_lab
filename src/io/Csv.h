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
 * IO-routines for writing a snapshot as Comma Separated Values (CSV).
 **/
#ifndef TSUNAMI_LAB_IO_CSV
#define TSUNAMI_LAB_IO_CSV

#include "../constants.h"
#include <cstring>
#include <iostream>
#include <vector>
#include <optional>

namespace tsunami_lab {
  namespace io {
    class Csv;
  }
}

class tsunami_lab::io::Csv {
  public:
    /**
     * Writes the data as CSV to the given stream.
     *
     * @param i_dxy cell width in x- and y-direction.
     * @param i_nx number of cells in x-direction.
     * @param i_ny number of cells in y-direction.
     * @param i_step only every step-th cell is written.
     * @param i_stride stride of the data arrays in y-direction (x is assumed to be stride-1).
     * @param i_h water height of the cells; optional: use nullptr if not required.
     * @param i_hu momentum in x-direction of the cells; optional: use nullptr if not required.
     * @param i_hv momentum in y-direction of the cells; optional: use nullptr if not required.
     * @param io_stream stream to which the CSV-data is written.
     **/
    static void write( t_real               i_dxy,
                       t_idx                i_nx,
                       t_idx                i_ny,
                       t_idx                i_step,
                       t_idx                i_stride,
                       t_real       const * i_h,
                       t_real       const * i_hu,
                       t_real       const * i_hv,
                       std::ostream       & io_stream );
    
    /**
     * Writes the data as CSV to the given stream.
     *
     * @param i_dxy cell width in x- and y-direction.
     * @param i_nx number of cells in x-direction.
     * @param i_ny number of cells in y-direction.
     * @param i_stride stride of the data arrays in y-direction (x is assumed to be stride-1).
     * @param i_h water height of the cells.
     * @param i_hu momentum in x-direction of the cells; optional: use nullptr if not required.
     * @param i_hv momentum in y-direction of the cells; optional: use nullptr if not required.
     * @param i_b bathymetry data of the cells.
     * @param io_stream stream to which the CSV-data is written.
     **/
    static void write( t_real                      i_dxy,
                       t_idx                       i_nx,
                       t_idx                       i_ny,
                       t_idx                       i_step,
                       t_idx                       i_stride,
                       t_real              const * i_h,
                       t_real              const * i_hu,
                       t_real              const * i_hv,
                       t_real              const * i_b,
                       std::ostream              & io_stream );
    
    /**
     * Reads numeric data in CSV format from a given stream.
     *
     * @param io_stream stream from which data is read.
     * @return vector of attributes: first the name of the column, then all values.
     **/
    static std::vector<std::pair<std::string, std::vector<tsunami_lab::t_real>>> read( std::istream& io_stream );
    
    /**
     * Finds the column inside the csv data.
     *
     * @param i_csvData the csv to search in.
     * @param i_columnName the name of the column.
     * @return the data vector, if found; else an empty vector.
     **/
    static std::vector<tsunami_lab::t_real> findColumn(
      std::vector<std::pair<std::string, std::vector<tsunami_lab::t_real>>> &i_csvData, std::string i_columnName
    );
};

#endif