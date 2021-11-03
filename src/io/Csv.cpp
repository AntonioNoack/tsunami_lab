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
 * IO-routines for writing a snapshot as Comma Separated Values (CSV).
 **/
#include "Csv.h"
#include <sstream>

void tsunami_lab::io::Csv::write( t_real               i_dxy,
                                  t_idx                i_nx,
                                  t_idx                i_ny,
                                  t_idx                i_stride,
                                  t_real       const * i_h,
                                  t_real       const * i_hu,
                                  t_real       const * i_hv,
                                  std::ostream       & io_stream ) {
  // write the CSV header
  io_stream << "x,y";
  if( i_h  != nullptr ) io_stream << ",height";
  if( i_hu != nullptr ) io_stream << ",momentum_x";
  if( i_hv != nullptr ) io_stream << ",momentum_y";
  io_stream << "\n";

  // iterate over all cells
  for( t_idx l_iy = 0; l_iy < i_ny; l_iy++ ) {
    for( t_idx l_ix = 0; l_ix < i_nx; l_ix++ ) {
      // derive coordinates of cell center
      t_real l_posX = i_nx > 1 ? (l_ix + 0.5) * i_dxy : 0;
      t_real l_posY = i_ny > 1 ? (l_iy + 0.5) * i_dxy : 0;

      t_idx l_id = l_iy * i_stride + l_ix;

      // write data
      io_stream << l_posX << "," << l_posY;
      if( i_h  != nullptr ) io_stream << "," << i_h[l_id];
      if( i_hu != nullptr ) io_stream << "," << i_hu[l_id];
      if( i_hv != nullptr ) io_stream << "," << i_hv[l_id];
      io_stream << "\n";
    }
  }
  io_stream << std::flush;
}

std::vector<std::pair<std::string, std::vector<tsunami_lab::t_real>>> tsunami_lab::io::Csv::read( std::istream &io_stream ){
  
  std::vector<std::pair<std::string, std::vector<tsunami_lab::t_real>>> l_result;
  
  char l_comma;
  
  // first find the header and define all properties
  std::string l_line;
  while(getline(io_stream, l_line)){
    if(!l_line.empty() && l_line[0] != '#'){
      // read column names, separated by whitespace
      for(t_idx l_startIndex = 0, l_currentIndex = 0, l_size = l_line.size(); l_currentIndex <= l_size; l_currentIndex++){
        char l_currentChar = l_currentIndex < l_size ? l_line[l_currentIndex] : ' ';
        if(isspace(l_currentChar) || l_currentChar == ',' || l_currentChar == ';'){
          if(l_currentIndex > l_startIndex){
            // add new column to result
            std::string l_columnName = l_line.substr(l_startIndex, l_currentIndex-l_startIndex);
            l_result.push_back(std::pair<std::string, std::vector<t_real>>(l_columnName, std::vector<t_real>()));
          }
          l_startIndex = l_currentIndex + 1;
        }
      }
      break;
    }
  }
  
  t_idx l_numProperties = l_result.size();
  
  // read all data
  while(getline(io_stream, l_line)){
    if(!l_line.empty()){
      std::stringstream l_splitter(l_line);
      // split data by comma
      // not really ideal or secure, probably only works if all fields have values,
      // and there are no additional spaces
      for(t_idx l_i=0;l_i<l_numProperties;l_i++){
        t_real l_value = 0;
        l_splitter >> l_value;
        l_splitter >> l_comma;
        l_result[l_i].second.push_back(l_value);
      }
    }
  }
  
  return l_result;
  
}
