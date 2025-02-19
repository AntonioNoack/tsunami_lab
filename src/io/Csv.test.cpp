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
 * Unit tests for the CSV-interface.
 **/
#include <catch2/catch.hpp>
#include "../constants.h"
#include <sstream>
#define private public
#include "Csv.h"
#undef public

#define t_idx  tsunami_lab::t_idx
#define t_real tsunami_lab::t_real

TEST_CASE( "Test the CSV-writer for 1D settings.", "[CsvWrite1d]" ) {
  // define a simple example
  t_real l_h [7] = { 0, 1, 2, 3, 4, 5, 6 };
  t_real l_hu[7] = { 6, 5, 4, 3, 2, 1, 0 };

  std::stringstream l_stream0;
  tsunami_lab::io::Csv::write( 0.5,
                               5,
                               1,
                               1,
                               7,
                               l_h+1,
                               l_hu+1,
                               nullptr,
                               l_stream0 );

  std::string l_ref0 = R"V0G0N(x,y,height,momentum_x
0.25,0,1,5
0.75,0,2,4
1.25,0,3,3
1.75,0,4,2
2.25,0,5,1
)V0G0N";

  REQUIRE( l_stream0.str().size() == l_ref0.size() );
  REQUIRE( l_stream0.str() == l_ref0 );
}

TEST_CASE( "Test the CSV-writer for 2D settings.", "[CsvWrite2d]" ) {
  // define a simple example
  t_real l_h[16]  = {  0,  1,  2,  3,
                                    4,  5,  6,  7,
                                    8,  9, 10, 11,
                                   12, 13, 14, 15 };
  t_real l_hu[16] = { 15, 14, 13, 12,
                                   11, 10,  9,  8,
                                    7,  6,  5,  4,
                                    3,  2,  1,  0 };
  t_real l_hv[16] = {  0,  4,  8, 12,
                                    1,  5,  9, 13,
                                    2,  6, 10, 14,
                                    3,  7, 11, 15 };

  std::stringstream l_stream1;
  tsunami_lab::io::Csv::write( 10,
                               2,
                               2,
                               1,
                               4,
                               l_h+4+1,
                               l_hu+4+1,
                               l_hv+4+1,
                               l_stream1 );

  std::string l_ref1 = R"V0G0N(x,y,height,momentum_x,momentum_y
5,5,5,10,5
15,5,6,9,9
5,15,9,6,6
15,15,10,5,10
)V0G0N";

  REQUIRE( l_stream1.str().size() == l_ref1.size() );
  REQUIRE( l_stream1.str() == l_ref1 );
}


TEST_CASE( "Test the CSV Reader", "[CSVReader]" ) {
  
  // define a sample csv file
  std::string l_source = R"V0G0N(x,y,height,momentum_x
0.25,0,1,5
0.75,0,2,4
1.25,0,3,3
1.75,0,4,2
2.25,0,5,1
)V0G0N";
  std::stringstream l_inputStream(l_source);
  auto l_read = tsunami_lab::io::Csv::read(l_inputStream);
  
  // check if the header was read correctly
  REQUIRE( l_read.size() == 4 ); // 4 columns
  REQUIRE( l_read[0].first == "x" );
  REQUIRE( l_read[1].first == "y" );
  REQUIRE( l_read[2].first == "height" );
  REQUIRE( l_read[3].first == "momentum_x" );
  
  // check if the length of read data is correct
  for(t_idx l_i=0;l_i<4;l_i++){
    REQUIRE( l_read[l_i].second.size() == 5 );
  }
  
  t_real l_data[20] = {
    0.25,0,1,5,
    0.75,0,2,4,
    1.25,0,3,3,
    1.75,0,4,2,
    2.25,0,5,1,
  };
  
  // check if all data was read correctly
  for(t_idx l_j=0;l_j<4;l_j++){
    for(t_idx l_i=0;l_i<5;l_i++){
      REQUIRE( l_read[l_j].second[l_i] == Approx(l_data[l_i*4+l_j]) );
    }
  }
  
}
