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
 * Entry-point for simulations.
 **/
#include "simulation/Simulation.h"
#include "setups/Discontinuity1d.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <fstream>
#include <limits>
#include <algorithm> // std::max
#include <cstdio> // delete old files

#define t_real tsunami_lab::t_real
#define t_idx  tsunami_lab::t_idx

int main( int   i_argc,
          char *i_argv[] ) {
              
  // this sample once expected a size, and then implemented the dam problem
  // instead, I'll change it to use my Discontinuity1d setup to test the shock-shock and rare-rare problems.
  
  // arguments: <number of cells> <left height> <right height> <left impulse> <right impulse>
  
  // number of cells in x- and y-direction
  t_idx  l_nx = 0;
  t_idx  l_ny = 1;
  t_real l_heightLeft = 10;
  t_real l_heightRight = 5;
  t_real l_impulseLeft = 0;
  t_real l_impulseRight = 0;
  t_real l_cellSizeMeters = 1;
  t_idx  l_numTimesteps = 1000;
  t_idx  l_numOutputSteps = 100;

  std::cout << "####################################" << std::endl;
  std::cout << "### Tsunami Lab                  ###" << std::endl;
  std::cout << "###                              ###" << std::endl;
  std::cout << "### https://scalable.uni-jena.de ###" << std::endl;
  std::cout << "####################################" << std::endl;

  if( i_argc < 2 ) {
    std::cerr << "invalid number of arguments, usage:" << std::endl;
    std::cerr << "  ./build/tsunami_lab N_CELLS_X HEIGHT_LEFT HEIGHT_RIGHT IMPULSE_LEFT IMPULSE_RIGHT" << std::endl;
    std::cerr << "    where" << std::endl;
    std::cerr << "      N_CELLS_X is the number of cells in x-direction," << std::endl;
    std::cerr << "      HEIGHT_LEFT is the height of the water on the left half, default = 10," << std::endl;
    std::cerr << "      HEIGHT_RIGHT is the height of the water on the right half, default = 5," << std::endl;
    std::cerr << "      IMPULSE_LEFT is the impulse on the left side, default = 0," << std::endl;
    std::cerr << "      IMPULSE_RIGHT is the impulse on the right side, default = 0." << std::endl;
    return EXIT_FAILURE;
  } else {
    l_nx = atoi( i_argv[1] );
    if( l_nx < 1 ) {
      std::cerr << "invalid number of cells" << std::endl;
      return EXIT_FAILURE;
    }
    if(i_argc > 2) l_heightLeft   = atof(i_argv[2]);
    if(i_argc > 3) l_heightRight  = atof(i_argv[3]);
    if(i_argc > 4) l_impulseLeft  = atof(i_argv[4]);
    if(i_argc > 5) l_impulseRight = atof(i_argv[5]);
  }
  
  // construct setup
  tsunami_lab::setups::Discontinuity1d l_setup( l_heightLeft, l_heightRight, l_impulseLeft, l_impulseRight, l_nx / 2 );
  
  // run simulation
  tsunami_lab::simulation::Simulation::run( l_nx, l_ny, l_setup, 1, l_cellSizeMeters, l_numTimesteps, 1, l_numOutputSteps );
  
  std::cout << "finished, exiting" << std::endl;
  return EXIT_SUCCESS;
}
