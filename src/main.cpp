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
 * Entry-point for simulations.
 **/
#include "patches/WavePropagation1d.h"
// #include "setups/DamBreak1d.h"
#include "setups/Discontinuity1d.h"
#include "io/Csv.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <fstream>
#include <limits>

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
  t_real l_numTimesteps = 1000;
  t_real l_numOutputSteps = 100;

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
  
  std::cout << "runtime configuration" << std::endl;
  std::cout << "  number of cells in x-direction: " << l_nx << std::endl;
  std::cout << "  number of cells in y-direction: " << l_ny << std::endl;
  std::cout << "  cell size (meters):             " << l_cellSizeMeters << std::endl;

  // construct setup
  tsunami_lab::setups::Setup *l_setup = new tsunami_lab::setups::Discontinuity1d( l_heightLeft, l_heightRight, l_impulseLeft, l_impulseRight, l_nx / 2 );
  // construct solver
  tsunami_lab::patches::WavePropagation1d *l_waveProp;
  l_waveProp = new tsunami_lab::patches::WavePropagation1d( l_nx, l_setup, 1.0 );

  // set up print control
  t_idx  l_nOut = 0;
  t_real l_timestep;
  
  // currently, this is often a pretty good amount of timesteps
  // it may be changed in the future
  t_real l_simTime = 0;

  std::cout << "entering time loop" << std::endl;

  // iterate over time
  t_idx l_lastOutputIndex = -1;// -1 to always print the inital state, 0 to skip it
  for(t_idx l_timeStepIndex=0; l_timeStepIndex < l_numTimesteps; l_timeStepIndex++ ){
    // index, which frame we'd need to print theoretically
    // if there are more frames requested than simulated, we just skip some
    t_idx l_outputIndex = l_timeStepIndex * l_numOutputSteps / l_numTimesteps;
    if(l_lastOutputIndex != l_outputIndex) {
      l_lastOutputIndex = l_outputIndex;
      
      std::cout << "  simulation time: " << l_simTime << ", #time steps: "<< l_timeStepIndex << std::endl;

      std::string l_path = "solution_" + std::to_string(l_nOut) + ".csv";
      std::cout << "  writing wave field to " << l_path << std::endl;

      std::ofstream l_file;
      l_file.open( l_path  );

      tsunami_lab::io::Csv::write( l_cellSizeMeters, l_nx, l_ny, l_waveProp->getStride(), l_waveProp->getHeight(), l_waveProp->getMomentumX(), nullptr, l_file );
      l_file.close();
      l_nOut++;
    }

    l_timestep = l_waveProp->computeMaxTimestep(l_cellSizeMeters);
    l_waveProp->setGhostOutflow();
    
    t_real l_scaling = l_timestep / l_cellSizeMeters;
    l_waveProp->timeStep(l_scaling);

    l_simTime += l_timestep;
  }
  
  std::cout << "finished time loop" << std::endl;
  
  // free memory
  std::cout << "freeing memory" << std::endl;
  delete l_setup;
  delete l_waveProp;
  
  std::cout << "finished, exiting" << std::endl;
  return EXIT_SUCCESS;
}
