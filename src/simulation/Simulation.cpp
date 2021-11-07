/**
 * @author Antonio Noack
 * 
 * @section DESCRIPTION
 * Simulation functions.
 **/

#include "Simulation.h"
#include "../patches/WavePropagation1d.h"
#include "../io/Csv.h"

#include <iostream>
#include <fstream>
#include <algorithm> // std::max
#include <cstdio> // delete old files

void tsunami_lab::simulation::Simulation::run( tsunami_lab::t_idx           l_nx,
                                               tsunami_lab::t_idx           l_ny,
                                               tsunami_lab::setups::Setup & l_setup,
                                               tsunami_lab::t_real          l_setupScale,
                                               tsunami_lab::t_real          l_cellSizeMeters,
                                               tsunami_lab::t_idx           l_numTimesteps,
                                               tsunami_lab::t_idx           l_outputStepSize,
                                               tsunami_lab::t_idx           l_numOutputSteps) {
   
  std::cout << "runtime configuration" << std::endl;
  std::cout << "  number of cells in x-direction: " << l_nx << std::endl;
  std::cout << "  number of cells in y-direction: " << l_ny << std::endl;
  std::cout << "  cell size (meters):             " << l_cellSizeMeters << std::endl;
  
  for(tsunami_lab::t_idx l_i=0;l_i<1000;l_i++){
    // delete all old files;
    // ParaView caused me enough headaches
    std::string l_path = "solution_" + std::to_string(l_i) + ".csv";
    remove(l_path.c_str());
  }

  // construct solver
  tsunami_lab::patches::WavePropagation1d l_waveProp( l_nx, &l_setup, l_setupScale );

  // set up print control
  tsunami_lab::t_idx  l_nOut = 0;
  tsunami_lab::t_real l_timestep;
  
  // currently, this is often a pretty good amount of timesteps
  // it may be changed in the future
  tsunami_lab::t_real l_simTime = 0;

  std::cout << "entering time loop" << std::endl;

  // iterate over time
  tsunami_lab::t_idx l_lastOutputIndex = -1;// -1 to always print the initial state, 0 to skip it
  for(tsunami_lab::t_idx l_timeStepIndex=0; l_timeStepIndex < l_numTimesteps; l_timeStepIndex++ ){
    // index, which frame we'd need to print theoretically
    // if there are more frames requested than simulated, we just skip some
    tsunami_lab::t_idx l_outputIndex = l_timeStepIndex * (l_numOutputSteps-1) / std::max(l_numTimesteps-1, (t_idx) 1);
    if(l_lastOutputIndex != l_outputIndex) {
      l_lastOutputIndex = l_outputIndex;
      
      std::cout << "  simulation time: " << l_simTime << ", #time steps: "<< l_timeStepIndex << std::endl;

      std::string l_path = "solution_" + std::to_string(l_nOut) + ".csv";
      std::cout << "  writing wave field to " << l_path << std::endl;

      std::ofstream l_file;
      l_file.open(l_path);

      tsunami_lab::io::Csv::write( l_cellSizeMeters, l_nx, l_ny, l_outputStepSize, l_waveProp.getStride(), l_waveProp.getHeight(), l_waveProp.getMomentumX(), nullptr, l_waveProp.getBathymetry(), l_file );
      l_file.close();
      l_nOut++;
    }

    l_timestep = l_waveProp.computeMaxTimestep(l_cellSizeMeters);
    l_waveProp.setGhostOutflow();
    
    tsunami_lab::t_real l_scaling = l_timestep / l_cellSizeMeters;
    l_waveProp.timeStep(l_scaling);

    l_simTime += l_timestep;
  }
  
  std::cout << "finished time loop" << std::endl;
  
}
