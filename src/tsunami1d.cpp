/**
 * @author Antonio Noack
 * @section DESCRIPTION
 * Entry-point for tsunami simulation from a csv file in 1d.
 **/

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <fstream>
#include <limits>
#include <algorithm> // std::max
#include <cstdio> // delete old files
 
#include "patches/WavePropagation1d.h"
#include "setups/TsunamiEvent1d.h"
#include "io/Csv.h"

#define t_real tsunami_lab::t_real
#define t_idx  tsunami_lab::t_idx

int main( int   i_argc, char *i_argv[] ) {
  
  // this sample simulates a tsunami event from a CSV file
  
  t_real l_cellSizeMeters = 1;// is computed from the file: (last.x - first.x) / (last.index - first.index)
  t_idx  l_numTimesteps = 1000;
  t_idx  l_numOutputSteps = 100;
  t_real l_scale = 1;
  
  t_real l_displacement = -10;
  
  std::string l_fileName = "./data/fukushima.csv";

  std::cout << "####################################" << std::endl;
  std::cout << "###   Tsunami Event Simulation   ###" << std::endl;
  std::cout << "####################################" << std::endl;

  if( i_argc < 2 ) {
    std::cerr << "invalid number of arguments, usage:" << std::endl;
    std::cerr << "  ./build/tsunami1d NUM_TIMESTEPS NUM_OUTPUT_STEPS INPUT_FILE" << std::endl;
    std::cerr << "    where" << std::endl;
    std::cerr << "      NUM_TIMESTEPS is the number of time steps," << std::endl;
    std::cerr << "      NUM_OUTPUT_STEPS is the number of time steps, that are printed to csv files," << std::endl;
    std::cerr << "      INPUT_FILE is the csv file containing the event data; default = ./data/fukushima.csv ." << std::endl;
    return EXIT_FAILURE;
  } else {
    if(i_argc > 1) l_numTimesteps = atol( i_argv[1] );
    if(i_argc > 2) l_numOutputSteps = atol(i_argv[2]);
    if(i_argc > 3) l_fileName = i_argv[3];
  }
  
  // todo load the data from the csv file
  auto l_inputStream = std::ifstream(l_fileName, std::ios::in);
  if(!l_inputStream){
    std::cerr << "could not open file '" << l_fileName << "'" << std::endl;
    return EXIT_FAILURE;
  }
  
  auto l_loadedData = tsunami_lab::io::Csv::read(l_inputStream);
  l_inputStream.close();// close the file
  
  // here we probably copy; this potentially could be optimized
  auto l_xs      = tsunami_lab::io::Csv::findColumn(l_loadedData, "track_location");
  auto l_heights = tsunami_lab::io::Csv::findColumn(l_loadedData, "height");
  
  if(l_xs.size() < 1){
    std::cerr << "did not find position data in the file" << std::endl;
    return EXIT_FAILURE;
  } else if(l_heights.size() < 1){
    std::cerr << "did not find bathymetry data in the file" << std::endl;
    return EXIT_FAILURE;
  } else if(l_xs.size () != l_heights.size()){
    std::cerr << "the sizes of positions and heights don't match" << std::endl;
    return EXIT_FAILURE;
  }
  
  l_cellSizeMeters = (l_xs.back() - l_xs.front()) / (l_xs.size() - 1) / l_scale;
  
  t_idx l_nx = l_xs.size() * l_scale - 2;// 2 = ghost cells
  t_idx l_ny = 1;// it's just a 1d simulation
  
  std::cout << "runtime configuration" << std::endl;
  std::cout << "  number of cells in x-direction: " << l_nx << std::endl;
  std::cout << "  cell size (meters):             " << l_cellSizeMeters << std::endl;
  
  for(t_idx l_i=0;l_i<1000;l_i++){
	// delete all old files;
	// paraview caused me enough headaches
	std::string l_path = "solution_" + std::to_string(l_i) + ".csv";
	remove(l_path.c_str());
  }
  
  t_real l_displacementStart = 175000 / l_cellSizeMeters;
  t_real l_displacementEnd   = 250000 / l_cellSizeMeters;

  // construct setup
  tsunami_lab::setups::Setup *l_setup = new tsunami_lab::setups::TsunamiEvent1d( l_heights.data(), l_heights.size(), 1/l_scale, l_displacementStart, l_displacementEnd, l_displacement );
  
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
  t_idx l_lastOutputIndex = -1;// -1 to always print the initial state, 0 to skip it
  for(t_idx l_timeStepIndex=0; l_timeStepIndex < l_numTimesteps; l_timeStepIndex++ ){
    // index, which frame we'd need to print theoretically
    // if there are more frames requested than simulated, we just skip some
    t_idx l_outputIndex = l_timeStepIndex * (l_numOutputSteps-1) / std::max(l_numTimesteps-1, (t_idx) 1);
    if(l_lastOutputIndex != l_outputIndex) {
      l_lastOutputIndex = l_outputIndex;
      
      std::cout << "  simulation time: " << l_simTime << ", #time steps: "<< l_timeStepIndex << std::endl;

      std::string l_path = "solution_" + std::to_string(l_nOut) + ".csv";
      std::cout << "  writing wave field to " << l_path << std::endl;

      std::ofstream l_file;
      l_file.open( l_path  );

      tsunami_lab::io::Csv::write( l_cellSizeMeters, l_nx, l_ny, l_scale, l_waveProp->getStride(), l_waveProp->getHeight(), l_waveProp->getMomentumX(), nullptr, l_waveProp->getBathymetry(), l_file );
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
