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
 
#include "simulation/Simulation.h"
#include "setups/TsunamiEvent1d.h"
#include "io/Csv.h"

#define t_real tsunami_lab::t_real
#define t_idx  tsunami_lab::t_idx

int main( int i_argc, char *i_argv[] ) {
  
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
  
  t_real l_displacementStart = 175000 / l_cellSizeMeters;
  t_real l_displacementEnd   = 250000 / l_cellSizeMeters;

  // construct setup
  tsunami_lab::setups::TsunamiEvent1d l_setup( l_heights.data(), l_heights.size(), 1/l_scale, l_displacementStart, l_displacementEnd, l_displacement );
  
  // run the simulation
  tsunami_lab::simulation::Simulation::run( l_nx, l_ny, l_setup, 1, l_cellSizeMeters, l_numTimesteps, l_scale, l_numOutputSteps );
  
  std::cout << "finished, exiting" << std::endl;
  return EXIT_SUCCESS;
}
