/**
 * @author Antonio Noack
 *
 * @section DESCRIPTION
 * Subcritical & Supercritical simulation runs.
 **/
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <fstream>
#include <limits>
#include <algorithm> // std::max
#include <cstdio> // delete old files

#include "simulation/Simulation.h"
#include "setups/SubcriticalFlow1d.h"
#include "setups/SupercriticalFlow1d.h"

#define t_real tsunami_lab::t_real
#define t_idx  tsunami_lab::t_idx

int main( int   i_argc,
          char *i_argv[] ) {
  
  t_idx  l_ny = 1;
  t_idx  l_type;
  t_idx  l_numTimesteps = 1000;
  t_idx  l_numOutputSteps = 100;
  t_real l_scale = 10;

  if( i_argc < 2 ) {
    std::cerr << "invalid number of arguments, usage:" << std::endl;
    std::cerr << "  ./build/tsunami_lab TYPE TIMESTEPS OUTPUT_STEPS" << std::endl;
    std::cerr << "    where" << std::endl;
    std::cerr << "      TYPE: 0 = subcritical, 1 = supercritical," << std::endl;
    std::cerr << "      TIMESTEPS number of time steps, default = 1000," << std::endl;
    std::cerr << "      OUTPUT_STEPS number of written time steps, default = 100." << std::endl;
    return EXIT_FAILURE;
  } else {
    l_type = atoi( i_argv[1] );
    if( l_type != 0 && l_type != 1 ) {
      std::cerr << "invalid type" << std::endl;
      return EXIT_FAILURE;
    }
    if(i_argc > 2) l_numTimesteps   = atoi(i_argv[2]);
    if(i_argc > 3) l_numOutputSteps = atoi(i_argv[3]);
  }
  
  t_idx l_nx = (t_idx) (25 * l_scale);
  t_real l_cellSizeMeters = 1 / l_scale;
  
  // construct setups
  tsunami_lab::setups::Setup* l_setup;
  if( l_type == 0 ) l_setup = new tsunami_lab::setups::SubcriticalFlow1d();
  else              l_setup = new tsunami_lab::setups::SupercriticalFlow1d();
  
  // run simulation
  tsunami_lab::simulation::Simulation::run( l_nx, l_ny, *l_setup, 1/l_scale, l_cellSizeMeters, l_numTimesteps, 1, l_numOutputSteps );
  
  delete l_setup;
  
  std::cout << "finished, exiting" << std::endl;
  return EXIT_SUCCESS;
}
