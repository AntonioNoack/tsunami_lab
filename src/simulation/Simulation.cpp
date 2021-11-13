/**
 * @author Antonio Noack
 * 
 * @section DESCRIPTION
 * Simulation functions.
 **/

#include "Simulation.h"
#include "../io/Csv.h"
#include "../patches/WavePropagation.h"
#include "../patches/WavePropagation1d.h"
#include "../patches/WavePropagation2d.h"

#include <iostream>
#include <fstream>
#include <algorithm> // std::max
#include <cstdio> // delete old files
#include <vector>
#include <cmath> // ceil

void tsunami_lab::simulation::Simulation::run( tsunami_lab::t_idx                      i_nx,
                                               tsunami_lab::t_idx                      i_ny,
                                               tsunami_lab::setups::Setup            & i_setup,
                                               tsunami_lab::t_real                     i_setupScale,
                                               tsunami_lab::t_real                     i_cellSizeMeters,
                                               tsunami_lab::t_idx                      i_maxTimesteps,
											   tsunami_lab::t_real                     i_maxDuration,
                                               tsunami_lab::t_idx                      i_outputStepSize,
                                               tsunami_lab::t_idx                      i_numOutputSteps,
                                               std::vector<tsunami_lab::io::Station> & i_stations ) {
   
  std::cout << "runtime configuration" << std::endl;
  std::cout << "  max timesteps:                  " << i_maxTimesteps << std::endl;
  std::cout << "  max duration:                   " << i_maxDuration << std::endl;
  std::cout << "  number of cells in x-direction: " << i_nx << std::endl;
  std::cout << "  number of cells in y-direction: " << i_ny << std::endl;
  std::cout << "  cell size (meters):             " << i_cellSizeMeters << std::endl;
  std::cout << "  number of stations:             " << i_stations.size() << std::endl;
  
  for(tsunami_lab::t_idx l_i=0;l_i<1000;l_i++){
    // delete all old files;
    // ParaView caused me enough headaches
    std::string l_path = "solution_" + std::to_string(l_i) + ".csv";
    remove(l_path.c_str());
  }

  // construct solver
  tsunami_lab::patches::WavePropagation* l_waveProp;
  if(i_ny <= 1){
    l_waveProp = new tsunami_lab::patches::WavePropagation1d( i_nx, &i_setup, i_setupScale );
  } else {
    l_waveProp = new tsunami_lab::patches::WavePropagation2d( i_nx, i_ny, &i_setup, i_setupScale, i_setupScale );
  }

  // set up print control
  tsunami_lab::t_idx  l_nOut = 0;
  tsunami_lab::t_real l_timestep;
  
  double l_time = 0;
  
  t_idx l_maxTimestepsForPrinting = i_maxTimesteps;
  if(l_maxTimestepsForPrinting >= std::numeric_limits<t_idx>::max()){
	// if no timestep limit is set, then guess a limit for writing files to disk,
	// otherwise nearly no files would be written
	t_real l_guessedDt = l_waveProp->computeMaxTimestep( i_cellSizeMeters );
	l_maxTimestepsForPrinting = (t_idx) std::ceil(i_maxDuration / l_guessedDt);
  }

  std::cout << "entering time loop" << std::endl;

  // iterate over time
  tsunami_lab::t_idx l_lastOutputIndex = -1;
  for(tsunami_lab::t_idx l_timeStepIndex = 0; l_timeStepIndex < i_maxTimesteps && l_time < i_maxDuration; l_timeStepIndex++ ){

    // index, which frame we'd need to print theoretically
    // if there are more frames requested than simulated, we just skip some
    tsunami_lab::t_idx l_outputIndex = l_timeStepIndex * (i_numOutputSteps-1) / std::max(l_maxTimestepsForPrinting-1, (t_idx) 1);
    if(l_lastOutputIndex != l_outputIndex) {
      l_lastOutputIndex = l_outputIndex;
      
      std::cout << "  simulation time: " << l_time << ", #time steps: "<< l_timeStepIndex << std::endl;

      std::string l_path = "solution_" + std::to_string(l_nOut) + ".csv";
      std::cout << "  writing wave field to " << l_path << std::endl;

      std::ofstream l_file;
      l_file.open( l_path );

      tsunami_lab::io::Csv::write( i_cellSizeMeters, i_nx, i_ny, i_outputStepSize, l_waveProp->getStride(), l_waveProp->getHeight(), l_waveProp->getMomentumX(), nullptr, l_waveProp->getBathymetry(), l_file );
      l_file.close();
      l_nOut++;
    }
    
    // update recording stations, if there are any
    if(!i_stations.empty() && i_stations[0].needsUpdate( l_time )) {
      for(auto &l_station : i_stations) {
        l_station.recordState( *l_waveProp, l_time );
      }
    }

    l_timestep = l_waveProp->computeMaxTimestep( i_cellSizeMeters );
    l_waveProp->setGhostOutflow();
    
    tsunami_lab::t_real l_scaling = l_timestep / i_cellSizeMeters;
    l_waveProp->timeStep( l_scaling );

    l_time += l_timestep;
  }
  
  std::cout << "finished time loop" << std::endl;
  
  for(auto &l_station : i_stations){
    l_station.write();
  }
  
  std::cout << "finished writing stations" << std::endl;
  
  delete l_waveProp;
  
}
