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
 
#include <vector>
#include <map>
#include <fstream>
#include <limits> // infinity, max int

#include <yaml-cpp/yaml.h>

#include "patches/WavePropagation1d.h"
#include "patches/WavePropagation2d.h"
#include "setups/Discontinuity1d.h"
#include "setups/DamBreak1d.h"
#include "setups/DamBreak2d.h"
#include "setups/SubcriticalFlow1d.h"
#include "setups/SupercriticalFlow1d.h"
#include "setups/TsunamiEvent1d.h"
#include "io/Station.h"
#include "io/Csv.h"

#define t_real tsunami_lab::t_real
#define t_idx  tsunami_lab::t_idx

std::string trim(std::string l_s) {// somehow not part of the C++ standard
  auto l_whitespace = " \t\n\r\f\v";
  auto l_i0 = l_s.find_first_not_of(l_whitespace);
  if(l_i0 == std::string::npos) return "";
  auto l_i1 = l_s.find_last_not_of(l_whitespace) + 1;
  return l_s.substr(l_i0, l_i1 - l_i0);
}

bool readBoolean(YAML::Node &i_config, std::string i_key, bool i_defaultValue) {
  if(!i_config[i_key]) return i_defaultValue;
  return i_config[i_key].as<bool>();
}

t_idx readInt(YAML::Node &i_config, std::string i_key, t_idx i_defaultValue) {
  if(!i_config[i_key]) return i_defaultValue;
  return i_config[i_key].as<t_idx>();
}

t_real readFloat(YAML::Node &i_config, std::string i_key, t_real i_defaultValue) {
  if(!i_config[i_key]) return i_defaultValue;
  return i_config[i_key].as<t_real>();
}

int main( int i_argc, char *i_argv[] ) {
  
  std::cout << "####################################" << std::endl;
  std::cout << "### Tsunami Lab                  ###" << std::endl;
  std::cout << "###                              ###" << std::endl;
  std::cout << "### https://scalable.uni-jena.de ###" << std::endl;
  std::cout << "####################################" << std::endl;
  
  if( i_argc != 2 ) {
    std::cerr << "invalid number of arguments, usage:" << std::endl;
    std::cerr << "  ./build/tsunami_lab YAML_CONFIG_FILE" << std::endl;
    std::cerr << "    where YAML_CONFIG_FILE is a file containing the runtime configuration as yaml" << std::endl;
    return EXIT_FAILURE;
  }
  
  auto l_config = YAML::LoadFile(i_argv[1]);
  // if there are more params from the console, we could add them to the config as well
  
  std::vector<tsunami_lab::io::Station> l_stations;
  tsunami_lab::setups::Setup* l_setup = nullptr;
  
  // many setups share similar setup information, so just request them all at once;
  // it should not matter for performance, as this all will be done within 1ms.
  t_idx  l_nx              =   readInt(l_config, "nx",  1);
  t_idx  l_ny              =   readInt(l_config, "ny",  1);
  t_real l_heightLeft      = readFloat(l_config, "hl", 10);
  t_real l_heightRight     = readFloat(l_config, "hr",  5);
  t_real l_impulseLeft     = readFloat(l_config, "hul", 0);
  t_real l_impulseRight    = readFloat(l_config, "hur", 0);
  t_real l_bathymetryLeft  = readFloat(l_config, "bl", -1);
  t_real l_bathymetryRight = readFloat(l_config, "br", -1);
  t_real l_splitPositionX  = readFloat(l_config, "splitPositionX", l_nx * 0.5);
  t_real l_splitPositionY  = readFloat(l_config, "splitPositionY", l_ny * 0.5);
  t_real l_damRadius       = readFloat(l_config, "damRadius", l_nx / 4);
  t_real l_cellSizeMeters  = readFloat(l_config, "cellSize", 1);// size of a single cell in meters
  t_idx  l_maxTimesteps    =   readInt(l_config, "maxSteps", std::numeric_limits<t_idx>::max());// max number of simulation timesteps
  t_real l_maxDuration     = readFloat(l_config, "maxDuration", std::numeric_limits<t_real>::infinity());// max simulation time in seconds
  t_real l_outputPeriod    = readFloat(l_config, "outputPeriod", 1);// every n th second, a result file will be written
  t_idx  l_outputStepSize  =   readInt(l_config, "outputStepSize", 1);// n = every nth field is written to disk; saves storage space
  t_real l_scale           = readFloat(l_config, "scale", 1);// scales the setup for accuracy or performance
  
  // we could/should write a few warnings to console, if any of these happen
  if(l_nx < 1) l_nx = 1;
  if(l_ny < 1) l_ny = 1;
  if(l_outputStepSize < 1) l_outputStepSize = 1;
  if(l_heightLeft < 0) l_heightLeft = 0;
  if(l_heightRight < 0) l_heightRight = 0;
  if(l_cellSizeMeters < 0) l_cellSizeMeters = 1;
  
  if(l_maxTimesteps == std::numeric_limits<t_idx>::max() && l_maxDuration == std::numeric_limits<t_real>::infinity()){
    std::cerr << "you need to specify a time limit or timestep limit (maxDuration, maxSteps)" << std::endl;
    return EXIT_FAILURE;
  }
  
  bool l_printStationComments = readBoolean(l_config, "printStationComments", true);
  if(l_config["stations"]) {
    auto   l_stationData = l_config["stations"].as<std::vector<YAML::Node>>();
    t_real l_delayBetweenRecords = readFloat(l_config, "delayBetweenRecords", 1);
    for(t_idx i=0;i<l_stationData.size();i++){
      auto l_station = l_stationData[i];
      std::string l_name = l_station["name"].as<std::string>();
      t_idx       l_x    = l_station["x"].as<t_idx>();
      t_idx       l_y    = l_station["y"].as<t_idx>();
      tsunami_lab::io::Station l_station1(l_x, l_y, l_name, l_delayBetweenRecords);
      l_stations.push_back(l_station1);
    }
  }
  
  std::string l_setupName = "Discontinuity1d";
  if(l_config["setup"]) l_setupName = l_config["setup"].as<std::string>();
  else std::cout << "using default setup: " << l_setupName << std::endl;
  
  // must stay in scope; I don't have a better solution currently; maybe the value could be moved into setup
  // a copy shouldn't be expensive compared to our simulation and writing the results to disk
  std::vector<t_real> l_heights;
  
  // C++ doesn't have switch-case for strings? :/
  if(
    l_setupName == "DamBreak" ||
    l_setupName == "DamBreak1d"
  ) {
    l_setup = new tsunami_lab::setups::DamBreak1d(l_heightLeft, l_heightRight, l_splitPositionX);
  } else if(
    l_setupName == "DamBreakCircle" ||
    l_setupName == "DamBreak2d"
  ) {
    auto l_dambreak2d = new tsunami_lab::setups::DamBreak2d(l_heightLeft, l_heightRight, l_splitPositionX, l_splitPositionY, l_damRadius);
    if(l_config["obstacleBathymetry"]) {
      // an obstacle should be defined
      l_dambreak2d->setObstacle(
        readFloat(l_config, "obstacleX0", 0),
        readFloat(l_config, "obstacleX1", 0), 
        readFloat(l_config, "obstacleY0", 0),
        readFloat(l_config, "obstacleY1", 0), 
        readFloat(l_config, "obstacleBathymetry", 0)
      );
    }
    l_setup = l_dambreak2d;
  } else if(
    l_setupName == "Discontinuity1d" ||
    l_setupName == "Discontinuity"
  ){
    l_setup = new tsunami_lab::setups::Discontinuity1d(l_heightLeft, l_heightRight, l_impulseLeft, l_impulseRight, l_bathymetryLeft, l_bathymetryRight, l_splitPositionX);
  } else if(
    l_setupName == "Track" || 
    l_setupName == "Tsunami"
  ) {
    if(l_config["setupFile"]){
      
      // load the data from the csv file
      auto l_loadedData = tsunami_lab::io::Csv::read(l_config["setupFile"].as<std::string>());
      
      // here we probably copy; this potentially could be optimized
      auto l_xs = tsunami_lab::io::Csv::findColumn(l_loadedData, "track_location");
      l_heights = tsunami_lab::io::Csv::findColumn(l_loadedData, "height");
      
      if(l_xs.size() < 1){
        std::cerr << "did not find position data in track file" << std::endl;
        return EXIT_FAILURE;
      } else if(l_heights.size() < 1){
        std::cerr << "did not find bathymetry data in track file" << std::endl;
        return EXIT_FAILURE;
      }
      
      l_cellSizeMeters = (l_xs.back() - l_xs.front()) / (l_xs.size() - 1) / l_scale;
      
      l_nx = l_xs.size() * l_scale - 2;// 2 = ghost cells
      l_ny = 1;// it's just a 1d simulation
      
      t_real l_displacementStart = readFloat(l_config, "displacementStart", 175000) / l_cellSizeMeters;
      t_real l_displacementEnd   = readFloat(l_config, "displacementEnd",   250000) / l_cellSizeMeters;
      t_real l_displacement      = readFloat(l_config, "displacement", 10);
      
      // construct setup
      l_setup = new tsunami_lab::setups::TsunamiEvent1d(l_heights.data(), l_heights.size(), 1/l_scale, l_displacementStart, l_displacementEnd, l_displacement);
      
    } else {
      std::cerr << "missing parameter 'setupFile' for setup type Track" << std::endl;
      return EXIT_FAILURE;
    }
  } else if( 
    l_setupName == "Subcritical" ||
    l_setupName == "SubcriticalFlow" ||
    l_setupName == "SubcriticalFlow1d"
  ) {
    l_nx = (t_idx) (25 * l_scale);
    l_cellSizeMeters = 1 / l_scale;
    l_setup = new tsunami_lab::setups::SubcriticalFlow1d();
  } else if(
    l_setupName == "Supercritical" ||
    l_setupName == "SupercriticalFlow" ||
    l_setupName == "SupercriticalFlow1d"
  ) {
    l_nx = (t_idx) (25 * l_scale);
    l_cellSizeMeters = 1 / l_scale;
    l_setup = new tsunami_lab::setups::SupercriticalFlow1d();
  } else {
    std::cerr << "unknown/invalid setup type \"" << l_setupName << "\"" << std::endl;
    return EXIT_FAILURE;
  }
  
  // run simulation
  std::cout << "runtime configuration" << std::endl;
  std::cout << "  max timesteps:                  " << l_maxTimesteps << std::endl;
  std::cout << "  max duration:                   " << l_maxDuration << std::endl;
  std::cout << "  number of cells in x-direction: " << l_nx << std::endl;
  std::cout << "  number of cells in y-direction: " << l_ny << std::endl;
  std::cout << "  cell size (meters):             " << l_cellSizeMeters << std::endl;
  std::cout << "  number of stations:             " << l_stations.size() << std::endl;
  
  for(t_idx l_i=0;l_i<1000;l_i++){
    // delete all old files;
    // ParaView caused me enough headaches
    std::string l_path = "solution_" + std::to_string(l_i) + ".csv";
    remove(l_path.c_str());
  }

  // construct solver
  tsunami_lab::patches::WavePropagation* l_waveProp;
  if(l_ny <= 1){
    l_waveProp = new tsunami_lab::patches::WavePropagation1d(l_nx, l_setup, l_scale);
  } else {
    l_waveProp = new tsunami_lab::patches::WavePropagation2d(l_nx, l_ny, l_setup, l_scale, l_scale, readFloat(l_config, "cflFactor", 0.45));
  }

  // set up print control
  t_idx  l_nOut = 0;
  t_real l_timestep;
  
  double l_time = 0;

  std::cout << "entering time loop" << std::endl;

  // iterate over time
  t_idx l_lastOutputIndex = -1;
  t_idx l_timeStepIndex = 0;
  for(; l_timeStepIndex < l_maxTimesteps && l_time < l_maxDuration; l_timeStepIndex++ ){

    // index, which frame we'd need to print theoretically
    // if there are more frames requested than simulated, we just skip some
    t_idx l_outputIndex = (t_idx) (l_time / l_outputPeriod);
    if(l_lastOutputIndex != l_outputIndex) {
      l_lastOutputIndex = l_outputIndex;
      
      std::cout << "  simulation time: " << l_time << ", #time steps: "<< l_timeStepIndex << std::endl;

      std::string l_path = "solution_" + std::to_string(l_nOut) + ".csv";
      std::cout << "  writing wave field to " << l_path << std::endl;

      std::ofstream l_file(l_path, std::ios::out);
      tsunami_lab::io::Csv::write(l_cellSizeMeters, l_nx, l_ny, l_outputStepSize, l_waveProp->getStride(), l_waveProp->getHeight(), l_waveProp->getMomentumX(), l_waveProp->getMomentumY(), l_waveProp->getBathymetry(), l_file);
      l_file.close();
      l_nOut++;
    }
    
    // update recording stations, if there are any
    if(!l_stations.empty() && l_stations[0].needsUpdate(l_time)) {
      for(auto &l_station : l_stations) {
        l_station.recordState(*l_waveProp, l_time);
      }
    }

    l_timestep = l_waveProp->computeMaxTimestep(l_cellSizeMeters);
    l_waveProp->setGhostOutflow();
    
    t_real l_scaling = l_timestep / l_cellSizeMeters;
    l_waveProp->timeStep(l_scaling);

    l_time += l_timestep;
  }
  
  std::cout << "finished time loop" << std::endl;
  
  // print last state
  std::cout << "  simulation time: " << l_time << ", #time steps: "<< l_timeStepIndex << std::endl;

  std::string l_path = "solution_" + std::to_string(l_nOut) + ".csv";
  std::cout << "  writing wave field to " << l_path << std::endl;

  std::ofstream l_file(l_path, std::ios::out);
  tsunami_lab::io::Csv::write( l_cellSizeMeters, l_nx, l_ny, l_outputStepSize, l_waveProp->getStride(), l_waveProp->getHeight(), l_waveProp->getMomentumX(), l_waveProp->getMomentumY(), l_waveProp->getBathymetry(), l_file );
  l_file.close();
  
  for(auto &l_station : l_stations){
    l_station.write(l_printStationComments);
  }
  
  std::cout << "finished writing last state" << std::endl;
  
  delete l_waveProp;
  
  delete l_setup;
  
  return EXIT_SUCCESS;
}
