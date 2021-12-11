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
#include <chrono> // measuring performance
#include <sys/stat.h> // check whether a file exists
#include <omp.h> // for max threads
#include <cmath> // std::sqrt

#include <yaml-cpp/yaml.h>

#include "io/Csv.h"
#include "io/NetCdf.h"
#include "io/Station.h"
#include "patches/WavePropagation1d.h"
#include "patches/WavePropagation2d.h"
#include "setups/ArtificialTsunami2d.h"
#include "setups/DamBreak1d.h"
#include "setups/DamBreak2d.h"
#include "setups/Discontinuity1d.h"
#include "setups/SubcriticalFlow1d.h"
#include "setups/SupercriticalFlow1d.h"
#include "setups/TsunamiEvent1d.h"
#include "setups/TsunamiEvent2d.h"

#define t_real tsunami_lab::t_real
#define t_idx  tsunami_lab::t_idx

std::string trim(std::string l_s) {// somehow not part of the C++ standard
  auto l_whitespace = " \t\n\r\f\v";
  auto l_i0 = l_s.find_first_not_of(l_whitespace);
  if(l_i0 == std::string::npos) return "";
  auto l_i1 = l_s.find_last_not_of(l_whitespace) + 1;
  return l_s.substr(l_i0, l_i1 - l_i0);
}

template<typename T>
T readOrDefault(YAML::Node &i_config, std::string i_key, T i_defaultValue) {
  try {
    if(!i_config[i_key]) return i_defaultValue;
    return i_config[i_key].as<T>();
  } catch(const YAML::Exception &ex) {
    std::cout << ex.what() << " for key '" << i_key << "', value '" << i_config[i_key].as<std::string>() << "', and type " << typeid(T).name() << std::endl;
    return i_defaultValue;
  }
}

bool fileExists(std::string i_fileName){
  struct stat buffer;   
  return stat(i_fileName.c_str(), &buffer) == 0;
}

std::string readFileAsString(std::string i_fileName){
  std::ifstream stream(i_fileName);
  std::stringstream buffer;
  buffer << stream.rdbuf();
  return buffer.str();
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
  
  char* l_configPath = i_argv[1];
  if(!fileExists(l_configPath)){
    std::cerr << "Configuration file could not be found!" << std::endl;
    return EXIT_FAILURE;
  }
  
  auto l_config = YAML::LoadFile(l_configPath);
  
  std::vector<tsunami_lab::io::Station> l_stations;
  tsunami_lab::setups::Setup* l_setup = nullptr;
  
  t_idx  l_nx, l_ny;
  t_real l_scale, l_cellSizeMeters;
  t_real l_cflFactor;
  
  t_idx  l_timeStepIndex;
  double l_simulationTime;
  
  // output parameters / time limits
  t_idx  l_maxTimesteps = readOrDefault<t_idx >(l_config, "maxSteps", std::numeric_limits<t_idx>::max());// max number of simulation timesteps
  t_real l_maxDuration     = readOrDefault<t_real>(l_config, "maxDuration", std::numeric_limits<t_real>::infinity());// max simulation time in seconds
  
  t_real l_outputPeriod = readOrDefault<t_real>(l_config, "outputPeriod", 1);// every n th second, a result file will be written
  t_idx  l_outputStepSize  = readOrDefault<t_idx >(l_config, "outputStepSize", 1);// n = every nth field is written to disk; saves storage space

  t_real l_gridOffsetX, l_gridOffsetY;
  
  // if there are more params from the console, we could add them to the config as well
  // we would need to include them in the hashing function for checkpoints then as well
  
  ///////////////////////////////////////
  // check whether a checkpoint exists //
  ///////////////////////////////////////
  size_t l_configHash = std::hash<std::string>{}(readFileAsString(l_configPath));
  std::string l_configPathStr(l_configPath);
  std::string l_configPathName = l_configPathStr.substr(l_configPathStr.find_last_of("/\\") + 1);
  std::string l_checkpointPathDefault = "checkpoint-" + l_configPathName + "." + std::to_string(l_configHash) + ".nc";
  std::string l_checkpointPath = readOrDefault<std::string>(l_config, "checkpointFile", l_checkpointPathDefault);
  t_idx l_checkpointingPeriod = readOrDefault<t_real>(l_config, "checkpointPeriod", 10 * 60);// default: create a checkpoint every 10 mins
  
  std::cout << "checkpoint file: " << l_checkpointPath << ", interval: " << l_checkpointingPeriod << "s" << std::endl;
  
  bool l_printStationComments = readOrDefault(l_config, "printStationComments", true);
  
  if(readOrDefault(l_config, "readCheckpoints", true) && fileExists(l_checkpointPath)){
    // h, hu, hv, b
    l_scale = 1;
    // todo create setup
    l_setup = tsunami_lab::io::NetCDF::loadCheckpoint(l_checkpointPath, l_nx, l_ny, l_cellSizeMeters, l_cflFactor, l_simulationTime, l_timeStepIndex, l_stations);
    // only required for the first frame
    l_gridOffsetX = 0;
    l_gridOffsetY = 0;
	// remove ghost cells
	l_nx -= 2;
	if(l_ny > 2) l_ny -= 2;
    if(l_setup == nullptr){
      std::cerr << "warn: checkpoint could not be loaded!, starting from zero!" << std::endl;
    } else {
      std::cout << "loaded checkpoint: " << l_simulationTime << "s, frameIndex: " << l_timeStepIndex << std::endl;
    }
  }
  
  // must stay in scope; I don't have a better solution currently; maybe the value could be moved into setup
  // a copy shouldn't be expensive compared to our simulation and writing the results to disk
  std::vector<t_real> l_bathymetry;
  std::vector<t_real> l_displacement;
  
  // setup is null if no checkpoint was found or it could not be loaded
  if(l_setup == nullptr){
    
    std::cout << "loading setup from config file" << std::endl;
    
    // many setups share similar setup information, so just request them all at once;
    // it should not matter for performance, as this all will be done within 1ms.
           l_nx              = readOrDefault<t_idx >(l_config, "nx",  1);
           l_ny              = readOrDefault<t_idx >(l_config, "ny",  1);
    t_real l_heightLeft      = readOrDefault<t_real>(l_config, "hl", 10);
    t_real l_heightRight     = readOrDefault<t_real>(l_config, "hr",  5);
    t_real l_impulseLeft     = readOrDefault<t_real>(l_config, "hul", 0);
    t_real l_impulseRight    = readOrDefault<t_real>(l_config, "hur", 0);
    t_real l_bathymetryLeft  = readOrDefault<t_real>(l_config, "bl",  0);
    t_real l_bathymetryRight = readOrDefault<t_real>(l_config, "br",  0);
    t_real l_damBathymetry   = readOrDefault<t_real>(l_config, "bathymetry", 0);// bathymetry for 1d and 2d dam
    // (dam bathymetry) caution: values above zero mean stopped fluid, because they mean a dry area, and the solver does not support wetting/drying.
    t_real l_splitPositionX  = readOrDefault<t_real>(l_config, "splitPositionX", l_nx * 0.5);
    t_real l_splitPositionY  = readOrDefault<t_real>(l_config, "splitPositionY", l_ny * 0.5);
    t_real l_damRadius       = readOrDefault<t_real>(l_config, "damRadius", l_nx / 4);
           l_cellSizeMeters  = readOrDefault<t_real>(l_config, "cellSize", 1);// size of a single cell in meters
           l_scale           = readOrDefault<t_real>(l_config, "scale", 1);// scales the setup for accuracy or performance
           l_cflFactor       = readOrDefault<t_real>(l_config, "cflFactor", 0.45);
           l_simulationTime  = 0;
           l_timeStepIndex   = 0;
    
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
  
    ////////////////
    // read setup //
    ////////////////
    std::string l_setupName = readOrDefault<std::string>(l_config, "setup", "Discontinuity1d");
    if(!l_config["setup"]) std::cout << "using default setup: " << l_setupName << std::endl;
  
    auto l_performanceTime0 = std::chrono::high_resolution_clock::now();
  
    // coordinates of the top left corner in the loaded tsunami data; used for stations
    // if not defined, this is the center of that corner
    l_gridOffsetX = l_cellSizeMeters * 0.5;
    l_gridOffsetY = l_cellSizeMeters * 0.5;
  
    if(
      l_setupName == "DamBreak" ||
      l_setupName == "DamBreak1d"
    ) {
      l_setup = new tsunami_lab::setups::DamBreak1d(l_heightLeft, l_heightRight, l_splitPositionX, l_damBathymetry);
    } else if(
      l_setupName == "DamBreakCircle" ||
      l_setupName == "DamBreak2d"
    ) {
      auto l_dambreak2d = new tsunami_lab::setups::DamBreak2d(l_heightLeft, l_heightRight, l_splitPositionX, l_splitPositionY, l_damRadius, l_damBathymetry);
      if(l_config["obstacleBathymetry"]) {
        // an obstacle should be defined
        l_dambreak2d->setObstacle(
          readOrDefault(l_config, "obstacleX0", 0),
          readOrDefault(l_config, "obstacleX1", 0), 
          readOrDefault(l_config, "obstacleY0", 0),
          readOrDefault(l_config, "obstacleY1", 0), 
          readOrDefault(l_config, "obstacleBathymetry", 0)
        );
      }
      l_setup = l_dambreak2d;
    } else if(
      l_setupName == "Discontinuity1d" ||
      l_setupName == "Discontinuity"
    ){
      l_setup = new tsunami_lab::setups::Discontinuity1d(l_heightLeft, l_heightRight, l_impulseLeft, l_impulseRight, l_bathymetryLeft, l_bathymetryRight, l_splitPositionX);
    } else if(
      l_setupName == "GMTTrack" || 
      l_setupName == "Tsunami"
    ) {
      if(l_config["setupFile"]){
      
        // load the data from the csv file
        auto l_fileName = l_config["setupFile"].as<std::string>();
        if(!fileExists(l_fileName)){
          std::cerr << "could not find setup '" << l_fileName << "'" << std::endl;
          return EXIT_FAILURE;
        }
      
        auto l_loadedData = tsunami_lab::io::Csv::read(l_fileName);
      
        // here we probably copy; this potentially could be optimized
        auto l_xs = tsunami_lab::io::Csv::findColumn(l_loadedData, "track_location");
        l_bathymetry = tsunami_lab::io::Csv::findColumn(l_loadedData, "height");
      
        if(l_xs.size() < 1){
          std::cerr << "did not find position data in track file" << std::endl;
          return EXIT_FAILURE;
        } else if(l_bathymetry.size() < 1){
          std::cerr << "did not find bathymetry data in track file" << std::endl;
          return EXIT_FAILURE;
        }
      
        l_cellSizeMeters = (l_xs.back() - l_xs.front()) / (l_xs.size() - 1) / l_scale;
      
        l_nx = l_xs.size() * l_scale - 2;// 2 = ghost cells
        l_ny = 1;// it's just a 1d simulation
      
        t_real l_displacementStart  = readOrDefault<t_real>(l_config, "displacementStart", 175000) / l_cellSizeMeters;
        t_real l_displacementEnd    = readOrDefault<t_real>(l_config, "displacementEnd",   250000) / l_cellSizeMeters;
        t_real l_displacementHeight = readOrDefault<t_real>(l_config, "displacement", 10);
      
        // construct setup
        l_setup = new tsunami_lab::setups::TsunamiEvent1d(l_bathymetry.data(), l_bathymetry.size(), 1.0 / l_scale, l_displacementStart, l_displacementEnd, l_displacementHeight);
      
      } else {
        std::cerr << "missing parameter 'setupFile' for setup type Tsunami1d/GMTTrack" << std::endl;
        return EXIT_FAILURE;
      }
    } else if(
      l_setupName == "Tsunami2d" ||
      l_setupName == "NetCDF"
    ){
    
      if(!l_config["bathymetryFile"]){
        std::cerr << "missing parameter 'bathymetryFile' for setup type Tsunami2d/NetCDF" << std::endl;
        return EXIT_FAILURE;
      }
      if(!l_config["displacementFile"]){
        std::cerr << "missing parameter 'displacementFile' for setup type Tsunami2d/NetCDF" << std::endl;
        return EXIT_FAILURE;
      }
    
      auto l_bathymetryFileName = l_config["bathymetryFile"].as<std::string>();
      auto l_displacementFileName = l_config["displacementFile"].as<std::string>();
      if(l_bathymetryFileName == l_displacementFileName){// as long as the property is not yet customizable
        std::cout << "---------------------------------------------------------" << std::endl;
        std::cout << "warning: bathymetry and displacement share the same data!" << std::endl;
        std::cout << "---------------------------------------------------------" << std::endl;
      }
      if(!fileExists(l_bathymetryFileName)){
        std::cerr << "could not find bathymetry file '" << l_bathymetryFileName << "'" << std::endl;
        return EXIT_FAILURE;
      } else if(!fileExists(l_displacementFileName)){
        std::cerr << "could not find displacement file '" << l_displacementFileName << "'" << std::endl;
        return EXIT_FAILURE;
      }
    
      // load data
      t_idx l_nx2, l_ny2;
      t_real l_cellSizeMeters2;
      t_real l_tmp;
      tsunami_lab::io::NetCDF::load2dArray(l_bathymetryFileName, "z", l_nx, l_ny, l_cellSizeMeters2, l_gridOffsetX, l_gridOffsetY, l_bathymetry);
      tsunami_lab::io::NetCDF::load2dArray(l_displacementFileName, "z", l_nx2, l_ny2, l_cellSizeMeters2, l_tmp, l_tmp, l_displacement);
      if(l_cellSizeMeters == 1.0){
        l_cellSizeMeters = l_cellSizeMeters2 / l_scale;// cell size depends on data & applied scale
      } else std::cout << "used cell size override from config. Cell size from file: " << l_cellSizeMeters2 << std::endl;
      t_real l_sideRatio = (l_nx2 * l_ny) / (t_real) (l_ny2 * l_nx); // ideally 1
      if(l_sideRatio < 0.99 || l_sideRatio > 1.01){
        std::cerr << "warning: aspect ratio from bathymetry and displacement are different!" << std::endl;
      }
      if(l_nx < 2 || l_ny < 2 || l_nx2 < 2 || l_ny2 < 2){
        std::cerr << "data must have at least 2 x 2 fields, read " << l_nx << " x " << l_ny << ", " << l_nx2 << " x " << l_ny2 << std::endl;
        return EXIT_FAILURE;
      }
      // create setup
	  t_real l_scaleBath = 1.0 / (l_scale * l_scale);
	  t_real l_scaleDisp = 1.0 / (l_scale * l_scale) * std::sqrt((l_nx2 * l_ny2)/(t_real)(l_nx * l_ny));
      l_setup = new tsunami_lab::setups::TsunamiEvent2d(
        l_bathymetry.data(), l_nx, l_ny, l_nx, l_scaleBath,
        l_displacement.data(), l_nx2, l_ny2, l_nx2, l_scaleDisp
      );
      // scale it up
      l_nx = (l_nx-2) * l_scale;// 2 for ghost cells; those cannot be scaled, and are re-added by the WavePropagation class
      l_ny = (l_ny-2) * l_scale;
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
    } else if(
      l_setupName == "ArtificialTsunami2d"
    ){
      l_setup = new tsunami_lab::setups::ArtificialTsunami2d();
    } else {
      std::cerr << "unknown/invalid setup type \"" << l_setupName << "\"" << std::endl;
      return EXIT_FAILURE;
    }
  
    ///////////////////
    // read stations //
    ///////////////////
    if(l_config["stations"]) {
      auto   l_stationData = l_config["stations"].as<std::vector<YAML::Node>>();
      t_real l_delayBetweenRecords = readOrDefault(l_config, "delayBetweenRecords", 1);
      for(t_idx i=0;i<l_stationData.size();i++){
        
        auto l_station = l_stationData[i];
        std::string l_name = l_station["name"].as<std::string>();
        
        int64_t l_x, l_y;
        if(l_station["x"]){
          
          l_x = l_station["x"].as<int64_t>();
          l_y = l_station["y"].as<int64_t>();
          
        } else if(l_station["gridX"]){
          
          t_real l_x0 = l_station["gridX"].as<double>();
          t_real l_y0 = l_station["gridY"].as<double>();
          
          l_x = (l_x0 - l_gridOffsetX) / l_cellSizeMeters;
          l_y = (l_y0 - l_gridOffsetY) / l_cellSizeMeters;
          
          std::cout << "Station " << i << ", '" << l_name << "' with input grid coordinates (" << l_x0 << ", " << l_y0 << ") was mapped to simulation grid coordinates (" << l_x << ", " << l_y << ")" << std::endl;
          
        } else {
          
          std::cerr << "Missing location for station '" << l_name << "'" << std::endl;
          return EXIT_FAILURE;
          
        }
        
        if(l_x >= 0 && (t_idx) l_x < l_nx && l_y >= 0 && (t_idx) l_y < l_ny){
          t_real l_stationBathymetry = l_setup->getBathymetry(l_x * l_scale, l_y * l_scale);
          if(l_stationBathymetry <= 0){
            tsunami_lab::io::Station l_station1((t_idx) l_x, (t_idx) l_y, l_name, l_delayBetweenRecords);
            l_stations.push_back(l_station1);
            std::cout << "Station " << i << ", '" << l_name << "' was placed on bathymetry " << l_stationBathymetry << std::endl;
          } else {
            std::cerr << "Station " << i << ", '" << l_name << "' is above sea level: " << l_stationBathymetry << std::endl;
            return EXIT_FAILURE;
          }
        } else {
          std::cerr << "Station " << i << ", '" << l_name << "' is out of bounds: " << l_x << "," << l_y << " !in 0.." << l_nx << ",0.." << l_ny << std::endl;
          std::cerr << "This invalid entry was found in the file '" << i_argv[1] << "'" << std::endl;
          return EXIT_FAILURE;
        }
      }
    }
  
    auto l_performanceTime1 = std::chrono::high_resolution_clock::now();
    double l_dur0 = std::chrono::duration<double>(l_performanceTime1-l_performanceTime0).count();
    std::cout << "used " << l_dur0 << "s to load setup" << std::endl;
  
    // run simulation
    std::cout << "runtime configuration" << std::endl;
    std::cout << "  configuration file:             " << i_argv[1] << std::endl;
    std::cout << "  setup type:                     " << l_setupName << std::endl;
    std::cout << "  max timesteps:                  " << l_maxTimesteps << std::endl;
    std::cout << "  max duration:                   " << l_maxDuration << std::endl;
    std::cout << "  number of cells in x-direction: " << l_nx << std::endl;
    std::cout << "  number of cells in y-direction: " << l_ny << std::endl;
    std::cout << "  cell size (meters):             " << l_cellSizeMeters << std::endl;
    std::cout << "  number of stations:             " << l_stations.size() << std::endl;
    std::cout << "  split position x:               " << l_splitPositionX << std::endl;
    std::cout << "  split position y:               " << l_splitPositionY << std::endl;
    std::cout << "  applied scale:                  " << l_scale << std::endl;
    std::cout << "  cfl factor:                     " << l_cflFactor << std::endl;
    std::cout << "  OpenMP threads:                 " << omp_get_max_threads() << std::endl;
  }
  
  auto l_performanceTime1 = std::chrono::high_resolution_clock::now();
  
  // construct solver
  tsunami_lab::patches::WavePropagation* l_waveProp;
  if(l_ny <= 1){
    l_waveProp = new tsunami_lab::patches::WavePropagation1d(l_nx, l_setup, l_scale);
  } else {
    auto l_waveProp2 = new tsunami_lab::patches::WavePropagation2d(l_nx, l_ny, l_setup, l_scale, l_scale);
	l_waveProp2->setCflFactor(l_cflFactor);
	l_waveProp = l_waveProp2;
  }
  
  // no longer needed
  // l_bathymetry.resize(0);
  
  // set up print control
  t_idx  l_nOut = 0;
  t_real l_timestep;

  std::cout << "entering time loop" << std::endl;
  
  std::string l_netCdfPath = readOrDefault<std::string>(l_config, "outputFile", "solution.nc");
  int    l_deflateLevel = readOrDefault<int>(l_config, "outputCompression", 5);
  bool   l_exportCSV = readOrDefault(l_config, "exportCSV", false);
  double l_debugPrintPerformanceInterval = readOrDefault<double>(l_config, "debugPrintPerformanceInterval", 1.0);
  
  auto l_performanceTimeDebug0 = std::chrono::high_resolution_clock::now();
  auto l_checkpointingTime0 = l_performanceTimeDebug0;
  
  // iterate over time
  t_idx l_lastOutputIndex = -1;
  t_idx l_timeStepIndexPerf = l_timeStepIndex;
  for(; l_timeStepIndex < l_maxTimesteps && l_simulationTime < l_maxDuration; l_timeStepIndex++ ){
    
    auto l_stepTime = std::chrono::high_resolution_clock::now();
    double l_durI = std::chrono::duration<double>(l_stepTime-l_performanceTimeDebug0).count();
    if(l_durI >= l_debugPrintPerformanceInterval){
      double l_stepsPerSecond = (l_timeStepIndex - l_timeStepIndexPerf) / l_durI;
      std::cout << "  step: " << l_timeStepIndex << ", simulation time: " << l_simulationTime << ", steps per second: " << l_stepsPerSecond << std::endl;
      l_performanceTimeDebug0 = l_stepTime;
      l_timeStepIndexPerf = l_timeStepIndex;
    }
    
    double l_durI2 = std::chrono::duration<double>(l_stepTime-l_checkpointingTime0).count();
    if(l_durI2 >= l_checkpointingPeriod){
      // create a new checkpoint
      std::cout << "  saving checkpoint" << std::endl;
      tsunami_lab::io::NetCDF::storeCheckpoint(l_checkpointPath, l_nx, l_ny, l_cellSizeMeters, l_cflFactor, l_simulationTime, l_timeStepIndex, l_stations, l_waveProp);
      std::cout << "  finished saving checkpoint" << std::endl;
      l_checkpointingTime0 = std::chrono::high_resolution_clock::now();// reset the timer for the next checkpoint
    }
    
    // index, which frame we'd need to print theoretically
    // if there are more frames requested than simulated, we just skip some
    t_idx l_outputIndex = (t_idx) (l_simulationTime / l_outputPeriod);
    if(l_timeStepIndex == 0 || l_lastOutputIndex != l_outputIndex) {
      l_lastOutputIndex = l_outputIndex;
      
      std::cout << "  simulation time: " << l_simulationTime << ", #time steps: "<< l_timeStepIndex << std::endl;
      
      if(l_exportCSV){
        
        std::string l_path = "solution_" + std::to_string(l_nOut) + ".csv";
        std::cout << "  writing wave field to " << l_path << std::endl;
        
        std::ofstream l_file(l_path, std::ios::out);
        tsunami_lab::io::Csv::write(l_cellSizeMeters, l_nx, l_ny, l_outputStepSize, l_waveProp->getStride(), l_waveProp->getHeight(), l_waveProp->getMomentumX(), l_waveProp->getMomentumY(), l_waveProp->getBathymetry(), l_file);
        l_file.close();
      
        l_nOut++;
        
      } else {
        if(tsunami_lab::io::NetCDF::appendTimeframe( l_cellSizeMeters, l_nx, l_ny, l_gridOffsetX, l_gridOffsetY, l_outputStepSize, l_waveProp->getStride(), l_waveProp->getHeight(), l_waveProp->getMomentumX(), l_waveProp->getMomentumY(), l_waveProp->getBathymetry(), l_setup, l_simulationTime, l_deflateLevel, l_netCdfPath)) return EXIT_FAILURE;
      }
	  
	  // only needed for file export
	  // and we may need the memory
	  // l_displacement.resize(0);
      
    }
    
    // update recording stations, if there are any
    if(!l_stations.empty() && l_stations[0].needsUpdate(l_simulationTime)) {
      for(auto &l_station : l_stations) {
        l_station.recordState(*l_waveProp, l_simulationTime);
      }
    }

    l_waveProp->setGhostOutflow();
    l_timestep = l_waveProp->computeMaxTimestep(l_cellSizeMeters);
    if(!std::isfinite(l_timestep)){
      std::cerr << "  there no longer is any valid fluid in the simulation! Stopping." << std::endl;
      break;// NaN or Infinite timestep -> illegal -> stop simulation
    }
    
    t_real l_scaling = l_timestep / l_cellSizeMeters;
    l_waveProp->timeStep(l_scaling);

    l_simulationTime += l_timestep;
  }
  
  std::cout << "finished time loop" << std::endl;
  
  // print last state
  std::cout << "  simulation end time: " << l_simulationTime << ", #total time steps: "<< l_timeStepIndex << std::endl;
  
  auto l_performanceTimeN = std::chrono::high_resolution_clock::now();
  double l_durN = std::chrono::duration<double>(l_performanceTimeN-l_performanceTime1).count();
  double l_stepsPerSecond = l_timeStepIndex / l_durN;
  std::cout << "average steps per second: " << l_stepsPerSecond << ", total simulation time: " << l_durN << std::endl;
  
  if(l_exportCSV){
      
    std::string l_path = "solution_" + std::to_string(l_nOut) + ".csv";
    std::cout << "  writing wave field to " << l_path << std::endl;
    
    std::ofstream l_file(l_path, std::ios::out);
    tsunami_lab::io::Csv::write( l_cellSizeMeters, l_nx, l_ny, l_outputStepSize, l_waveProp->getStride(), l_waveProp->getHeight(), l_waveProp->getMomentumX(), l_waveProp->getMomentumY(), l_waveProp->getBathymetry(), l_file );
    l_file.close();
    
  } else {
    if(tsunami_lab::io::NetCDF::appendTimeframe( l_cellSizeMeters, l_nx, l_ny, l_gridOffsetX, l_gridOffsetY, l_outputStepSize, l_waveProp->getStride(), l_waveProp->getHeight(), l_waveProp->getMomentumX(), l_waveProp->getMomentumY(), l_waveProp->getBathymetry(), l_setup, l_simulationTime, l_deflateLevel, l_netCdfPath)) return EXIT_FAILURE;
  }
  
  // todo init files once, then only append the measurements
  for(auto &l_station : l_stations){
    l_station.write(l_printStationComments);
  }
  
  std::cout << "finished writing last state" << std::endl;
  
  delete l_waveProp;
  
  delete l_setup;
  
  return EXIT_SUCCESS;
}
