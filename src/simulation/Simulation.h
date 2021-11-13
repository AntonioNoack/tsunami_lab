/**
 * @author Antonio Noack
 *
 * @section DESCRIPTION
 * Simulation functions.
 **/
#ifndef TSUNAMI_LAB_SIMULATION_SIMULATION_H
#define TSUNAMI_LAB_SIMULATION_SIMULATION_H

#include <vector>

#include "../constants.h"
#include "../setups/Setup.h"
#include "../io/Station.h"

namespace tsunami_lab {
  namespace simulation {
    class Simulation;
  }
}

class tsunami_lab::simulation::Simulation {

  public:
    /**
     * Runs a simulation created from setup, and writes the result into csv files.
     *
     * @param i_nx number of fields on the x axis.
     * @param i_ny number of fields on the y axis.
     * @param i_setup initial state of simulation.
     * @param i_setupScale scale for initial state.
     * @param i_cellSizeMeters size of cells in meters.
     * @param i_maxTimesteps how many timesteps will be computed at maximum.
     * @param i_maxDuration how many seconds will be computed at maximum.
     * @param i_outputStepSize every step-th cell will be printed to save storage space.
     * @param i_numOutputSteps how many time steps will be written to disk.
	 * @param i_stations recording wave-measuring stations
     **/
    static void run( tsunami_lab::t_idx                      i_nx,
                     tsunami_lab::t_idx                      i_ny,
                     tsunami_lab::setups::Setup            & i_setup,
                     tsunami_lab::t_real                     i_setupScale,
                     tsunami_lab::t_real                     i_cellSizeMeters,
                     tsunami_lab::t_idx                      i_maxTimesteps,
					 tsunami_lab::t_real                     i_maxDuration,
                     tsunami_lab::t_idx                      i_outputStepSize,
                     tsunami_lab::t_idx                      i_numOutputSteps,
                     std::vector<tsunami_lab::io::Station> & i_stations );

};

#endif