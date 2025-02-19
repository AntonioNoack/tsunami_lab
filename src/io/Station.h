/**
 * @author Antonio Noack
 *
 * @section DESCRIPTION
 * Wave-recording station, which watches height and momentum.
 **/
#ifndef TSUNAMI_LAB_IO_STATION
#define TSUNAMI_LAB_IO_STATION

#include "../constants.h"
#include "../patches/WavePropagation.h"

#include <cstring>
#include <iostream>
#include <vector>

namespace tsunami_lab {
  namespace io {
    class Station;
    struct RecordedValue {
      t_real time;
      t_real height;
      t_real momentumX;
      t_real momentumY;
    };
  }
}

class tsunami_lab::io::Station {
  private:
    t_idx m_positionX;
    t_idx m_positionY;
    
    //! station name, e.g. for saving the recorded data
    std::string m_name;
    
    t_real m_nextRecordTime = 0;
    t_real m_delayBetweenRecords = 0;
    
    //! all recorded values
    std::vector<RecordedValue> m_records;
    
  public:
  
    /**
     * Constructor.
     *
     * @param i_positionX position in x-direction on the simulated grid.
     * @param i_positionX position in y-direction on the simulated grid.
     * @param i_name station name.
     * @param i_delayBetweenRecords minimum delay between records in seconds.
     **/
    Station( t_idx i_positionX, t_idx i_positionY, std::string i_name, t_real i_delayBetweenRecords):
      m_positionX(i_positionX), m_positionY(i_positionY), m_name(i_name), m_delayBetweenRecords(i_delayBetweenRecords) {}
    
    /**
     * Returns whether a new record should be written at the given point in time.
     *
     * @param i_time queried time in seconds.
     * @return whether a new record should be added.
     **/
    bool needsUpdate( t_real l_time );
    
    /**
     * Records the current wave state to memory on the station.
     * Then sets the time for the next record.
     *
     * @param i_waveProp running simulation.
     * @param i_time current time in seconds.
     **/
    void recordState( tsunami_lab::patches::WavePropagation & i_waveProp, t_real i_time );
    
    /**
     * Records the current wave state to memory on the station.
     * Then sets the time for the next record.
     *
     * @param i_time current time in seconds.
     * @param i_height water height in meters.
     * @param i_momentumX water momentum in x direction.
     * @param i_momentumY water momentum in y direction.
     **/
    void recordState( t_real i_time, t_real i_height, t_real i_momentumX, t_real i_momentumY );
    
    /**
     * Writes all recorded data to a csv file.
     *
     * @param i_withComment whether the comment should be printed.
     **/
    void write( bool i_withComment );
    
	// just a few getters; is there a better way in C++? (besided making the properties public)
    /** self explanatory */
    std::string getName(){
      return m_name;
    }
	
	/** self explanatory */
    void getPosition(t_idx &o_x, t_idx &o_y){
      o_x = m_positionX;
	  o_y = m_positionY;
    }
    
    /** self explanatory */
    std::vector<RecordedValue>& getRecords(){
      return m_records;
    }
    
    /** self explanatory */
    t_real getNextRecordTime(){
      return m_nextRecordTime;
    }
    
    /** self explanatory */
    t_real getDelayBetweenRecords(){
      return m_delayBetweenRecords;
    }
    
};

#endif