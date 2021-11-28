/**
 * @author Antonio Noack
 * 
 * @section DESCRIPTION
 * Wave-recording station, which watches height and momentum.
 **/
#include "Station.h"
#include <sstream>
#include <fstream> // file streams

bool tsunami_lab::io::Station::needsUpdate( t_real i_time ) {
  return m_nextRecordTime <= i_time;
}

void tsunami_lab::io::Station::recordState( tsunami_lab::patches::WavePropagation &i_waveProp, t_real i_time ) {
  t_idx l_index = m_positionX + m_positionY * i_waveProp.getStride();
  auto l_hPtr = i_waveProp.getHeight();// somehow in the sanitizer this has become null... how?
  t_real l_h  = l_hPtr ? l_hPtr[l_index] : 0;
  auto l_huPtr = i_waveProp.getMomentumX();
  auto l_hvPtr = i_waveProp.getMomentumY();
  t_real l_hu = l_huPtr ? l_huPtr[l_index] : 0;
  t_real l_hv = l_hvPtr ? l_hvPtr[l_index] : 0;
  m_records.push_back( { i_time, l_h, l_hu, l_hv } );
  m_nextRecordTime = i_time + m_delayBetweenRecords;
}

void tsunami_lab::io::Station::write( bool i_withComment ) {
  
  std::string l_fileName = "station_" + m_name + ".csv";
  
  std::ofstream l_stream;
  l_stream.open( l_fileName );
  if(!l_stream.is_open()) {
    std::cerr << "Error opening file " << l_fileName << std::endl;
    return;
  }
  
  // print information about the station in some comments in the header
  if(i_withComment){
    l_stream << "# Station " << m_name << "\n";
    l_stream << "# Location (Grid) " << m_positionX << "," << m_positionY << "\n";
  }
  
  // write the CSV table header
  l_stream << "time,height,momentumX,momentumY\n";

  // iterate over all records
  for(auto &l_record : m_records) {
    l_stream << l_record.time << "," << l_record.height << "," << l_record.momentumX << "," << l_record.momentumY << "\n";
  }
  l_stream << std::flush;
  l_stream.close();
  
}