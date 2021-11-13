/**
 * @author Antonio Noack
 * 
 * @section DESCRIPTION
 * Unit tests for the wave-recording stations.
 **/
#include <catch2/catch.hpp>
#include "../constants.h"
#define private public
#include "Station.h"
#undef public

#define t_idx  tsunami_lab::t_idx
#define t_real tsunami_lab::t_real

TEST_CASE( "Test the stations.", "[Station]" ) {
  tsunami_lab::io::Station station( 10, 5, "test_station", 1.5 );
  station.m_nextRecordTime = 10 + 1.5;
  REQUIRE( station.needsUpdate(10 + 1.4) == false );
  REQUIRE( station.needsUpdate(10 + 1.6) == true );
}

#undef t_idx
#undef t_real
