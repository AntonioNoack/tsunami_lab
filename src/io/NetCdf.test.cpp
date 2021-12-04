/**
 * @author Antonio Noack
 * 
 * @section DESCRIPTION
 * Unit tests for the NetCDF-interface.
 **/
#include <catch2/catch.hpp>
#include "../constants.h"
#include "../patches/WavePropagation1d.h"

#include <sstream>
#include <iostream>
#include <cstdio>
#include <netcdf.h>

#define private public
#include "NetCdf.h"
#undef public

#define t_idx  tsunami_lab::t_idx
#define t_real tsunami_lab::t_real

#ifdef check
#error "already defined check()"
#endif

#define check(error) {\
  l_err = error;\
  if(l_err != NC_NOERR){\
    std::cerr << "NetCDF-Error occurred: " << nc_strerror(l_err) << " (Code " << l_err << "), line " << __LINE__ << std::endl;\
    REQUIRE(l_err == NC_NOERR);\
  }\
}

TEST_CASE( "Test the NetCDF-writer to append data.", "[NetCDF][AppendData]" ) {

  int l_deflateLevel = 5;
  
  // define a simple example
  t_real l_h [16]  = { 0,  1,  2,  3,
                       4,  5,  6,  7,
                       8,  9, 10, 11,
                      12, 13, 14, 15 };
  t_real l_hu[16] = { 15, 14, 13, 12,
                      11, 10,  9,  8,
                       7,  6,  5,  4,
                       3,  2,  1,  0 };
  t_real l_hv[16] = {  0,  4,  8, 12,
                       1,  5,  9, 13,
                       2,  6, 10, 14,
                       3,  7, 11, 15 };

  t_real l_b [16] = {  0, -1, -2, -6,
                       0, -1, -3, -4,
                       0, -2, -3, -3,
                       1,  0, -1, -3 };

  t_real l_cellSizeMeters = 10;
  t_real l_gridOffsetX = l_cellSizeMeters * 0.5;
  t_real l_gridOffsetY = l_cellSizeMeters * 0.5;
  t_idx l_nx = 4, l_ny = 4, l_nt = 3;
  t_idx l_step = 1, l_stride = l_nx;
  
  std::string l_fileName = "tmp.nc";
  
  for(t_idx i=0;i<l_nt;i++){
    tsunami_lab::io::NetCDF::appendTimeframe(l_cellSizeMeters, l_nx, l_ny, l_gridOffsetX, l_gridOffsetY, l_step, l_stride, l_h, l_hu, l_hv, l_b, nullptr, (t_real) (i * 0.5), i, l_deflateLevel, l_fileName);
  }
  
  int l_err;
  
  int l_handle; // file handle
  
  check(nc_open(l_fileName.c_str(), NC_NOWRITE, &l_handle));
  
  int l_xDimId, l_yDimId, l_tDimId;
  check(nc_inq_dimid(l_handle, "x",    &l_xDimId));
  check(nc_inq_dimid(l_handle, "y",    &l_yDimId));
  check(nc_inq_dimid(l_handle, "time", &l_tDimId));
  
  t_idx l_nx2, l_ny2, l_nt2;
  
  char l_tmpDimName[NC_MAX_NAME+1];// name, ignored
  check(nc_inq_dim(l_handle, l_xDimId, l_tmpDimName, &l_nx2));
  check(nc_inq_dim(l_handle, l_yDimId, l_tmpDimName, &l_ny2));
  check(nc_inq_dim(l_handle, l_tDimId, l_tmpDimName, &l_nt2));
  
  REQUIRE(l_nx2 == l_nx);
  REQUIRE(l_ny2 == l_ny);
  REQUIRE(l_nt2 == l_nt);
  
  check(nc_close(l_handle));
  
}

TEST_CASE( "Test the NetCDF-reading functionality for a 2d field named z.", "[NetCDF][Read2d]" ) {
  // for the baseline see the file ../../data/netcdf-test.nc
  t_idx l_sizeX = 0, l_sizeY = 0;
  t_real l_cellSizeMeters = 0;
  std::vector<t_real> l_data;
  t_real l_tmp;
  tsunami_lab::io::NetCDF::load2dArray("data/netcdf-test.nc", "z", l_sizeX, l_sizeY, l_cellSizeMeters, l_tmp, l_tmp, l_data);
  REQUIRE(l_sizeX == 4);
  REQUIRE(l_sizeY == 3);
  REQUIRE(l_cellSizeMeters == Approx(2.0));
  REQUIRE(l_data.size() == 12);
  for(int i=0;i<12;i++) REQUIRE(l_data[i] == i);
}

TEST_CASE( "Test sampling down an image for coarse output", "[NetCDF][downsample]" ) {
  const t_idx l_sizeIn = 5, l_sizeOut = 2, l_step = 3;
  const t_real l_dataIn[25] = {
      5, 5, 5,  3, 4,
      6, 4, 5,  2, 2,
      4, 6, 5,  4, 3,
      
      3, 2, 1,  7, 8,
      3, 1, 2,  8, 9
  };
  t_real l_dataOut[5] = { 0, 0, 0, 0, 19 };
  tsunami_lab::io::NetCDF::downsample(l_sizeIn, l_sizeIn, l_sizeIn, l_dataIn, l_sizeOut, l_sizeOut, l_sizeOut, l_dataOut, l_step);
  REQUIRE(l_dataOut[0] == 5);
  REQUIRE(l_dataOut[1] == 3);
  REQUIRE(l_dataOut[2] == 2);
  REQUIRE(l_dataOut[3] == 8);
  REQUIRE(l_dataOut[4] == 19);
}

TEST_CASE( "Test checkpointing", "[NetCDF][Checkpointing]" ) {
  
  /*
   * Test case from WavePropagation1d:
   *
   *   Single dam break problem between cell 49 and 50.
   *     left | right
   *       10 | 8
   *        0 | 0
   *
   *   Elsewhere steady state.
   *
   * The net-updates at the respective edge are given as
   * (see derivation in Roe solver):
   *    left          | right
   *      9.394671362 | -9.394671362
   *    -88.25985     | -88.25985
   */

  // define variables before & after serialization
  t_idx l_nx = 100, l_ny = 1, l_timeStepIndex = 21;
  t_real l_cellSizeMeters = 0.345, l_cflFactor = 0.25;
  double l_simulationTime = 0.123;
  
  // sample stations
  std::vector<tsunami_lab::io::Station> l_stations;
  l_stations.push_back(tsunami_lab::io::Station(2, 5, "Bikini Bottom", 0.17));
  l_stations[0].recordState(0,1,2,3);
  l_stations[0].recordState(4,5,6,7);
  l_stations[0].recordState(8,9,1,2);
  
  l_stations.push_back(tsunami_lab::io::Station(9, 13, "Atlantis", 0.17));
  l_stations[1].recordState(2,1,1,2);
  l_stations[1].recordState(7,8,8,7);
  l_stations[1].recordState(1,2,2,1);
  
  t_idx l_nx2, l_ny2, l_timeStepIndex2;
  t_real l_cellSizeMeters2, l_cflFactor2;
  double l_simulationTime2;
  std::vector<tsunami_lab::io::Station> l_stations2;
  
  
  // construct solver and setup a dam break problem
  tsunami_lab::patches::WavePropagation1d l_waveProp(l_nx);
  
  for( std::size_t l_ce = 0; l_ce < 50; l_ce++ ) {
    l_waveProp.setHeight( l_ce, 0, 10 );
    l_waveProp.setMomentumX( l_ce, 0, 0 );
  }
  for( std::size_t l_ce = 50; l_ce < 100; l_ce++ ) {
    l_waveProp.setHeight( l_ce, 0, 8 );
    l_waveProp.setMomentumX( l_ce, 0, 0 );
  }

  // set outflow boundary condition
  l_waveProp.setGhostOutflow();

  // perform a time step
  l_waveProp.timeStep( 0.1 );
  
  std::string l_fileName = "tmp-cp.nc";
  std::remove(l_fileName.c_str()); // delete just in case it exists
  
  // save the checkpoint
  auto l_error = tsunami_lab::io::NetCDF::storeCheckpoint( l_fileName, l_nx, l_ny, l_cellSizeMeters,
                                                           l_cflFactor, l_simulationTime, l_timeStepIndex,
                                                           l_stations, &l_waveProp );
  REQUIRE(l_error == 0);
  
  // load the checkpoint
  auto l_setup = tsunami_lab::io::NetCDF::loadCheckpoint( l_fileName, l_nx2, l_ny2, l_cellSizeMeters2,
														  l_cflFactor2, l_simulationTime2, l_timeStepIndex2,
                                                          l_stations2 );
														  
  REQUIRE(l_nx2 == l_nx);
  REQUIRE(l_ny2 == l_ny);
  REQUIRE(l_cellSizeMeters2 == l_cellSizeMeters);
  REQUIRE(l_cflFactor2 == l_cflFactor);
  REQUIRE(l_simulationTime2 == l_simulationTime);
  REQUIRE(l_timeStepIndex2 == l_timeStepIndex);
  
  // test the stations (all properties at least once)
  REQUIRE(l_stations2.size() == 2);
  REQUIRE(l_stations2[0].getName() == "Bikini Bottom");
  t_idx px,py;l_stations2[0].getPosition(px,py);
  REQUIRE(px == 2);
  REQUIRE(py == 5);
  REQUIRE(l_stations2[0].getDelayBetweenRecords() == Approx(0.17));
  REQUIRE(l_stations2[0].getRecords().size() == 3);
  REQUIRE(l_stations2[0].getRecords()[2].time == 8);
  REQUIRE(l_stations2[0].getRecords()[2].height == 9);
  REQUIRE(l_stations2[0].getRecords()[2].momentumX == 1);
  REQUIRE(l_stations2[0].getRecords()[2].momentumY == 2);
  REQUIRE(l_stations2[1].getName() == "Atlantis");
  REQUIRE(l_stations2[1].getRecords().size() == 3);
  REQUIRE(l_stations2[1].getRecords()[2].time == 1);
  REQUIRE(l_stations2[1].getRecords()[2].height == 2);
  REQUIRE(l_stations2[1].getRecords()[2].momentumX == 2);
  REQUIRE(l_stations2[1].getRecords()[2].momentumY == 1);
  
  // create the new scenario from the loaded data
  tsunami_lab::patches::WavePropagation1d l_waveProp2(l_nx);
  l_waveProp2.initWithSetup(l_setup, 1.0);

  // steady state
  for( std::size_t l_ce = 0; l_ce < 49; l_ce++ ) {
    REQUIRE( l_waveProp2.getHeight()   [l_ce] == Approx(10) );
    REQUIRE( l_waveProp2.getMomentumX()[l_ce] == Approx(0) );
  }

  // dam-break
  // margin added for support of different gravity constants (9.81 vs ~9.8066)
  REQUIRE( l_waveProp2.getHeight()   [49] == Approx(10 - 0.1 * 9.394671362).margin(0.001) );
  REQUIRE( l_waveProp2.getMomentumX()[49] == Approx( 0 + 0.1 * 88.25985).margin(0.01) );

  REQUIRE( l_waveProp2.getHeight()   [50] == Approx(8 + 0.1 * 9.394671362).margin(0.001) );
  REQUIRE( l_waveProp2.getMomentumX()[50] == Approx(0 + 0.1 * 88.25985).margin(0.01) );

  // steady state
  for( std::size_t l_ce = 51; l_ce < 100; l_ce++ ) {
    REQUIRE( l_waveProp2.getHeight()   [l_ce] == Approx(8) );
    REQUIRE( l_waveProp2.getMomentumX()[l_ce] == Approx(0) );
  }
  
  delete l_setup;
  
}
