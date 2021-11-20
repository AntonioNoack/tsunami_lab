/**
 * @author Antonio Noack
 * 
 * @section DESCRIPTION
 * Unit tests for the NetCDF-interface.
 **/
#include <catch2/catch.hpp>
#include "../constants.h"
#include <sstream>
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
    printf("NetCDF-Error occurred: %s (Code %d), line %d\n", nc_strerror(l_err), l_err, __LINE__);\
    REQUIRE(l_err == NC_NOERR);\
  }\
}

TEST_CASE( "Test the NetCDF-writer to append data.", "[NetCDF][AppendData]" ) {

  
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
  t_idx l_nx = 4, l_ny = 4, l_nt = 3;
  t_idx l_step = 1, l_stride = l_nx;
  
  std::string l_fileName = "tmp.nc";
  
  for(t_idx i=0;i<l_nt;i++){
    tsunami_lab::io::NetCDF::appendTimeframe(l_cellSizeMeters, l_nx, l_ny, l_step, l_stride, l_h, l_hu, l_hv, l_b, nullptr, (t_real) (i * 0.5), i, l_fileName);
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
  tsunami_lab::io::NetCDF::load2dArray("data/netcdf-test.nc", "z", l_sizeX, l_sizeY, l_cellSizeMeters, l_data);
  REQUIRE(l_sizeX == 4);
  REQUIRE(l_sizeY == 3);
  REQUIRE(l_cellSizeMeters == Approx(2.0));
  REQUIRE(l_data.size() == 12);
  for(int i=0;i<12;i++) REQUIRE(l_data[i] == i);
}
// we could create a test, where we read and write data; issue: we only read in z, and output height/bathymetry/momentum currently

#undef t_idx
#undef t_real
#undef check
