/**
 * @author Antonio Noack
 * @section DESCRIPTION
 * Tests the checkpoint setup.
 **/
#include <catch2/catch.hpp>

#include "CheckPoint.h"
#include "ArtificialTsunami2d.h"
#include "../constants.h"
#include "../patches/WavePropagation2d.h"
#include "../io/NetCdf.h"

#define t_real tsunami_lab::t_real
#define t_idx  tsunami_lab::t_idx

/*TEST_CASE( "Tests for 2d tsunami event.", "[CheckPoint]" ) {

  t_idx l_nx = 1000, l_ny = 1000;
  
  tsunami_lab::setups::ArtificialTsunami2d l_baselineSetup;
  
  t_idx l_nx2 = 0, l_ny2 = 0, l_nx3 = 0, l_ny3 = 0;
  t_real l_cellSizeMeters2 = 0, l_cellSizeMeters3 = 0, l_tmp;
  
  std::vector<t_real> l_bath;
  std::vector<t_real> l_disp;
  
  tsunami_lab::io::NetCDF::load2dArray("data/artificialtsunami_bathymetry_1000.nc", "z", l_nx2, l_ny2, l_cellSizeMeters2, l_tmp, l_tmp, l_bath);// 100 x 100 grid
  tsunami_lab::io::NetCDF::load2dArray("data/artificialtsunami_displ_1000.nc", "z", l_nx3, l_ny3, l_cellSizeMeters3, l_tmp, l_tmp, l_disp);// 1000 x 1000 grid
  
  tsunami_lab::setups::CheckPoint l_testedSetup(l_bath.data(), l_nx2, l_ny2, l_nx2, 10, l_disp.data(), l_nx3, l_ny3, l_nx3, 1);
  
  tsunami_lab::patches::WavePropagation2d l_baselineProp(l_nx, l_ny, &l_baselineSetup, 1, 1);
  tsunami_lab::patches::WavePropagation2d l_testProp(l_nx, l_ny, &l_testedSetup, 1, 1);
  
  const t_real* l_baselineH = l_baselineProp.getHeight();
  const t_real* l_testH = l_testProp.getHeight();
  
  const t_real* l_baselineB = l_baselineProp.getBathymetry();
  const t_real* l_testB = l_testProp.getBathymetry();
  
  t_idx l_stride = l_testProp.getStride();
  
  for(t_idx y=0;y<l_ny-2;y++){// -2 for the ghost zone
    for(t_idx x=0;x<l_nx-2;x++){
      t_idx l_index = l_nx + l_ny * l_stride;
	  // there is a small difference, but probably only because the grid is defined slightly differently (e.g. using corners instead of cell centers, or filling the ghost zones with distinct data/just copying it)
      REQUIRE(l_baselineH[l_index] == Approx(l_testH[l_index]).margin(0.01));
      REQUIRE(l_baselineB[l_index] == Approx(l_testB[l_index]).margin(0.01));
    }
  }
  
}*/
