/**
 * @author Antonio Noack
 * 
 * @section DESCRIPTION
 * IO-routines for writing a snapshot as a NetCDF file.
 **/

#include <cmath> // isnan
#include <iostream> // std::cerr
#include <fstream>
#include <cstdio> // std::rename
#include <stdexcept>
#include <netcdf.h>
#include <sys/stat.h> // check whether a file exists

#include "NetCdf.h"
#include "../setups/CheckPoint.h"

#ifdef check
#error "already defined check()"
#endif

#define check(error) {\
  l_err = error;\
  if(l_err != NC_NOERR){\
    std::cerr << "NetCDF-Error occurred: " << nc_strerror(l_err) << " (Code " << l_err << "), line " << __LINE__ << std::endl;\
    nc_close(l_handle);\
    return -1;\
  }\
}

#define check2(error) {\
  l_err = error;\
  if(l_err != NC_NOERR){\
    std::cerr << "NetCDF-Error occurred: " << nc_strerror(l_err) << " (Code " << l_err << "), line " << __LINE__ << std::endl;\
    nc_close(l_handle);\
    return nullptr;\
  }\
}

#define check3(error) {\
  l_err = error;\
  if(l_err != NC_NOERR){\
    std::cerr << "NetCDF-Error occurred: " << nc_strerror(l_err) << " (Code " << l_err << "), line " << __LINE__ << std::endl;\
    nc_close(l_handle);\
    return;\
  }\
}

// functions for t_real
#define put_vara(handle, varId, start, count, data)\
  sizeof(t_real) == 4 ? nc_put_vara_float(handle, varId, start, count, (float*) data) : nc_put_vara_double(handle, varId, start, count, (double*) data)

#define get_vara(handle, varId, start, count, data)\
  sizeof(t_real) == 4 ? nc_get_vara_float(handle, varId, start, count, (float*) data) : nc_get_vara_double(handle, varId, start, count, (double*) data)

#define put_var1(handle, varId, start, data)\
  sizeof(t_real) == 4 ? nc_put_var1_float(handle, varId, start, (float*) data) : nc_put_var1_double(handle, varId, start, (double*) data)

#define get_var1(handle, varId, start, data)\
  sizeof(t_real) == 4 ? nc_get_var1_float(handle, varId, start, (float*) data) : nc_get_var1_double(handle, varId, start, (double*) data)

#ifndef CEIL_DIV
#define CEIL_DIV(a,div) (a+div-1)/(div)
#endif

int tsunami_lab::io::NetCDF::storeRow( int l_handle,
                                       int i_varId,
                                       int i_timeIndex,
                                       t_idx i_yOut,
                                       t_idx i_sizeXOut,
                                       t_real const* i_dataOut ){
  if(i_timeIndex < 0){// no time axis present
    size_t l_start[2] = { i_yOut, 0 };
    size_t l_count[2] = { 1, i_sizeXOut };
    return put_vara(l_handle, i_varId, l_start, l_count, i_dataOut);
  } else {
    size_t l_start[3] = { (size_t) i_timeIndex, i_yOut, 0 };
    size_t l_count[3] = { 1, 1, i_sizeXOut };
    return put_vara(l_handle, i_varId, l_start, l_count, i_dataOut);
  }
}

int tsunami_lab::io::NetCDF::downsample( int l_handle,
                                         int i_varId,
                                         int i_timeIndex,
                                         t_idx         i_sizeXIn,
                                         t_idx         i_sizeYIn,
                                         t_idx         i_strideIn,
                                         t_real const* i_dataIn,
                                         t_idx         i_sizeXOut,
                                         t_idx         i_sizeYOut,
                                         t_real*       i_dataOut,
                                         t_idx         i_step ){
  int l_err = 0;
  if(i_step > 1){
    std::cout << "    writing " << i_sizeXIn << " x " << i_sizeYIn << " -> " << i_sizeXOut << " x " << i_sizeYOut << " with interpolation" << std::endl;
    for(t_idx l_yOut=0;l_yOut<i_sizeYOut;l_yOut++){
      const t_idx l_yIn0 = l_yOut * i_step;
      const t_idx l_yIn1 = std::min(l_yIn0 + i_step, i_sizeYIn);
      // the outer loop cannot be parallelized that easily, because NetCDF may not be thread safe
      // this inner loop, as small as it may seem, might be pretty large (e.g. 54k elements)
      #pragma omp parallel for
      for(t_idx l_xOut=0;l_xOut<i_sizeXOut;l_xOut++){
        const t_idx l_xIn0 = l_xOut * i_step;
        const t_idx l_xIn1 = std::min(l_xIn0 + i_step, i_sizeXIn);
        t_real l_sum = 0;
        for(t_idx l_yIn=l_yIn0;l_yIn<l_yIn1;l_yIn++){
          t_idx l_indexIn = l_xIn0 + l_yIn * i_strideIn;
          for(t_idx l_xIn=l_xIn0;l_xIn<l_xIn1;l_xIn++){
            l_sum += i_dataIn[l_indexIn++];
          }
        }
        i_dataOut[l_xOut] = l_sum / (t_real)((l_xIn1-l_xIn0)*(l_yIn1-l_yIn0));
      }
      l_err = storeRow(l_handle, i_varId, i_timeIndex, l_yOut, i_sizeXOut, i_dataOut);
      if(l_err) return l_err;
    }
  } else {// just copy the stripes
    if(i_strideIn == i_sizeXOut){// just a pure copy is required
      std::cout << "    writing " << i_sizeXOut << " x " << i_sizeYOut << " in one block" << std::endl;
      if(i_timeIndex < 0){// no time axis present
        size_t l_start[2] = { 0, 0 };
        size_t l_count[2] = { i_sizeYOut, i_sizeXOut };
        l_err = put_vara(l_handle, i_varId, l_start, l_count, i_dataIn);
      } else {
        size_t l_start[3] = { (size_t) i_timeIndex, 0, 0 };
        size_t l_count[3] = { 1, i_sizeYOut, i_sizeXOut };
        l_err = put_vara(l_handle, i_varId, l_start, l_count, i_dataIn);
      }
    } else {
      std::cout << "    writing " << i_sizeXOut << " x " << i_sizeYOut << " in stripes" << std::endl;
      if(i_timeIndex < 0){
        size_t l_start[2] = { 0, 0 };
        size_t l_count[2] = { i_sizeYOut, i_sizeXOut };
        ptrdiff_t l_stride[2] = { 1, 1 };
        ptrdiff_t l_map[2] = { (ptrdiff_t) i_strideIn, 1 };
        l_err = nc_put_varm_float(l_handle, i_varId, l_start, l_count, l_stride, l_map, i_dataIn);
      } else {
        size_t l_start[3] = { (size_t) i_timeIndex, 0, 0 };
        size_t l_count[3] = { 1, i_sizeYOut, i_sizeXOut };
        ptrdiff_t l_stride[3] = { 1, 1, 1 };
        ptrdiff_t l_map[3] = { 0, (ptrdiff_t) i_strideIn, 1 };
        l_err = nc_put_varm_float(l_handle, i_varId, l_start, l_count, l_stride, l_map, i_dataIn);
      }
    }
    // old code, that may be used for performance comparisons
    /*for(t_idx l_yOut=0;l_yOut<i_sizeYOut;l_yOut++){
      // std::cout << "writing " << l_yOut << "/" << i_sizeYOut << ", " << i_sizeXOut << " values" << std::endl;
      t_idx l_indexIn  = l_yOut * i_strideIn;
      l_err = storeRow(l_handle, i_varId, i_timeIndex, l_yOut, i_sizeXOut, i_dataIn + l_indexIn);
      if(l_err) return l_err;
    }*/
  }
  std::cout << "    done writing chunk of memory" << std::endl;
  return EXIT_SUCCESS;
}

tsunami_lab::setups::Setup* tsunami_lab::io::NetCDF::loadCheckpoint( std::string i_fileName, t_idx &o_nx, t_idx &o_ny, t_real &o_cellSizeMeters, t_real &o_cflFactor, double &o_simulationTime, t_idx &o_timeStepIndex, std::vector<tsunami_lab::io::Station> &o_stations ){
  
  bool f = sizeof(t_real) == 4;

  int l_err, l_handle;
  check2(nc_open(i_fileName.c_str(), NC_NOWRITE, &l_handle));
  
  // query the dimension ids
  int l_xDimId, l_yDimId, l_iDimId;
  check2(nc_inq_dimid(l_handle, "x", &l_xDimId));
  check2(nc_inq_dimid(l_handle, "y", &l_yDimId));
  check2(nc_inq_dimid(l_handle, "i", &l_iDimId));
  
  // query dimension sizes
  char l_tmpDimName[NC_MAX_NAME+1];// name, ignored
  check2(nc_inq_dim(l_handle, l_xDimId, l_tmpDimName, &o_nx));
  check2(nc_inq_dim(l_handle, l_yDimId, l_tmpDimName, &o_ny));
  
  // query the variable ids, for which the values are vectors
  int l_hVarId, l_bVarId, l_huVarId, l_hvVarId;
  check2(nc_inq_varid(l_handle, "height",     &l_hVarId));
  check2(nc_inq_varid(l_handle, "bathymetry", &l_bVarId));
  check2(nc_inq_varid(l_handle, "momentumX",  &l_huVarId));
  if(o_ny > 1) check2(nc_inq_varid(l_handle, "momentumY",  &l_hvVarId));

  // query the variable ids, for which the values are scalars
  int l_simTimeVarId, l_simIdxVarId, l_cellSizeVarId, l_cflFactorVarId;
  check2(nc_inq_varid(l_handle, "simulationTime",  &l_simTimeVarId));
  check2(nc_inq_varid(l_handle, "timeStepIndex",   &l_simIdxVarId));
  check2(nc_inq_varid(l_handle, "cflFactor",       &l_cflFactorVarId));
  check2(nc_inq_varid(l_handle, "cellSizeMeters",  &l_cellSizeVarId));
  
  long long int o_timeStepIndex2;
  check2(nc_get_var1_double(  l_handle, l_simTimeVarId, nullptr, (double*)        &o_simulationTime));
  check2(nc_get_var1_longlong(l_handle, l_simIdxVarId,  nullptr, (long long int*) &o_timeStepIndex2));
  o_timeStepIndex = o_timeStepIndex2;// just in case the size is different
  
  // read all scalar properties
  if(f){
    check2(nc_get_var1_float( l_handle, l_cflFactorVarId, nullptr, (float*)  &o_cflFactor));
    check2(nc_get_var1_float( l_handle, l_cellSizeVarId,  nullptr, (float*)  &o_cellSizeMeters));
  } else {
    check2(nc_get_var1_double(l_handle, l_cflFactorVarId, nullptr, (double*) &o_cflFactor));
    check2(nc_get_var1_double(l_handle, l_cellSizeVarId,  nullptr, (double*) &o_cellSizeMeters));
  }
  
  size_t l_size = sizeof(t_real) * o_nx * o_ny;
  t_real* l_h  = new t_real[l_size];
  t_real* l_b  = new t_real[l_size];
  t_real* l_hu = new t_real[l_size];
  t_real* l_hv = new t_real[l_size];
  
  // we cannot use check() in this section, because we would leak memory; therefore just collect the errors,
  // and hope that it's the same error, and then print it
  
  // read all vector values
  size_t l_startVec[2] = { 0, 0 };
  size_t l_countVec[2] = { o_ny, o_nx };// fastest dimensions are last
  check2(get_vara(l_handle, l_hVarId,  l_startVec, l_countVec, l_h));
  check2(get_vara(l_handle, l_bVarId,  l_startVec, l_countVec, l_b));
  check2(get_vara(l_handle, l_huVarId, l_startVec, l_countVec, l_hu));
  if(o_ny > 1)
    check2(get_vara( l_handle, l_hvVarId, l_startVec, l_countVec, l_hv));
  
  // todo query stations:
  // todo query all variables starting with "station-", then read in their values
  // for better serialization, we could use netcdf compounds in the future
  (void) o_stations;
  int l_numVariables;
  check2(nc_inq(l_handle, nullptr, &l_numVariables, nullptr, nullptr));
  size_t l_dimILength;
  check2(nc_inq_dimlen(l_handle, l_iDimId, &l_dimILength));
  std::vector<t_real> l_stationData(l_dimILength);
  for(int l_var=0;l_var<l_numVariables;l_var++){
    char l_name[NC_MAX_NAME+1];
    check2(nc_inq_var(l_handle, l_var, l_name, nullptr, nullptr, nullptr, nullptr));
    if(strncmp(l_name, "station-", strlen("station-")) == 0){// starts with "station-", so it's a station
      
      std::string l_stationName(l_name + strlen("station-"));// variable name without prefix
      
      // read in all station data
      size_t l_startVec0[1] = { 0 };
      size_t l_countVec0[1] = { l_dimILength };
      check2(get_vara(l_handle, l_var, l_startVec0, l_countVec0, l_stationData.data()));
      
      t_real l_delayBetweenRecords = l_stationData[0];
      t_idx  l_x = (t_idx) l_stationData[1];
      t_idx  l_y = (t_idx) l_stationData[2];
      tsunami_lab::io::Station l_station(l_x, l_y, l_stationName, l_delayBetweenRecords);
      for(size_t i=3,l=l_dimILength-3;i<l;){
        t_real l_t  = l_stationData[i++];
        t_real l_h  = l_stationData[i++];
        t_real l_hu = l_stationData[i++];
        t_real l_hv = l_stationData[i++];
        l_station.recordState(l_t, l_h, l_hu, l_hv);
      }
      o_stations.push_back(l_station);
      
    }
  }
  
  if(l_err){
    delete[] l_h;
    delete[] l_b;
    delete[] l_hu;
    delete[] l_hv;
    check2(l_err);
  }
  
  check2(nc_close(l_handle));
  
  return new tsunami_lab::setups::CheckPoint(l_h, l_b, l_hu, l_hv, o_nx, o_ny, o_nx);
  
}

int tsunami_lab::io::NetCDF::storeCheckpoint( std::string i_fileName, t_idx i_nx, t_idx i_ny, t_real i_cellSizeMeters, t_real i_cflFactor, double i_simulationTime, t_idx i_timeStepIndex, std::vector<tsunami_lab::io::Station> &i_stations, tsunami_lab::patches::WavePropagation* i_waveProp ){
  
  #ifdef MEMORY_IS_SCARCE
  int l_deflateLevel = 0;
  #else
  int l_deflateLevel = 2;
  #endif
  bool l_deflate = l_deflateLevel > 0;
  bool l_shuffle = l_deflate;
  
  // for the ghost cells
  size_t offset;
  
  if(i_ny > 1){
    // 2d
    i_nx += 2;
    i_ny += 2;
    offset = i_nx + 1;
  } else {
    // 1d
    i_nx += 2;
    i_ny = 1;
    offset = 1;
  }
  
  int l_err, l_handle;
  
  // if file exists, create tmp file
  struct stat buffer;
  bool l_fileExisted = stat(i_fileName.c_str(), &buffer) == 0;
  std::string l_tmpFileName = i_fileName + ".tmp";
  
  if(l_fileExisted) std::cout << "  writing to " << l_tmpFileName << std::endl;
  check(nc_create((l_fileExisted ? l_tmpFileName : i_fileName).c_str(), NC_CLOBBER | NC_NETCDF4, &l_handle));
  
  int l_xDimId, l_yDimId, l_iDimId;
  int l_stationDimSize = i_stations.size() > 0 ? i_stations[0].getRecords().size() * 4 + 3 : 1;
  check(nc_def_dim(l_handle, "x", i_nx, &l_xDimId));
  check(nc_def_dim(l_handle, "y", i_ny, &l_yDimId));
  check(nc_def_dim(l_handle, "i", l_stationDimSize, &l_iDimId));// meaningless dimension; used for easy serialization of stations
  // die Koordinatenachsen brauchen keine Werte (Variablen)
  
  auto l_type = sizeof(t_real) == 4 ? NC_FLOAT : NC_DOUBLE;
  
  int dims2[2] = { l_yDimId, l_xDimId };
  int l_hVarId, l_bVarId, l_huVarId, l_hvVarId;
  
  check(nc_def_var(l_handle, "height",     l_type, 2, dims2, &l_hVarId));
  check(nc_def_var_deflate(l_handle, l_hVarId,  l_shuffle, l_deflate, l_deflateLevel));
  
  check(nc_def_var(l_handle, "bathymetry", l_type, 2, dims2, &l_bVarId));
  check(nc_def_var_deflate(l_handle, l_bVarId,  l_shuffle, l_deflate, l_deflateLevel));
  
  check(nc_def_var(l_handle, "momentumX",  l_type, 2, dims2, &l_huVarId));
  check(nc_def_var_deflate(l_handle, l_huVarId, l_shuffle, l_deflate, l_deflateLevel));
  
  if(i_waveProp->getMomentumY()){
    check(nc_def_var(l_handle, "momentumY",  l_type, 2, dims2, &l_hvVarId));
    check(nc_def_var_deflate(l_handle, l_hvVarId, l_shuffle, l_deflate, l_deflateLevel));
  }
  
  // für die 0d-Variablen hingegen bringt es nichts
  int l_simTimeVarId, l_simIdxVarId, l_cellSizeVarId, l_cflFactorVarId;
  check(nc_def_var(l_handle, "simulationTime", NC_DOUBLE, 0, nullptr, &l_simTimeVarId));
  check(nc_def_var(l_handle, "timeStepIndex",  NC_INT64,  0, nullptr, &l_simIdxVarId));
  check(nc_def_var(l_handle, "cellSizeMeters", l_type,    0, nullptr, &l_cellSizeVarId));
  check(nc_def_var(l_handle, "cflFactor",      l_type,    0, nullptr, &l_cflFactorVarId));
  
  // station variables; one per station (not the best, but maybe the easiest)
  std::vector<int> l_stationVarIds(i_stations.size());
  int l_i = 0;
  for(auto &l_station : i_stations){
    int l_stationVarId;
    std::string l_joinedName = "station-" + l_station.getName();// station prefix, in case there is a station called like a variable
    check(nc_def_var(l_handle, l_joinedName.c_str(), l_type, 1, &l_iDimId, &l_stationVarId));
    l_stationVarIds[l_i++] = l_stationVarId;
  }
  
  check(nc_enddef(l_handle));
  
  t_idx l_stride = i_waveProp->getStride();
  
  // 2d-Variablen
  check(downsample(l_handle,   l_hVarId,  -1, i_nx, i_ny, l_stride, i_waveProp->getHeight()     - offset, i_nx, i_ny, nullptr, 1));
  check(downsample(l_handle,   l_bVarId,  -1, i_nx, i_ny, l_stride, i_waveProp->getBathymetry() - offset, i_nx, i_ny, nullptr, 1));
  check(downsample(l_handle,   l_huVarId, -1, i_nx, i_ny, l_stride, i_waveProp->getMomentumX()  - offset, i_nx, i_ny, nullptr, 1));
  if(i_waveProp->getMomentumY()){
    check(downsample(l_handle, l_hvVarId, -1, i_nx, i_ny, l_stride, i_waveProp->getMomentumY()  - offset, i_nx, i_ny, nullptr, 1));
  }
  
  // 0d-Variablen
  long long int i_timeStepIndex2 = i_timeStepIndex;// in case sizeof(long long int) != 8
  check(nc_put_var1_double(  l_handle, l_simTimeVarId, nullptr, (double*)        &i_simulationTime));
  check(nc_put_var1_longlong(l_handle, l_simIdxVarId,  nullptr, (long long int*) &i_timeStepIndex2));
  
  check(put_var1(l_handle, l_cellSizeVarId,  nullptr, &i_cellSizeMeters));
  check(put_var1(l_handle, l_cflFactorVarId, nullptr, &i_cflFactor));
  
  // write stations as heterogenous arrays
  // (could be done with multiple variables per station, but I'm lazy for this one)
  l_i = 0;
  for(auto &l_station : i_stations){
    int l_var = l_stationVarIds[l_i++];
    size_t l_index = 0;
    size_t l_size = 1;// one value at a time
    t_real l_value = 0;
    
    #define putValue(value)\
      l_value = (t_real) value;\
      check(put_vara(l_handle, l_var, &l_index, &l_size, &l_value));\
      l_index++;\
    
    putValue(l_station.getDelayBetweenRecords());
    
    t_idx l_posX, l_posY;
    l_station.getPosition(l_posX, l_posY);
    putValue(l_posX);
    putValue(l_posY);
    
    auto l_records = l_station.getRecords();
    for(auto &l_record : l_records){
      putValue(l_record.time);
      putValue(l_record.height);
      putValue(l_record.momentumX);
      putValue(l_record.momentumY);
    }
    
  }
  
  check(nc_close(l_handle));
  
  // if file existed, rm original, mv file
  if(l_fileExisted){
    l_err = std::remove(i_fileName.c_str());
    if(l_err) std::cerr << "Remove of " << i_fileName << " failed!" << std::endl;
    l_err = std::rename(l_tmpFileName.c_str(), i_fileName.c_str());
    if(l_err) std::cerr << "Rename from " << l_tmpFileName << " to " << i_fileName << " failed!" << std::endl;
  }
  
  return 0;
}

// inspect a file: ncdump -h fileName, -h to only see the header, -c to see header + axis data
int tsunami_lab::io::NetCDF::appendTimeframe( t_real                       i_cellSizeMeters,
                                              t_idx                        i_nx,
                                              t_idx                        i_ny,
                                              t_real                       i_gridOffsetX,
                                              t_real                       i_gridOffsetY,
                                              t_idx                        i_step,
                                              t_idx                        i_stride,
                                              t_real               const * i_h,
                                              t_real               const * i_hu,
                                              t_real               const * i_hv,
                                              t_real               const * i_b,
                                              tsunami_lab::setups::Setup * i_setup,
                                              t_real                       i_time,
                                              int                          i_deflateLevel,
                                              std::string                  i_fileName ) {

  int l_err;
  
  bool l_deflate = i_deflateLevel > 0;
  
  // reorders high and low bytes such that first all high bytes are written, then the low bytes.
  // is said to be useless for deflate level = 9. In my test, testing config/tsunami2d-tohoku.yaml, this was incorrect
  bool l_shuffle = l_deflate;
  
  bool l_isFirstFrame = i_time <= 0;
  
  int l_handle;// file handle
  
  // coarse output field size x/y if i_step > 1
  t_idx l_nx = CEIL_DIV(i_nx, i_step), l_ny = CEIL_DIV(i_ny, i_step);
  
  // dimension ids
  int l_xDimId, l_yDimId, l_tDimId;
  int l_heightId, l_momentumXId, l_momentumYId, l_bathymetryId, l_displacementId;
  int l_xVarId, l_yVarId; // values on x and y axes (e.g. 0 .. 50 x 0 .. 50)
  int l_tVarId;// time values in seconds
  
  size_t i_timeStepIndex = 0;
  
  if(l_isFirstFrame){
    check(nc_create(i_fileName.c_str(), NC_CLOBBER | NC_NETCDF4, &l_handle));
  } else {
    check(nc_open(i_fileName.c_str(), NC_WRITE, &l_handle));
  }
  
  std::vector<float> l_dataWithoutStride(std::max(l_nx, l_ny));
  
  // todo units, other axis descriptions
  // define dimensions
  if(l_isFirstFrame){
    
    check(nc_put_att_text(l_handle, NC_GLOBAL, "Conventions", 6, "COARDS"));
    // when we have meaningful x/y coordinates, we can define a long_name for them (https://ferret.pmel.noaa.gov/Ferret/documentation/coards-netcdf-conventions)
    
    check(nc_def_dim(l_handle, "x",    l_nx,         &l_xDimId));
    check(nc_def_dim(l_handle, "y",    l_ny,         &l_yDimId));
    check(nc_def_dim(l_handle, "time", NC_UNLIMITED, &l_tDimId));
    
    check(nc_def_var(l_handle, "x",    NC_FLOAT, 1, &l_xDimId, &l_xVarId));
    check(nc_def_var(l_handle, "y",    NC_FLOAT, 1, &l_yDimId, &l_yVarId));
    check(nc_def_var(l_handle, "time", NC_FLOAT, 1, &l_tDimId, &l_tVarId));
    
    check(nc_def_var_deflate(l_handle, l_xVarId, l_shuffle, l_deflate, i_deflateLevel));
    check(nc_def_var_deflate(l_handle, l_yVarId, l_shuffle, l_deflate, i_deflateLevel));
    check(nc_def_var_deflate(l_handle, l_tVarId, l_shuffle, l_deflate, i_deflateLevel));
    
    check(nc_put_att_text(l_handle, l_xVarId, "units", 1, "m"));// strlen(DEGREES_NORTH), DEGREES_NORTH
    check(nc_put_att_text(l_handle, l_yVarId, "units", 1, "m"));
    check(nc_put_att_text(l_handle, l_tVarId, "units", 1, "s"));
    
    // define variables
    int dims3[3] = { l_tDimId, l_yDimId, l_xDimId };// fastest dimensions are last
    int dims2[2] = { l_yDimId, l_xDimId };// fastest dimensions are last
    
    if(i_h){
      check(nc_def_var(l_handle, "height",     NC_FLOAT, 3, dims3, &l_heightId));
      check(nc_put_att_text(l_handle, l_heightId, "units", 1, "m"));
      check(nc_def_var_deflate(l_handle, l_heightId, l_shuffle, l_deflate, i_deflateLevel));
    }
    if(i_hu){
      check(nc_def_var(l_handle, "momentumX",  NC_FLOAT, 3, dims3, &l_momentumXId));
      check(nc_put_att_text(l_handle, l_momentumXId, "units", 5, "m*m/s"));// height * velocity
      check(nc_def_var_deflate(l_handle, l_momentumXId, l_shuffle, l_deflate, i_deflateLevel));
    }
    if(i_hv){
      check(nc_def_var(l_handle, "momentumY",  NC_FLOAT, 3, dims3, &l_momentumYId));
      check(nc_put_att_text(l_handle, l_momentumYId, "units", 5, "m*m/s"));
      check(nc_def_var_deflate(l_handle, l_momentumYId, l_shuffle, l_deflate, i_deflateLevel));
    }
    if(i_b){
      check(nc_def_var(l_handle, "bathymetry", NC_FLOAT, 2, dims2, &l_bathymetryId));
      check(nc_put_att_text(l_handle, l_bathymetryId, "units", 1, "m"));
      check(nc_def_var_deflate(l_handle, l_bathymetryId, l_shuffle, l_deflate, i_deflateLevel));
    }
    if(i_setup){
      check(nc_def_var(l_handle, "displacement", NC_FLOAT, 2, dims2, &l_displacementId));
      check(nc_put_att_text(l_handle, l_displacementId, "units", 1, "m"));
      check(nc_def_var_deflate(l_handle, l_displacementId, l_shuffle, l_deflate, i_deflateLevel));
    }
    
    // end definition mode
    check(nc_enddef(l_handle));
    
    // nx*ny >= max(nx,ny) (because nx >= 1, ny >= 1), so we can reuse l_dataWithoutStride
    // once + 0.5 * cellSizeMeters, but since gridOffset now will contain the values from the loaded bathymetry for tsunamis,
    // we will just return the same coordinates as the input file
    for(int i=0,j=0,l=std::max(l_nx,l_ny);j<l;i+=i_step,j++) l_dataWithoutStride[j] = i * i_cellSizeMeters + i_gridOffsetX;
    check(nc_put_var_float(l_handle, l_xVarId, l_dataWithoutStride.data()));
    
    for(int i=0,j=0,l=std::max(l_nx,l_ny);j<l;i+=i_step,j++) l_dataWithoutStride[j] = i * i_cellSizeMeters + i_gridOffsetY;
    check(nc_put_var_float(l_handle, l_yVarId, l_dataWithoutStride.data()));
    
  } else {
    
    check(nc_inq_dimid(l_handle, "x",    &l_xDimId));
    check(nc_inq_dimid(l_handle, "y",    &l_yDimId));
    check(nc_inq_dimid(l_handle, "time", &l_tDimId));
    check(nc_inq_varid(l_handle, "time", &l_tVarId));
    
    // query maximum time, and just return, if the maximum value is already bigger
    // also query the time index of this current frame
    check(nc_inq_dimlen(l_handle, l_tDimId, &i_timeStepIndex));
    
    float l_lastFloatTime;
    size_t l_index[1] = { i_timeStepIndex-1 };
    check(nc_get_var1_float(l_handle, l_tVarId, l_index, &l_lastFloatTime));
    
    if((float) i_time <= l_lastFloatTime){
      std::cout << "time was already found in file: " << i_time << " <= " << l_lastFloatTime << "; skipping writing frame" << std::endl;
      check(nc_close(l_handle));
      return EXIT_SUCCESS;
    }
    
    if(i_h ) check(nc_inq_varid(l_handle, "height",     &l_heightId));
    if(i_hu) check(nc_inq_varid(l_handle, "momentumX",  &l_momentumXId));
    if(i_hv) check(nc_inq_varid(l_handle, "momentumY",  &l_momentumYId));
    
  }
  
  size_t l_timeStepIndex = i_timeStepIndex;// index for that dimension
  float l_floatTime = (float) i_time;// time may be a double
  check(nc_put_var1_float(l_handle, l_tVarId, &l_timeStepIndex, &l_floatTime));
  
  if(i_h){  check(downsample(l_handle, l_heightId,    i_timeStepIndex, i_nx, i_ny, i_stride, i_h,  l_nx, l_ny, l_dataWithoutStride.data(), i_step)); }
  if(i_hu){ check(downsample(l_handle, l_momentumXId, i_timeStepIndex, i_nx, i_ny, i_stride, i_hu, l_nx, l_ny, l_dataWithoutStride.data(), i_step)); }
  if(i_hv){ check(downsample(l_handle, l_momentumYId, i_timeStepIndex, i_nx, i_ny, i_stride, i_hv, l_nx, l_ny, l_dataWithoutStride.data(), i_step)); }
  
  if(l_isFirstFrame){
    if(i_b){ check(downsample(l_handle, l_bathymetryId, -1, i_nx, i_ny, i_stride, i_b, l_nx, l_ny, l_dataWithoutStride.data(), i_step)); }
    if(i_setup){
      // this is sub-ideal, and would ideally require downsampling
      // however, for that we'd need to allocate a new, large buffer
      // or create a second downsample function
      t_real l_scaleX, l_scaleY;
      i_setup->getInitScale(l_scaleX, l_scaleY);
      for(t_idx y=0;y<l_ny;y++){
        t_real l_y = ((y * i_step) + (t_real) 0.5) * l_scaleY;
        for(t_idx x=0;x<l_nx;x++){
          t_real l_x = ((x * i_step) + (t_real) 0.5) * l_scaleX;
          l_dataWithoutStride[x] = i_setup->getDisplacement(l_x, l_y);
        }
        check(storeRow(l_handle, l_displacementId, -1, y, l_nx, l_dataWithoutStride.data()));
      }
    }
  }
  
  check(nc_close(l_handle));
  
  return EXIT_SUCCESS;
  
}

int tsunami_lab::io::NetCDF::load2dArray( std::string           i_fileName,
                                          std::string           i_variableName,
                                          t_idx               & o_sizeX,
                                          t_idx               & o_sizeY,
                                          t_real              & o_cellSizeMeters,
                                          t_real              & o_gridOffsetX,
                                          t_real              & o_gridOffsetY,
                                          std::vector<t_real> & o_data ) {
  
  int l_err;
  
  int l_handle; // file handle
  
  check(nc_open(i_fileName.c_str(), NC_NOWRITE, &l_handle));
  
  // query the dimension ids
  int l_xDimId, l_yDimId;
  check(nc_inq_dimid(l_handle, "x", &l_xDimId));
  check(nc_inq_dimid(l_handle, "y", &l_yDimId));
  
  int l_xVarId, l_yVarId, l_zVarId;
  check(nc_inq_varid(l_handle, i_variableName.c_str(), &l_zVarId)); 
  
  check(nc_inq_varid(l_handle, "x", &l_xVarId));
  check(nc_inq_varid(l_handle, "y", &l_yVarId));
  
  // query dimension sizes
  char l_tmpDimName[NC_MAX_NAME+1];// name, ignored
  check(nc_inq_dim(l_handle, l_xDimId, l_tmpDimName, &o_sizeX));
  check(nc_inq_dim(l_handle, l_yDimId, l_tmpDimName, &o_sizeY));
  
  // allocate space for data
  o_data.resize(o_sizeX * o_sizeY);
  
  // compute the cell size by using the values in the x variable
  // non-square cells are not supported currently
  size_t l_startVec0[1] = { 0 };
  size_t l_countVec0[1] = { o_sizeX };
  check(nc_get_vara_float(l_handle, l_xVarId, l_startVec0, l_countVec0, (float*) o_data.data()));
  
  o_cellSizeMeters = (o_data[o_sizeX-1] - o_data[0]) / (o_sizeX - 1);
  o_gridOffsetX = o_data[1];// the first one is used as ghost zone
  
  l_countVec0[0] = 2;
  check(nc_get_vara_float(l_handle, l_yVarId, l_startVec0, l_countVec0, (float*) o_data.data()));
  o_gridOffsetY = o_data[1];// the first one is used as ghost zone
  
  // read in the main variable
  size_t l_startVec[2] = { 0, 0 };
  size_t l_countVec[2] = { o_sizeY, o_sizeX };// fastest dimensions are last
  if(sizeof(t_real) == 4){// this check may be doable in the precompiler somehow as well
    check(nc_get_vara_float( l_handle, l_zVarId, l_startVec, l_countVec, (float*)  o_data.data()));
  } else {
    check(nc_get_vara_double(l_handle, l_zVarId, l_startVec, l_countVec, (double*) o_data.data()));
  }
  
  // check for NaNs, just so we know what's up
  for(t_idx i=0,l=o_sizeX*o_sizeY;i<l;i++){
    if(std::isnan(o_data[i])){
      t_idx j=i+1;
      for(;j<l && std::isnan(o_data[j]);j++){}
      std::cerr << "data contains NaN from " << i << " to " << j << " from 0 .. " << l << " in field of size " << o_sizeX << " x " << o_sizeY << std::endl;
      i = j-1;
    }
  }
  
  check(nc_close(l_handle));
  
  return EXIT_SUCCESS;
  
}
