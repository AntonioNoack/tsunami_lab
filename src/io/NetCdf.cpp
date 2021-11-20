/**
 * @author Antonio Noack
 * 
 * @section DESCRIPTION
 * IO-routines for writing a snapshot as a NetCDF file.
 **/

#include "NetCdf.h"

#include <cmath> // isnan
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <netcdf.h>

#ifdef check
#error "already defined check()"
#endif

#define check(error) {\
  l_err = error;\
  if(l_err != NC_NOERR){\
    printf("NetCDF-Error occurred: %s (Code %d), line %d\n", nc_strerror(l_err), l_err, __LINE__);\
    return -1;\
  }\
}

// inspect a file: ncdump -h fileName, -h to only see the header
int tsunami_lab::io::NetCDF::appendTimeframe( t_real                       i_cellSizeMeters,
                                              t_idx                        i_nx,
                                              t_idx                        i_ny,
                                              t_idx                        i_step,
                                              t_idx                        i_stride,
                                              t_real               const * i_h,
                                              t_real               const * i_hu,
                                              t_real               const * i_hv,
                                              t_real               const * i_b,
                                              tsunami_lab::setups::Setup * i_setup,
                                              t_real                       i_time,
                                              t_idx                        i_frameIndex,
                                              std::string                  i_fileName ) {

  int l_err;
  
  bool l_isFirstFrame = i_frameIndex <= 0;
  
  int l_handle;// file handle
  
  if(l_isFirstFrame){
    check(nc_create(i_fileName.c_str(), NC_CLOBBER, &l_handle));
  } else {
    check(nc_open(i_fileName.c_str(), NC_WRITE, &l_handle));
  }
  
  // for output step (smaller output size)
  // in the future, we could average the value :)
  t_idx l_nx = i_nx / i_step, l_ny = i_ny / i_step;
  
  // dimension ids
  int l_xDimId, l_yDimId, l_tDimId;
  
  int l_heightId, l_momentumXId, l_momentumYId, l_bathymetryId, l_displacementId;
  int l_xVarId, l_yVarId; // values on x and y axes (e.g. 0 .. 50 x 0 .. 50)
  int l_tVarId;// time values in seconds
  
  std::vector<float> l_dataWithoutStride(l_nx * l_ny);
  
  // todo units, other axis descriptions
  // define dimensions
  if(l_isFirstFrame){
    check(nc_def_dim(l_handle, "x",    l_nx,         &l_xDimId));
    check(nc_def_dim(l_handle, "y",    l_ny,         &l_yDimId));
    check(nc_def_dim(l_handle, "time", NC_UNLIMITED, &l_tDimId));
    
    check(nc_def_var(l_handle, "x",    NC_FLOAT, 1, &l_xDimId, &l_xVarId));
    check(nc_def_var(l_handle, "y",    NC_FLOAT, 1, &l_yDimId, &l_yVarId));
    check(nc_def_var(l_handle, "time", NC_FLOAT, 1, &l_tDimId, &l_tVarId));
    
    check(nc_put_att_text(l_handle, l_xVarId, "units", 1, "m"));// strlen(DEGREES_NORTH), DEGREES_NORTH
    check(nc_put_att_text(l_handle, l_yVarId, "units", 1, "m"));
    check(nc_put_att_text(l_handle, l_tVarId, "units", 1, "s"));
    
    // define variables
    int dims3[3] = { l_tDimId, l_yDimId, l_xDimId };// fastest dimensions are last
    int dims2[2] = { l_yDimId, l_xDimId };// fastest dimensions are last
    
    if(i_h){
      check(nc_def_var(l_handle, "height",     NC_FLOAT, 3, dims3, &l_heightId));
      check(nc_put_att_text(l_handle, l_heightId, "units", 1, "m"));
    }
    if(i_hu){
      check(nc_def_var(l_handle, "momentumX",  NC_FLOAT, 3, dims3, &l_momentumXId));
      check(nc_put_att_text(l_handle, l_momentumXId, "units", 5, "m*m/s"));// height * velocity
    }
    if(i_hv){
      check(nc_def_var(l_handle, "momentumY",  NC_FLOAT, 3, dims3, &l_momentumYId));
      check(nc_put_att_text(l_handle, l_momentumYId, "units", 5, "m*m/s"));
    }
    if(i_b){
      check(nc_def_var(l_handle, "bathymetry", NC_FLOAT, 2, dims2, &l_bathymetryId));
      check(nc_put_att_text(l_handle, l_bathymetryId, "units", 1, "m"));
    }
    if(i_setup){
      check(nc_def_var(l_handle, "displacement", NC_FLOAT, 2, dims2, &l_displacementId));
      check(nc_put_att_text(l_handle, l_displacementId, "units", 1, "m"));
    }
    
    // end definition mode
    check(nc_enddef(l_handle));
    
    // nx*ny >= max(nx,ny) (because nx >= 1, ny >= 1), so we can reuse l_dataWithoutStride
    for(int i=0,j=0,l=std::max(l_nx,l_ny);j<l;i+=i_step,j++) l_dataWithoutStride[j] = (i+0.5) * i_cellSizeMeters;
    check(nc_put_var_float(l_handle, l_xVarId, l_dataWithoutStride.data()));// x and y axis share the same values, so we have to write them only once
    check(nc_put_var_float(l_handle, l_yVarId, l_dataWithoutStride.data()));
    
  } else {
    
    check(nc_inq_dimid(l_handle, "x",    &l_xDimId));
    check(nc_inq_dimid(l_handle, "y",    &l_yDimId));
    check(nc_inq_varid(l_handle, "time", &l_tVarId));
      
    if(i_h ) check(nc_inq_varid(l_handle, "height",     &l_heightId));
    if(i_hu) check(nc_inq_varid(l_handle, "momentumX",  &l_momentumXId));
    if(i_hv) check(nc_inq_varid(l_handle, "momentumY",  &l_momentumYId));
    
  }
  
  size_t l_zero = i_frameIndex;// index for that dimension (mmh)
  float l_floatTime = (float) i_time;// time may be a double
  check(nc_put_var1_float(l_handle, l_tVarId, &l_zero, &l_floatTime));
  
  size_t l_startVec[3] = { i_frameIndex, 0, 0 };
  size_t l_countVec[3] = { 1, l_ny, l_nx };// required when working with infinite dimensions
  
  if(i_h){
    for(t_idx y=0,o=0,my=i_ny-i_step+1;y<my;y+=i_step){
      for(t_idx x=0,mx=i_nx-i_step+1;x<mx;x+=i_step){
        t_idx l_idx = x + y * i_stride;
        l_dataWithoutStride[o++] = i_h[l_idx];
      }
    }
    check(nc_put_vara_float(l_handle, l_heightId, l_startVec, l_countVec, l_dataWithoutStride.data()));
  }
  
  if(i_hu){
    for(t_idx y=0,o=0,my=i_ny-i_step+1;y<my;y+=i_step){
      for(t_idx x=0,mx=i_nx-i_step+1;x<mx;x+=i_step){
        t_idx l_idx = x + y * i_stride;
        l_dataWithoutStride[o++] = i_hu[l_idx];
      }
    }
    check(nc_put_vara_float(l_handle, l_momentumXId, l_startVec, l_countVec, l_dataWithoutStride.data()));
  }
  
  if(i_hv){
    for(t_idx y=0,o=0,my=i_ny-i_step+1;y<my;y+=i_step){
      for(t_idx x=0,mx=i_nx-i_step+1;x<mx;x+=i_step){
        t_idx l_idx = x + y * i_stride;
        l_dataWithoutStride[o++] = i_hv[l_idx];
      }
    }
    check(nc_put_vara_float(l_handle, l_momentumYId, l_startVec, l_countVec, l_dataWithoutStride.data()));
  }
  
  if(l_isFirstFrame){
    if(i_b){
      for(t_idx y=0,o=0,my=i_ny-i_step+1;y<my;y+=i_step){
        for(t_idx x=0,mx=i_nx-i_step+1;x<mx;x+=i_step){
          t_idx l_idx = x + y * i_stride;
          l_dataWithoutStride[o++] = i_b[l_idx];
        }
      }
      // there is only a single frame of bathymetry, so we can use var instead of vara
      check(nc_put_var_float(l_handle, l_bathymetryId, l_dataWithoutStride.data()));
    }
    if(i_setup){
      t_real l_scaleX, l_scaleY;
      i_setup->getInitScale(l_scaleX, l_scaleY);
      for(t_idx y=0,o=0,my=i_ny-i_step+1;y<my;y+=i_step){
        t_real l_y = (y + (t_real) 0.5) * l_scaleY;
        for(t_idx x=0,mx=i_nx-i_step+1;x<mx;x+=i_step){
          t_real l_x = (x + (t_real) 0.5) * l_scaleX;
          l_dataWithoutStride[o++] = i_setup->getDisplacement(l_x, l_y);
        }
      }
      // there is only a single frame of displacement, so we can use var instead of vara
      check(nc_put_var_float(l_handle, l_displacementId, l_dataWithoutStride.data()));
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
                                          std::vector<t_real> & o_data) {
  
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
  
  // read in the main variable
  size_t l_startVec[2] = { 0, 0 };
  size_t l_countVec[2] = { o_sizeY, o_sizeX };// fastest dimensions are last
  if(sizeof(t_real) == 4){// this check may be doable in the pre-compiler somehow as well
    check(nc_get_vara_float(l_handle, l_zVarId, l_startVec, l_countVec, (float*) o_data.data()));
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

#undef check