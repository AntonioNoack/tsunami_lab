/**
 * @author Antonio Noack
 * @section DESCRIPTION
 * F-Wave solver for the shallow water equations.
 **/
#include "FWave.h"
#include <cmath>
#include <stdio.h>
#include <iostream>
#include <stdexcept>

void tsunami_lab::solvers::FWave::inverse2x2(t_real io_matrix[2][2]){
  //            |a b|    |+d -b|       |a b|
  // inverse of |c d| is |-c +a| / det(|c d|)
  t_real l_a = io_matrix[0][0];
  t_real l_b = io_matrix[0][1];
  t_real l_c = io_matrix[1][0];
  t_real l_d = io_matrix[1][1];
  t_real l_inverseDet = 1.0 / (l_a * l_d - l_b * l_c);
  io_matrix[0][0] = +l_inverseDet * l_d;
  io_matrix[0][1] = -l_inverseDet * l_b;
  io_matrix[1][0] = -l_inverseDet * l_c;
  io_matrix[1][1] = +l_inverseDet * l_a;
}

void tsunami_lab::solvers::FWave::transform2x2(t_real i_matrix[2][2], t_real i_vector[2], t_real o_vector[2]){
  // transforms the two-element i_vector by the 2x2 matrix i_matrix, writes result into o_vector
  t_real v0 = i_vector[0];
  t_real v1 = i_vector[1];
  o_vector[0] = i_matrix[0][0] * v0 + i_matrix[0][1] * v1;
  o_vector[1] = i_matrix[1][0] * v0 + i_matrix[1][1] * v1;
}

void tsunami_lab::solvers::FWave::scaleVector2(t_real i_vector[2], t_real i_scale, t_real o_vector[2] ){
  // scales a vector by i_scale
  o_vector[0] = i_vector[0] * i_scale;
  o_vector[1] = i_vector[1] * i_scale;
}

void tsunami_lab::solvers::FWave::netUpdates( t_real i_hL,
                                              t_real i_hR,
                                              t_real i_huL,
                                              t_real i_huR,
                                              t_real o_netUpdateL[2],
                                              t_real o_netUpdateR[2] ) {
  tsunami_lab::solvers::FWave::netUpdates(i_hL, i_hR, i_huL, i_huR, 0, 0, o_netUpdateL, o_netUpdateR);
}

void tsunami_lab::solvers::FWave::netUpdates( t_real i_hL,
                                              t_real i_hR,
                                              t_real i_huL,
                                              t_real i_huR,
                                              t_real i_bL,
                                              t_real i_bR,
                                              t_real o_netUpdateL[2],
                                              t_real o_netUpdateR[2] ) {
  
  // check if there is (valid) water at all
  if(i_hL + i_hR > 0.0 && i_hL >= 0.0 && i_hR >= 0.0){
    
    t_real l_roeHeight = 0.5 * (i_hL + i_hR);
    
    // sqrt is quite expensive, so reuse it
    t_real l_sqrtHLeft = sqrt(i_hL);
    t_real l_sqrtHRight = sqrt(i_hR);
    
    // velocities
    // tests, because if one side was empty, we would create NaN otherwise
    t_real l_uL = i_hL ? i_huL / i_hL : 0.0;
    t_real l_uR = i_hR ? i_huR / i_hR : 0.0;
    t_real l_roeVelocity = (l_uL * l_sqrtHLeft + l_uR * l_sqrtHRight) / (l_sqrtHLeft + l_sqrtHRight);// velcity mixed by sqrt(h)
    t_real l_gravityTerm = sqrt(m_gravity * l_roeHeight);
    t_real l_bathymetryTerm = m_gravity * (i_bR - i_bL) * l_roeHeight;// [m/s²*m²]
    
    // wave speeds
    t_real l_lambda1 = l_roeVelocity - l_gravityTerm;
    t_real l_lambda2 = l_roeVelocity + l_gravityTerm;
    
    // t_real l_r12[2][2] = { { 1.0, 1.0 }, { l_lambda1, l_lambda2 } };// [1, 1, m/s, m/s]
	// t_real l_inverseDet = 1.0 / (l_lambda2 - l_lambda1);
	t_real l_inverseDet = 0.5 / l_gravityTerm;
    
    t_real l_deltaField[2] = {
      i_huR - i_huL, // [m*m/s]
      i_huR * l_uR - i_huL * l_uL // [m*m/s * m/s]
      + (t_real) 0.5 * m_gravity * ( i_hR * i_hR - i_hL * i_hL ) // [m/s² * m²]
      + l_bathymetryTerm // = 0.5 * gravity * (bR - bL) * (hL + hR)
    };
    
	// inverse matrix * delta field
    t_real l_delta_hL = ( l_lambda2 * l_deltaField[0] - l_deltaField[1]) * l_inverseDet;
    t_real l_delta_hR = (-l_lambda1 * l_deltaField[0] + l_deltaField[1]) * l_inverseDet;
    
    // Z_p = wave1/2 = alpha_p * r_p + bathymetryTerm
    t_real l_delta_huL = l_delta_hL * l_lambda1;
    t_real l_delta_huR = l_delta_hR * l_lambda2;
    
    if(l_lambda1 < 0){// first wave, typically first branch is used
      // to the left
      o_netUpdateL[0] = l_delta_hL;
      o_netUpdateL[1] = l_delta_huL;
      o_netUpdateR[0] = o_netUpdateR[1] = 0;
    } else {
      // to the right
      o_netUpdateR[0] = l_delta_hL;
      o_netUpdateR[1] = l_delta_huL;
      o_netUpdateL[0] = o_netUpdateL[1] = 0;
    }
    
    // In the super-sonic case, both waves can have the same direction.
    // Therefore we need to use "+=" instead of "=".
    if(l_lambda2 < 0){// second wave, typically second branch is used
      // to the left
      o_netUpdateL[0] += l_delta_hR;
      o_netUpdateL[1] += l_delta_huR;
    } else {
      // to the right
      o_netUpdateR[0] += l_delta_hR;
      o_netUpdateR[1] += l_delta_huR;
    }
    
  } else {
    // else there is no water -> no update required
    o_netUpdateL[0] = o_netUpdateL[1] = o_netUpdateR[0] = o_netUpdateR[1] = 0;
  }
}