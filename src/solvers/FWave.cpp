/**
 * @author Antonio Noack
 * @section DESCRIPTION
 * F-Wave solver for the shallow water equations.
 **/
#include "FWave.h"
#include <cmath>
#include <stdio.h>

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
  
  // check if there is (valid) water at all
  if(i_hL + i_hR > 0.0 && i_hL >= 0 && i_hR >= 0){
	
    t_real l_roeHeight = 0.5 * (i_hL + i_hR);
	
	// sqrt is quite expensive, so reuse it
    t_real l_sqrtHLeft = sqrt(i_hL);
	t_real l_sqrtHRight = sqrt(i_hR);
	
	// velocities
	// tests, because if one side was empty, we would create NaN otherwise
	t_real l_uL = i_hL ? i_huL / i_hL : 0.0;
	t_real l_uR = i_hR ? i_huR / i_hR : 0.0;
	t_real l_roeVelocity = (l_uL * l_sqrtHLeft + l_uR * l_sqrtHRight) / (l_sqrtHLeft + l_sqrtHRight);
	t_real l_gravityTerm = sqrt(m_gravity * l_roeHeight);
	
	// wave speeds
	t_real l_lambda1 = l_roeVelocity - l_gravityTerm;
	t_real l_lambda2 = l_roeVelocity + l_gravityTerm;
	
	t_real l_r12[2][2] = { { 1.0, 1.0 }, { l_lambda1, l_lambda2 } };
	inverse2x2(l_r12);
	
	t_real l_deltaField[2] = { i_hR - i_hL, i_huR - i_huL };
	t_real l_alpha[2];
	transform2x2(l_r12, l_deltaField, l_alpha);
	
	// Z_p = wave1/2 = alpha_p * r_p
	t_real l_wave1[2] = { l_alpha[0] * l_lambda1 };
	t_real l_wave2[2] = { l_alpha[1] * l_lambda2 };
	l_wave1[1] = l_wave1[0] * l_lambda1;
	l_wave2[1] = l_wave2[0] * l_lambda2;
	
	if(l_lambda1 < 0){// first wave, typically first branch is used
	  // to the left
	  o_netUpdateL[0] = l_wave1[0];
	  o_netUpdateL[1] = l_wave1[1];
	  o_netUpdateR[0] = o_netUpdateR[1] = 0;
	} else {
	  // to the right
	  o_netUpdateR[0] = l_wave1[0];
	  o_netUpdateR[1] = l_wave1[1];
	  o_netUpdateL[0] = o_netUpdateL[1] = 0;
	}
	
	// In the super-sonic case, both waves can have the same direction.
	// Therefore we need to use "+=" instead of "=".
	if(l_lambda2 < 0){// second wave, typically second branch is used
	  // to the left
	  o_netUpdateL[0] += l_wave2[0];
	  o_netUpdateL[1] += l_wave2[1];
	} else {
	  // to the right
	  o_netUpdateR[0] += l_wave2[0];
	  o_netUpdateR[1] += l_wave2[1];
	}
	
  } else {
	// else there is no water -> no update required
	o_netUpdateL[0] = o_netUpdateL[1] = o_netUpdateR[0] = o_netUpdateR[1] = 0;
  }
  
}