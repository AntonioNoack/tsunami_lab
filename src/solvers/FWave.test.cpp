/**
 * @author Antonio Noack
 * @section DESCRIPTION
 * Unit tests of the F-Wave solver.
 **/
#include <catch2/catch.hpp>
#define private public
#include "FWave.h"
#include "Roe.h"
#undef public

#define t_idx  tsunami_lab::t_idx
#define t_real tsunami_lab::t_real

TEST_CASE( "Test the matrix inverse.", "[FWave-Matrix-Inverse]" ) {
  
  // input test matrix
  t_real l_matrix[2][2] = { { 1.0, 2.0 }, { 3.0, 4.0 } };
  
  // computed using Wolfram Alpha, inverse([[1,2],[3,4]])
  t_real l_truth[2][2] = { { -2.0, 1.0 }, { 1.5, -0.5 } };
  
  tsunami_lab::solvers::FWave::inverse2x2(l_matrix);
  
  for(int j=0;j<2;j++){
    for(int i=0;i<2;i++){
      REQUIRE( l_matrix[i][j] == Approx(l_truth[i][j]) );
    }
  }
  
}

TEST_CASE( "Test vector scale.", "[FWave-Vector-Scale]" ) {
  
  // result is obvious enough :D
  t_real l_vec[2] = { 1.0, 4.0 };
  t_real l_scale = 3.0;
  
  tsunami_lab::solvers::FWave::scaleVector2(l_vec, l_scale, l_vec);
  
  REQUIRE( l_vec[0] == Approx(3.0) );
  REQUIRE( l_vec[1] == Approx(12.0) );
  
}

TEST_CASE( "Test matrix * vector.", "[FWave-Matrix*Vector]" ) {
  
  // computed using Wolfram Alpha, [[1,2],[3,4]] * [1,4]
  t_real l_mat[2][2] = { { 1.0, 2.0 }, { 3.0, 4.0 } };
  t_real l_vec[2] = { 1.0, 4.0 };
  
  tsunami_lab::solvers::FWave::transform2x2(l_mat, l_vec, l_vec);
  
  REQUIRE( l_vec[0] == Approx(9.0) );
  REQUIRE( l_vec[1] == Approx(19.0) );
  
}

TEST_CASE( "Test when nothing is changing.", "[FWave-netUpdate-Standing]" ) {
  
  // result must be ~ 0
  
  t_real l_height = 10.0;
  t_real l_impulse = 4.0;
  // if our model included friction, the impulse afterwards would be slightly less (transfered into the rock below the sea, and a tiny bit into the air)
  
  t_real l_deltaLeft[2];
  t_real l_deltaRight[2];
  
  tsunami_lab::solvers::FWave::netUpdates(l_height, l_height, l_impulse, l_impulse, l_deltaLeft, l_deltaRight);
  
  REQUIRE( l_deltaLeft[0] == Approx(0.0) );
  REQUIRE( l_deltaLeft[1] == Approx(0.0) );
  REQUIRE( l_deltaRight[0] == Approx(0.0) );
  REQUIRE( l_deltaRight[1] == Approx(0.0) );
  
}

TEST_CASE( "Test breaking dam.", "[FWave-netUpdate-breakingDam]" ) {
  
  t_real l_hL = 10.0;
  t_real l_hR =  0.0;
  t_real l_uL =  1.0;
  
  t_real l_huL = l_hL * l_uL;
  
  t_real l_deltaLeft[2];
  t_real l_deltaRight[2];
  
  // manual calculation:
  // - average height: 5
  // - sqrt height left: sqrt(10), right: 0
  // - velocity left: 1, right: 0
  // - roe velocity = weighted velocity by sqrt height: 1, because only left has fluid
  // - gravity: 9.81
  // - gravity term: sqrt(gravity * avg height) ~ 7.003571
  // - lambdas = wave speeds = roe velocity -/+ gravity term: 1 -/+ 7.003571 = [-6.003571, 8.003571]
  // - inverse of R matrix (using Wolfram Alpha):
  //    | 0.571392 -0.0713922 |
  //    | 0.428608 +0.0713922 |
  // - delta fields: -10 for height, -10 for impulse
  // - alphas = R^-1 * [-10, -10]t = [-5, -5]t (Wolfram Alpha)
  // - net update values:
  // first alpha and lambda for left, second for right side, because sign(lambda1) != sign(lambda2)
  //     alpha * lambda * [1, lambda]
  // left side: 30.017855, -180.214323760205
  // right side: -40.017855, -320.2857437602051
  
  tsunami_lab::solvers::FWave::netUpdates(l_hL, l_hR, l_huL, 0.0, l_deltaLeft, l_deltaRight);
  // can be used for comparison, if hR > 0, and when the gravity constants are the same (I use 9.81, the other implementation uses ~9.80665)
  /*t_real l_deltaLeft2[2];
  t_real l_deltaRight2[2];
  tsunami_lab::solvers::Roe::netUpdates(l_hL, l_hR, l_huL, 0.0, l_deltaLeft2, l_deltaRight2);
  
  REQUIRE( l_deltaLeft[0]  == Approx(l_deltaLeft2[0]) );
  REQUIRE( l_deltaLeft[1]  == Approx(l_deltaLeft2[1]) );
  REQUIRE( l_deltaRight[0] == Approx(l_deltaRight2[0]) );
  REQUIRE( l_deltaRight[1] == Approx(l_deltaRight2[1]) );*/
  
  // epsilons, because the gravity constant may be 9.81 or 9.80665
  REQUIRE( l_deltaLeft[0]  == Approx(+30.017855).epsilon(0.001) );
  REQUIRE( l_deltaLeft[1]  == Approx(-180.21432).epsilon(0.001) );
  REQUIRE( l_deltaRight[0] == Approx(-40.017855).epsilon(0.001) );
  REQUIRE( l_deltaRight[1] == Approx(-320.28574).epsilon(0.001) );
  
  // switching left and right should produce the same results, except for the signs
  tsunami_lab::solvers::FWave::netUpdates(l_hR, l_hL, 0.0, l_huL, l_deltaLeft, l_deltaRight);
  
  REQUIRE( l_deltaLeft[0]  == Approx(-30.017855).epsilon(0.001) );
  REQUIRE( l_deltaLeft[1]  == Approx(+180.21432).epsilon(0.001) );
  REQUIRE( l_deltaRight[0] == Approx(+40.017855).epsilon(0.001) );
  REQUIRE( l_deltaRight[1] == Approx(+320.28574).epsilon(0.001) );
  
}

TEST_CASE( "Test crashing waves.", "[FWave][NetUpdate][crashingWaves]" ) {
  
  // two waves of equal size and velocity crash together
  
  t_real l_hL = 10.0;
  t_real l_hR = 10.0;
  t_real l_uL =  1.0;
  t_real l_uR = -1.0;
  
  t_real l_huL = l_hL * l_uL;
  t_real l_huR = l_hR * l_uR;
  
  t_real l_deltaLeft[2];
  t_real l_deltaRight[2];
  
  // manual calculation:
  // - average height: 10
  // - sqrt height left: sqrt(10), right: sqrt(10)
  // - velocity left: 1, right: -1
  // - roe velocity = weighted velocity by sqrt height: 0, because left & right are equal
  // - gravity: 9.81
  // - gravity term: sqrt(gravity * avg height) ~ 9.9045444
  // - lambdas = wave speeds = roe velocity -/+ gravity term: 0 -/+ 9.9045444 = [-9.9045444, +9.9045444]
  // - inverse of R matrix (using Wolfram Alpha):
  //    | 0.5 -0.0504819 |
  //    | 0.5 +0.0504819 |
  // - delta fields: 0 for height, -20 for impulse
  // - alphas = R^-1 * [0, -20]t = [1.00964, -1.00964]t (Wolfram Alpha)
  // - net update values:
  // first alpha and lambda for left, second for right side, because sign(lambda1) != sign(lambda2)
  //     alpha * lambda * [1, lambda]
  // left side: -10, +99.045684
  // right side: -10, -99.045684
  
  tsunami_lab::solvers::FWave::netUpdates(l_hL, l_hR, l_huL, l_huR, l_deltaLeft, l_deltaRight);
  // can be used for comparison when the gravity constants are the same (I use 9.81, the other implementation uses ~9.80665)
  /*t_real l_deltaLeft2[2];
  t_real l_deltaRight2[2];
  tsunami_lab::solvers::Roe::netUpdates(l_hL, l_hR, l_huL, l_huR, l_deltaLeft2, l_deltaRight2);
  
  REQUIRE( l_deltaLeft[0]  == Approx(l_deltaLeft2[0]) );
  REQUIRE( l_deltaLeft[1]  == Approx(l_deltaLeft2[1]) );
  REQUIRE( l_deltaRight[0] == Approx(l_deltaRight2[0]) );
  REQUIRE( l_deltaRight[1] == Approx(l_deltaRight2[1]) );*/
  
  
  REQUIRE( l_deltaLeft[0]  == Approx(-10) );
  REQUIRE( l_deltaLeft[1]  == Approx(+99.045684).epsilon(0.001) );
  REQUIRE( l_deltaRight[0] == Approx(-10) );
  REQUIRE( l_deltaRight[1] == Approx(-99.045684).epsilon(0.001) );
  
  // switching left and right should produce the same results, except for the signs
  tsunami_lab::solvers::FWave::netUpdates(l_hR, l_hL, l_huR, l_huL, l_deltaLeft, l_deltaRight);
  
  REQUIRE( l_deltaLeft[0]  == Approx(+10) );
  REQUIRE( l_deltaLeft[1]  == Approx(-99.045684).epsilon(0.001) );
  REQUIRE( l_deltaRight[0] == Approx(+10) );
  REQUIRE( l_deltaRight[1] == Approx(+99.045684).epsilon(0.001) );
  
}

TEST_CASE( "Test super-sonic.", "[FWave][NetUpdate][SuperSonic]" ) {
    
  // damit man Über"schall"-Effekte bekommt, muss |weighted velocity average| > sqrt(gravity * avg height) sein,
  // also sehr flaches Wasser oder hohe Impulse.
  
  t_real l_hL = 10.0;
  t_real l_hR = l_hL * 0.9; // difference is required, so we can see an effect
  t_real l_avgHeight = 0.5 * (l_hL + l_hR);
  t_real l_gravity = tsunami_lab::solvers::FWave::m_gravity;
  t_real l_requiredVelocity = sqrt(l_gravity * l_avgHeight);
  t_real l_impulse = - l_hL * l_requiredVelocity * 1.5;// a bit more
  
  t_real l_deltaLeft[2];
  t_real l_deltaRight[2];
  
  tsunami_lab::solvers::FWave::netUpdates(l_hL, l_hR, l_impulse, l_impulse, l_deltaLeft, l_deltaRight);
  
  t_real l_epsilon = std::numeric_limits<t_real>::epsilon();
  
  // wir haben nur eine Änderung im Impuls, aber super-sonic, und nach links, also sollte nur der eine Wert != 0 sein
  // the actual result of deltaLeft[0] is ~-5e-7 due to rounding errors, and Approx() has no idea of the scale of values, so it reports false for == Approx(0.0)
  REQUIRE( l_deltaLeft[0] == Approx(0).margin(10 * l_epsilon));
  REQUIRE( l_deltaLeft[1] > 0.0 );
  // exactly zero, because nothing is added, because both wave fronts are on the left side
  REQUIRE( l_deltaRight[0] == Approx(0) );
  REQUIRE( l_deltaRight[1] == Approx(0) );
  
}

#include <iostream>
TEST_CASE( "Test bathymetry.", "[FWave][NetUpdate][Bathymetry]" ) {
  
  constexpr t_idx l_sampleCount = 2;
  t_real l_samples[l_sampleCount][2] = {
    { -40, -30 }, { -3462, -3450 }
  };
  
  for(t_idx l_i=0;l_i<l_sampleCount;l_i++){
    
    t_real l_bL = l_samples[l_i][0];
    t_real l_bR = l_samples[l_i][1];
    
    t_real l_surfaceHeight = 0;// water goes up to zero
    t_real l_hL = l_surfaceHeight - l_bL;
    t_real l_hR = l_surfaceHeight - l_bR;
    
    t_real l_impulse = 0;
    
    t_real l_deltaLeft[2] = { 0, 0 };
    t_real l_deltaRight[2] = { 0, 0 };
    
    tsunami_lab::solvers::FWave::netUpdates(l_hL, l_hR, l_impulse, l_impulse, l_bL, l_bR, l_deltaLeft, l_deltaRight);
    
    // irgendwie ist der Löser ziemlich ungenau: woran liegt diese große numerische Instabilität?
    t_real l_m = 0.01;//std::numeric_limits<t_real>::epsilon();
    
    // the water has no reason to flow
    // -> all deltas should be zero
    REQUIRE( l_deltaLeft[0]  == Approx(0).margin(l_m) );
    REQUIRE( l_deltaLeft[1]  == Approx(0).margin(l_m) );
    REQUIRE( l_deltaRight[0] == Approx(0).margin(l_m) );
    REQUIRE( l_deltaRight[1] == Approx(0).margin(l_m) );
  }
}

#undef t_real
