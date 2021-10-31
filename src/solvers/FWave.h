
#ifndef TSUNAMI_LAB_SOLVERS_FWAVE
#define TSUNAMI_LAB_SOLVERS_FWAVE

#include "../constants.h"

namespace tsunami_lab {
  namespace solvers {
    class FWave;
  }
}

/**
 * @author Antonio Noack
 * @section DESCRIPTION
 * F-Wave solver for the one-dimensional shallow water equations.
 **/
class tsunami_lab::solvers::FWave {
	
    /**
     * Inverses a 2x2 matrix.
     *
     * @param io_matrix the matrix, which is inversed in-place.
     * If we need this function more often in the future, it may be better to place it in a linear-algebra-related file.
     **/
    static void inverse2x2( t_real io_matrix[2][2] );

    /**
     * Multiplies a 2x2 matrix with a 2 element vector
     *
     * @param i_matrix the matrix.
     * @param i_vector the input vector.
     * @param o_vector the output vector.
     * If we need this function more often in the future, it may be better to place it in a linear-algebra-related file.
     **/
    static void transform2x2( t_real i_matrix[2][2],
                              t_real i_vector[2], 
                              t_real o_vector[2] );

    /**
     * Multiplies a 2 element vector with a scalar.
     * @param i_vector input vector.
     * @param i_scale scalar.
     * @param o_vector result vector.
     **/
    static void scaleVector2( t_real i_vector[2],
                              t_real i_scale,
                              t_real o_vector[2] );

  public:
  
    //! gravity, in m/s²
	// all tests seem to depend on this constant, so I switch from 9.81 to 9.80665; 3.131557121*3.131557121;
    static t_real constexpr m_gravity = 9.80665;
	
    /**
     * Computes the net-updates.
     *
     * @param i_hL height of the left side.
     * @param i_hR height of the right side.
     * @param i_huL momentum of the left side.
     * @param i_huR momentum of the right side.
     * @param o_netUpdateL will be set to the net-updates for the left side; 0: height, 1: momentum.
     * @param o_netUpdateR will be set to the net-updates for the right side; 0: height, 1: momentum.
     **/
    static void netUpdates( t_real i_hL,
                            t_real i_hR,
                            t_real i_huL,
                            t_real i_huR,
                            t_real o_netUpdateL[2],
                            t_real o_netUpdateR[2] );
};

#endif // #ifndef TSUNAMI_LAB_SOLVERS_FWAVE
