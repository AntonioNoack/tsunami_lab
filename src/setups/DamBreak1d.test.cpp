/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de), Antonio Noack
 * 
 * @section LICENSE
 * Copyright 2020, Friedrich Schiller University Jena
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Tests the dam break setup.
 **/
#include <catch2/catch.hpp>
#include "DamBreak1d.h"

TEST_CASE( "Test the one-dimensional dam break setup.", "[DamBreak1d]" ) {
  tsunami_lab::setups::DamBreak1d l_damBreak( 25, 55, 3, 0 );

  // left side
  REQUIRE( l_damBreak.getHeight( 2, 0 ) == 25 );

  REQUIRE( l_damBreak.getMomentumX( 2, 0 ) == 0 );

  REQUIRE( l_damBreak.getMomentumY( 2, 0 ) == 0 );

  REQUIRE( l_damBreak.getHeight( 2, 5 ) == 25 );

  REQUIRE( l_damBreak.getMomentumX( 2, 5 ) == 0 );

  REQUIRE( l_damBreak.getMomentumY( 2, 2 ) == 0 );

  // right side
  REQUIRE( l_damBreak.getHeight( 4, 0 ) == 55 );

  REQUIRE( l_damBreak.getMomentumX( 4, 0 ) == 0 );

  REQUIRE( l_damBreak.getMomentumY( 4, 0 ) == 0 );

  REQUIRE( l_damBreak.getHeight( 4, 5 ) == 55 );

  REQUIRE( l_damBreak.getMomentumX( 4, 5 ) == 0 );

  REQUIRE( l_damBreak.getMomentumY( 4, 2 ) == 0 );  
}
