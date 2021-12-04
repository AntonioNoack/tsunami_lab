##
# @author Alexander Breuer (alex.breuer AT uni-jena.de), Antonio Noack
# 
# @section LICENSE
# Copyright 2020, Friedrich Schiller University Jena
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# @section DESCRIPTION
# Entry-point for builds.
##
import SCons
import os
import inspect

print( '###################################' )
print( '### Tsunami Lab                 ###' )
print( '###                             ###' )
print( '### http://scalable.uni-jena.de ###' )
print( '###################################' )
print()
print('runnning build script')

# configuration
vars = Variables()

vars.AddVariables(
  EnumVariable( 'mode',
                'compile modes, option \'san\' enables address and undefined behavior sanitizers',
                'release',
                allowed_values=('release', 'debug', 'release+san', 'debug+san' )
              )
)

# exit in the case of unknown variables
if vars.UnknownVariables():
  print( "build configuration corrupted, don't know what to do with: " + str(vars.UnknownVariables().keys()) )
  exit(1)

# create environment
env = Environment( variables = vars )

# generate help message
Help( vars.GenerateHelpText( env ) )

# add default flags
env.Append( CXXFLAGS = [ '-std=c++11',
                         '-Wall',
                         '-Wextra',
                         '-Wpedantic',
                         '-Werror' ] )

# set optimization mode
if 'debug' in env['mode']:
  env.Append( CXXFLAGS = [ '-g',
                           '-O0' ] )
else:
  env.Append( CXXFLAGS = [ '-O3' ] )

# add sanitizers
if 'san' in  env['mode']:
  env.Append( CXXFLAGS =  [ '-g',
                            '-fsanitize=float-divide-by-zero',
                            '-fsanitize=bounds',
                            '-fsanitize=address',
                            '-fsanitize=undefined',
                            '-fno-omit-frame-pointer' ] )
  env.Append( LINKFLAGS = [ '-g',
                            '-fsanitize=address',
                            '-fsanitize=undefined' ] )

env.Append( CXXFLAGS = [ '-g' ] )

# add Catch2
env.Append( CXXFLAGS = [ '-Isubmodules/Catch2/single_include' ] )

# add netCDF & YamlCpp
conf = Configure(env)
externalLibs = []
# packages: netcdf-bin, libnetcdf-dev
# the task says, that zlib and hdf5 are required as packages, but SCons tells me it cannot find them, but builds fine anyways
libs = ['netcdf', 'yaml-cpp'] # ,'zlib','hdf5'
for lib in libs:
  if conf.CheckLib(lib):
    externalLibs.append(lib)
  else:
    print(f'libary {lib} was not found!')

env.Append( CXXFLAGS = [ '-fopenmp' ] )
env.Append( LINKFLAGS = [ '-fopenmp' ] )

# get source files
VariantDir( variant_dir = 'build/src', src_dir     = 'src' )

env.sources = []
env.tests = []
env.tsunami1d = []

#src_file_path = inspect.getfile(lambda: None)
#yamlcpp_src = os.path.dirname(src_file_path) + "/submodules/YamlCpp2/src"
#for file in os.listdir(yamlcpp_src):
#  src_file = yamlcpp_src + "/" + file
#  if os.path.isfile(src_file):
#    if src_file.endswith(".cpp"):
#      env.sources.append(src_file)
#  else:
#    pass

Export('env')
SConscript( 'build/src/SConscript' )
Import('env')

env.Program( target = 'build/tsunami_lab',
             source = env.sources + env.standalone,
             LIBS = externalLibs)

env.Program( target = 'build/tests',
             source = env.sources + env.tests,
             LIBS = externalLibs)