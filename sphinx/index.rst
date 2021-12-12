.. Tsunami Simulation documentation master file

Welcome to Tsunami Simulation's documentation!
==============================================

.. toctree::
  reports/riemann_solver.rst
  reports/finite_volume_discretization.rst
  reports/bathymetry.rst
  reports/two_dimensions.rst
  reports/large_data_io.rst
  reports/tsunami_simulations.rst
  reports/coarse_output.rst
  reports/optimization.rst
  :maxdepth: 2
  :caption: Contents:


Building the project
--------------------

Install the packages "netcdf-bin" and "libnetcdf-dev" for NetCDF.

Install YAML-Cpp.

The project is using SCons, so you can build the code using the command "scons".



Building the docs
-----------------

Doxygen: "doxygen Doxyfile"

It will generate the files inside your current directory

Sphinx: "make html" inside sphinx folder. To list all available formats, use "make" without additional arguments inside ./sphinx



Running the actual program
--------------------------

You can run the program via "./build/tsunami_lab <ConfigurationFile>".
There are a lot of sample configurations in the config folder, you can use. From these samples, you can create your own scenarios.



Running the tests
-----------------

You can run all tests via "./build/tests". "-l" as an argument lists all tests. Additional commands can be listed with "-?".




Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
