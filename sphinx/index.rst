.. Tsunami Simulation documentation master file, created by
   sphinx-quickstart on Sun Oct 24 14:00:52 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Tsunami Simulation's documentation!
==============================================

.. toctree::
  README.rst
  reports/riemann_solver.rst
  reports/finite_volume_discretization.rst
  :maxdepth: 2
  :caption: Contents:


Building the project
--------------------

The project is using SCons, so you can build the code using the command "scons".



Building the docs
------------------

Doxygen: "doxygen Doxyfile"

It will generate the files inside your current directory

Sphinx: "make html" inside sphinx folder. To list all available formats, use "make" without additional arguments inside ./sphinx



Running the actual program
----------------------------

You can run the program via "./build/tsunami_lab <TimeSteps>".
It will generate a set of csv files containing fluid height and impulse, and the coordinates of the water secions for every 25th time step. Therefore, your number of time steps should be at least 25.



Running the tests
---------------------
You can run all tests via "./build/tests". "-l" as an argument lists all tests. Additional commands can be listed with "-?".




Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
