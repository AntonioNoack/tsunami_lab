
Two-Dimensional Solver & Stations
=================================

First I added the two dimensionality. I also split the apply-flux-part into a separate function, because it's needed twice now.
For the tests, I just repurposed the 1d tests, because most things stay the same and can be tested visually.


Merging the main()s
-------------------

As recommended, I then replaced my three main() functions by a single one, which loads the runtime configuration from a YAML-like file format.
First, I implemented the reader myself, because often it's easier than to learn a new library. It comes with the limitation that the current format ignores indentation.

Because of that, stations (as requested in the next section) were specified in a CSV file, and names were limited to numeric values (because of the existing CSV reading method).

This is a nice improvement, because testing and reproducibility has simplified :).
A slight issue is that I'm often tempted to change parameters inside a file instead of creating a copy of that file. This nullifies the reproducibility-advantage.


For usability, I added an additional time limit. The defaults for both are infinity (or int::max), so the user can decide whether he wants a timestep limit or a time limit.
Printing results as CSV now happens at regular simulation time intervals instead of regular step intervals.


3d Debugging in ParaView 
------------------------

The 2d tsunami-simulations can be visualized in ParaView using the action "Table to Points" (Extractors/Search.../Enter "table to points"). You need to have your dataset selected to apply it. You then have to make it visible and define the x,y,z components of these points. Then apply the changes.
You can add a second layer by applying this action again on the original source, and by changing the z-component.
However, multiple layers often are hard to separate, so it's not that useful.

I like that it is possible, but I don't like the Arcball controls (rotation matrices are accumulated). I'd much rather have Orbit controls (pitch & yaw are accumulated). I haven't found a way to change it. (Why do people need the "roll" axis anyways?)

I also found the option "Table to Structured Grid", but haven't got it to work yet.


Circular dam break
------------------

The circular dam break works, but it is not perfectly symmetrical. When the wave is progressed pretty far, the diagonal edges become more and more straight. The following image was created on a 200x200 grid with dam radius of 40, and was taken after 6 seconds.

The configuration for this is simulation is defined in config/dambreak2d.yaml, so you can reproduce this sample by calling ./build/tsunami_lab config/dambreak2d.yaml.

.. figure:: w4_circle2d_t6.0_200r40.png

Obstacle in 2d fluid simulation
-------------------------------

First I experimented with an obstacle, which was dry. It reflected the incoming wave as expected. The config for that example is config/dambreak2d5050obstacleDry.yaml.
In the following images, the rendered height is the surface, the color is the water depth. The shallow section is blue, and the rest is orange-ish.

.. figure:: w4_obstacle_reflection.png

When I tried to experiment with the wave flowing over a very shallow piece of land, the solver had issues (kind of as expected), and the timesteps were really small, so it needed a long time to solve it. The higher the water, the easier and faster it was solvable. With a height of 10cm on top of the obstacle, I even got a NaN in the result. I tracked it down, and found that the issue was in the FWave solver when combining the waves: it was beyond super-sonic speeds, and both waves had a computed height/impulse update of +/-Infinity. Adding +Infinity and -Infinity resulted in NaN as output.
Increasing the precision to double wasn't enough, so I wrote as a "fix", that every update, which results in NaN, is reset to the original value (time stops).
The config for this test can be found in config/dambreak2d5050obstacleShallowCritical.yaml.

When looking at the result of this test, you see a few peaks (orange color) in the blue section. These high peaks were obtained with 100k steps instead of 10k in the file. 100k steps currently take quite a while to process.

.. figure:: w4_extremely_shallow.png


When a less extreme shallowness is used, the water behaves properly. Interesting is the triangular wave form on the shallow section.

.. figure:: w4_shallow.png


Stations
========

For properly loading the configuration files, I didn't want to switch to XML files, because YAML files look user friendlier. So I tried to use `Yaml-Cpp <https://github.com/jbeder/yaml-cpp>`_. Unfortunately, it is a CMake project and ours is a SCons project. When I added the include folder to SCons, it was missing a mysterious "dll.h", which was not part of the repository. So I executed CMake and found out that it generates the dll.h dynamically as some form of configuration.

To fix that, I forked the original project, added the missing dll.h, and changed the submodule to `my fork <https://github.com/AntonioNoack/yaml-cpp>`_.
First, I have loaded the stations from a separate CSV file, but thanks to this library, I could relatively easily load them as a YAML list.
Each station is a list entry with the properties name, x, and y.

Comparison 1d vs 2d solver
--------------------------

To compare the 1d and 2d solver, I created the configurations config/stations1d.yaml and config/stations2d.yaml. The results are similar, but not the same. The 2d sample seems to oscillate longer. It calms at t=14s, while the 1d calms already at 8s. Anyways, this doesn't mean that my solver is incorrect. The waves flow out as expected, no side seems different from another, and neither 1d nor 2d solutions are correct, and result in a water height of 5m at the end. The 1d sample just is a little closer to the true result (4.93m instead of 4.89m).

I wanted to take a closer look at the issue, so I run the test on a much larger scene, 510x510m in size, and there the issues are obvious. The solver doesn't have circular symmetry, only every 90 degrees.

.. raw:: html

    <iframe width="696" height="573" src="http://www.youtube.com/embed/KctdpZPAciY?rel=0" frameborder="0" allowfullscreen></iframe></br></br>

Convergence
-----------

In the convergence test, I use the cell sizes 50, 25, 10, 5, 2 and 1. Indeed does the curve converge roughly at a cell size of 5m, and 2m is pretty much the same as with 1m after enough time. The following video shows the progression of the station records from 50m cell size down to 1m cell size.

.. raw:: html

    <iframe width="696" height="573" src="http://www.youtube.com/embed/XzWFdh4GkbM?rel=0" frameborder="0" allowfullscreen></iframe></br></br>

