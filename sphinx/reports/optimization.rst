
Optimization
============


Connecting to a Node
--------------------

I created a job-node with the command scalloc --partition=s_hadoop, and got assigned node 274. Unfortunately, it looks like the execution nodes don't have my favourite "console text editor" Nano installed.
Vim is installed there, so I might have to get used to it, or log in with two ssh shells.


Compiling on the ARA-cluster
----------------------------

With the installing-instructions, everything worked fine.

When compiling the tsunami simulation, the only issue was that it didn't know "pragma omp simd". I loaded the OpenMPI and GCC module for version 10.0.2.

When I then tried to run the compiled tsunami simulation, it reported that libnetcdf.so.7 was not found. With the help of `a forum post <https://code.mpimet.mpg.de/boards/2/topics/939>`_, I found that I just needed to add `LD_LIBRARY_PATH=./ara/software/lib:$LD_LIBRARY_PATH` to the environment variables, where ara/software/ is the path where I installed NetCDF, ZLib and HDF5. And then I copied `ara/software/lib/libnetcdf.so.18` to `ara/software/lib/libnetcdf.so.7`.

All existing tests passed.


Optimizations
-------------

To evaluate optimizations, I first added timers to the functions initWithSetup(), timeStep() and computeMaxTimestep().

Then I inlined the matrix multiplication, assumed the water height to be positive, and made the solver be decidable via preprocessor definition instead of a class member variable.

This brought the computation time for the 250m grid on 12 threads with my Ryzen 5 2600 from down 0.33s/timestep to 0.31s/timestep (~6% improvement).

Then I tried whether the flag `-g` has an influence on the performance, but it didn't. My idea was that the compiler might inline less functions, when this flag is enabled.


Scarce Memory Mode
------------------

The mode for when the memory is scarce iterates in the inner loop over the columns, while the data layout is row-major. This is pretty sub-optimal: for the 250m/cell grid, a single timestep needs 1.16s instead of 0.31s. The output additionally writes how much longer the y-part needed then the x-part: it's about 7.85 times slower than the iteration on the x-axis. Usually, they have the same performance. The total speedup/slow-down is roughly half of that, because the x-part uses cache-optimized accesses in both cases.

Additional slow-down comes from the slow functions to write line by line with NetCDF.


Comparing Runtimes
------------------

12 Threads, Ryzen 5 2600, 250m/cell Tohoku: 0.33s/timestep, with VTune
0.31s without VTune as well

6 Threads, 0.34s
3 Threads, 0.60s
1 Thread,  1.62s

Skylake processor with 18 cores, dual socket
72 Threads, 0.28-0.31s/timestep
36 Threads, 0.38-0.39s/timestep
18 Threads, 0.47-0.48s/timestep
8 Threads, 0.65-0.66s/timestep
4 Threads, 1.14-1.15s/timestep
2 Threads, 1.88-1.90s/timestep (0.97x slower for y)
1 Thread,  2.17-2.21s/timestep (0.92x slower for y)

Therefore, currently my six core processor is only a little slower than the 18 core skylake processor in when computing the 250m/cell grid.
Part of the difference is probably the core clock: my processor runs on 3.4GHz while the Xeon processor runs on 2.3GHz.
The Skylake node has the major advantage that it has much more memory (192GB instead of 32GB), and therefore can compute large tsunami simulations such as the 50m/cell grid with much higher speed.


GNU Compiler vs GCC
-------------------
For the comparison between compilers, I tested on 36 Threads.
Sadly, the SCons didn't find the path to ipcp itself, so I had to use long commands.
Additionally, GCC 10.0.2 uses a confusing path: `/cluster/spack/opt/spack/linux-centos7-broadwell/gcc-8.1.0/gcc-10.2.0-ru4xdhhkxnma5i727b7njtnjoh6kff3s/bin/g++`. I compiled with -O3 and with -march=skylake-avx512.

The Intel compiler gave out the following warnings: "Inlining inhibited by limit max-size" and "Inlining inhibited by limit max-total-size".

GCC 10.0.2, with AVX 512: 0.334s/timestep, 0.0082s/max
Intel 19.1.2.254, with AVX 512: 0.312s/timestep, 0.0076s/max
GCC 9.3.0, with AVX 512: 0.334s/timestep, 0.0082s/max


-Ofast
----------

When researching -Ofast, I found `an interesting StackOverflow article <https://stackoverflow.com/questions/45685487/what-does-clangs-ofast-option-do-in-practical-terms-especially-for-any-diffe>` about it.
There are flags, which disable handling of special cases like NaNs and infinities, or the sign of a zero. In our tsunami simulation, NaNs and infinities should not happen, and the sign of a zero plays no role.

Intel with AVX 512 and Ofast: 0.319s/timestep, 0.0077s/max


Installing Intel VTune
----------------------

I wanted to use VTune without the massive latency from ssh.
To install it, I followed `a tutorial to install intel-basekit <https://www.intel.com/content/www/us/en/develop/documentation/installation-guide-for-intel-oneapi-toolkits-linux/top/installation/install-using-package-managers/apt.html>`_. I added their key and then installed `intel-oneapi-vtune`, because I only wanted VTune (1.5GB), not their other stuff (15GB total). Unfortunately, the installation does not add their exectutable to the path automatically. I found the executable in `/opt/intel/oneapi/vtune/2022.0.0/bin64`.

.. figure:: w8_vtune_installed.png

VTune might need other packages (though they usually are installed automatically), but currently it says that it does not know my processor. It's a Ryzen 5 2600, so not a new one, just not an Intel one :/.

.. figure:: w8_vtune_cpu_unknown.png

Their constant tutorial-tries are pretty annoying.


Is FWave::netUpdates inlined/vectorized?
----------------------------------------

FWave::netUpdates appears normally in the Top-down Tree view, and when looking at the disassembled code of the function in VTune, it doesn't seem to be vectorized either. This is understandable, because the function has a lot of branches currently.


