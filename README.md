WARP
=======================
Ryan M. Bergmann, 2014.

WARP is a GPU-powered library for continuous energy Monte Carlo neutron transport in general 3D geometries.  The host-side code has been written in C++, CUDA, and Python.  The C++ API is exposed to Python via SWIG.  The Python C API is used within WARP to load ACE-formatted cross section libraries via the PyNE package.  Geometry routines are implemented using the NVIDIA OptiX ray tracing framework.  Currently, WARP is only able to run on a single NVIDIA GPU.

If the name "WARP" were an acronym, it would stand for "weaving all the random particles," with the word "weaving" referring to the lock-step way in which "all the random particles," i.e. the neutrons, are sorted into coherent bundles and transported.  Using the word "warp" is also a nod to NVIDIA's terminology for a group of 32 concurrent threads. 

Please refer to the usage guide for examples of how to use WARP, and to the doxygen documentation for detail descriptions of the API.

### Related Publications:

[Performance and accuracy of criticality calculations performed using WARP – A framework for continuous energy Monte Carlo neutron transport in general 3D geometries on GPU](https://doi.org/10.1016/j.anucene.2017.01.027)

[Algorithmic choices in WARP – A framework for continuous energy Monte Carlo neutron transport in general 3D geometries on GPUs](https://doi.org/10.1016/j.anucene.2014.10.039)


### Miscellaneous

Logo provided by designer extrodinaire Daniel Castro, http://www.danieldanieldaniel.com/

