# NBody
N Body problem

This project contains multiple implementations of the N-body problem for comparison between serial, parallel, and reduced versions of the algorithm. 

Using the generated data, a visualization tool is provided using OpenGL and GLFW to show the particles. It reads the generated files (for easier cluster offloading in the future) and can either play the simulation loop or load the next time frame one step at a time.

![visualization](images/vis.png)

## Building

The optional executable for visualisation has dependencies that are provided as git submodules, so make sure that you clone the repository with the `--recurse-submodules` option or use
```bash
git submodule update --init --recursive
```
on an already cloned repository.

There are two CMake options available:
 * `BUILD_SIMULATION` (default: `ON`) Build the simulation executables
 * `BUILD_VISUALISATION` (default: `OFF`) Build the visualisation executable

### Simulation

First, make sure that you have an environment with MPI installed. 
For instance, open the container used in the laboratories and load the `gcc-glibc` module.
Then, create and enter a build directory and run `cmake ..` to configure the build system. 
Finally, running `make` should build all simulation-related targets.

### Visualisation

The container used in the laboratories does not support OpenGL. 
Therefore, visualisation targets are set up to be optional.
To build them, first ensure that your environment supports OpenGL (i.e. you are not inside the laboratory container) and that you have fetched all git submodules.
Then, create a new build directory (different from the one configured for simulation) and enter it.
To configure the build environment only for the visualisation target run
```bash
cmake .. -DBUILD_SIMULATION=OFF -DBUILD_VISUALISATION=ON
```
and then to build the target, simply run `make`.
