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

## Sequential algorithms

For both the Basic and Reduced sequential algorithms, a **struct** including mass, position, and velocity was defined for the particles to increase data locality. Since the force matrix requires initialization to zero at every iteration, we chose to keep it separate from the struct so that all its elements are stored in contiguous memory cells.

As the Euler method does not enforce **energy conservation**, a total_energy function was implemented to verify that the system's total energy variation remains negligible. Tests confirmed the model's validity, showing that as dt approaches zero, the energy variation decreases as expected.

The **computing time** for the Basic version was verified to be significantly higher than the Reduced one. This is due to the Reduced version halving the number of arithmetic operations by exploiting Newton's third law, although both algorithms technically maintain $O(N^2)$ complexity. Both solvers share a similar implementation structure to facilitate performance comparison.

Extensive **stress tests** were performed in both 2D and 3D to ensure code robustness, including scenarios with high particle counts, star-planet systems, clusters with outliers, near-frontal collisions, small gravitational constants, and extreme mass/speed values. Both solvers handled all cases successfully. Furthermore, comparing the outputs of the two solvers across various inputs consistently showed identical results.

## MPI implementation

MPI implementation is done for both Basic and Reduced versions, following the algorithms described in the Pacheco book[1].

In the Basic version, the MPI root process (the one with rank 0) scatter the particle information in a blocked distribution. Moreover, the masses and initial positions are broadcasted to each MPI process. At each simulation step, the velocities and positions for local particles are updated, and then each process gathers all new position information using `MPI_Allgatherv`. In short, in the Basic version, after each simulation step, the MPI processes share the new positions with all the other processes.

In the Reduced version, a ringpass algorithm is used. Processes are connected in a ring: process `i` is connected to `i-1` and `i+1` (and largest ranked process is connected to process with rank 0). Initially, the processes know only the particle data for their local particles. At each simulation step, the algorithm performs P-1 phases (where P is the number of MPI processes). At each phase, each MPI process sends data to a lower-ranked neighbor, and receives from a higher-ranked data. With this new data, the process computes all inter-particle forces, and adds to a local force array, while the opposite forces are subtracted from the passing buffer of the forces. Lastly, one more exchange is done and then each process can sum all the forces that they got during this process. 

The Reduced version is faster as it avoids all-to-all gathering at every simulation step.

## Execution time comparison

The comparison made with generated random file using `utils::generateRandomToFile<3>("test-timer.in", 1000, 5000, 0.01, 42);`, i.e. 1000 particles and 5000 steps.

| **Algorithm** | Execution time (seconds) |
|--------------|----------|
| **Serial Basic** | 29.152 |
| **Serial Reduced** | 20.187 |
| **MPI Basic** (4 processes) | 10.639 |
| **MPI Reduced** (4 processes) | 3.492 |


## References

[1] Chapter 7 of Pacheco Book
