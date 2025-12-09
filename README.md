# NBody
N Body problem

## Building

The optional executable for visualisation has dependencies that are provided as git submodules, so make sure that you clone the repository with the `--recurse-submodules` option or use
```bash
git submodule update --init --recursive
```
on an already cloned repository.

There are two CMake options available:
 * `BUILD_SIMULATION` (default: `ON`) Build the simulation executables
 * `BUILD_VISUALISATION` (default: `OFF`) Build the visualisation executable

Then the specified executables can be built using the standard CMake procedure.
