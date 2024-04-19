2D Unstructured Finite Volume Solver

This project provides a Fortran implementation for solving 2D unstructured finite volume problems based on the following research paper https://hal.science/hal-03157500v1. It includes functionalities for mesh generation, boundary condition implementation, flux computation, and error analysis, among other solver routines.

*** Features ***

+ Mesh Generation: Prepare and manipulate unstructured grids for simulation.

+ Boundary Conditions: Apply various boundary conditions to simulate physical scenarios accurately.

+ Flux Calculation: Compute fluxes across cell faces to solve conservation laws.

+ Error Computation: Assess the accuracy of simulation results.

+ Solver Routines: Implement numerical methods to solve partial differential equations on unstructured meshes.

*** Getting Started ***

To use this solver, clone the repository and compile the source code using the provided makefiles. Ensure you have a Fortran compiler installed.

git clone https://github.com/youzaina001/2D_Unstructured_Finite_Volume.git
cd 2D_Unstructured_Finite_Volume
make

*** Contributing ***

Contributions to improve the solver or add new features are welcome. Please feel free to fork the repository and submit pull requests.
