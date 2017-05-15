Parallel Benchmark Cases

This directory contains cases that are used for measuring the parallel performance
of MFIX. All the cases have identical geometry, discretization, and physical models.
The total number of computational cells is a power of 2 (e.g., 2^18), to facilitate
the even distribution of cells on multiple processors. The initial and boundary 
conditions are similar wherever possible. An increasing degree of complexity then
introduced into the cases, ultimately to facilitate a gasifier simulation:

Case A: 
Solves only the hydrodynamic equations including the granular energy equation;

Case B:
Solves gas species equations to simulate catalytic decomposition of ozone;

Case C:
Solves gas and solids species equations and energy equations to simulate char combustion;(yet to be added)

Case D:
Solves all of the above equations and includes the complex gasification 
chemistry. (yet to be added)

