``
./reproduce_dynamicsConvergence.sh
``

runs the parametric case study. The user has to start the simulation manually, not automated in the Allrun currently.

Mesh convergence study for partial wetting dynamics validation case study.

This study requires the compilation of function object [wettedArea](https://github.com/CRC-1194/b01-wetting-benchmark/tree/master/src/functionObjects/wettedArea), the [coxVoinov](https://github.com/CRC-1194/b01-wetting-benchmark/tree/master/src/boundaryConditions/coxVoinov), and the [coxVoinovWithDissipation](https://github.com/CRC-1194/b01-wetting-benchmark/tree/master/src/boundaryConditions/coxVoinovWithDissipation) boundary conditions.

