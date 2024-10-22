``
./reproduce_zetaStudy.sh
``

runs the parametric case study based on the non-dimensional zeta parameter. The user has to start the simulation manually, not automated in the Allrun currently.

The folders CCA and DM1 has the case non-parametrized case setup using the constant contact angle (CCA) boundary condition and the Cox-Voinov dynamic contact angle boundary condition (DM1) which uses the cell-centred velocity field for calculation of Capillary number.

This study requires the compilation of function object [wettedArea](https://github.com/CRC-1194/b01-wetting-benchmark/tree/master/src/functionObjects/wettedArea), the [coxVoinov](https://github.com/CRC-1194/b01-wetting-benchmark/tree/master/src/boundaryConditions/coxVoinov), and the [coxVoinovWithDissipation](https://github.com/CRC-1194/b01-wetting-benchmark/tree/master/src/boundaryConditions/coxVoinovWithDissipation) boundary conditions.
