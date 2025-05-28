# Welcome to MSA

## About MSA
MSA is short for *Matrix Structural Analysis*, and referes to the analysis of skeletal structures using the direct stiffness method.

Python package for the direct stiffness method of matrix analysis of skeletal structures. At present a subpackage `pf` for analysis of plane frames has been implemented. Code is based on the flow charts in *Weaver and Gere* and an example is taken from *Hall and Kabaila*.

## Limitations
The package has the following limitations

1. Analyses only plane frame structures, with only plane frame elements.
2. Non-zero boundary conditions are currently not allowed. Only known zerodisplacement boundary condition is implemented.
3. Member loads must be applied as equivalent member end-forces. Equivalent member end-forces must be calculated manually - as the negative values from the reactions produces in the corresponding fixed-end member.
4. Only member end-forces are calculated. Forces at points within the member are not calculated.

## References
1. Weaver, W. and Gere, J.M., *Matrix Analysis of Framed Structures*, 2ed., CBS Publishers, New Delhi, 1986
2. Hall, A. and Kabaila, A.P., *Basic Concepts of Structural Analysis*, Pitman Publishing, London, 1977
