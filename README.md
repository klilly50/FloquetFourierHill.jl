# FloquetFourierHill.jl
Computational implementation of the Floquet-Fourier-Hill method described in https://www.sciencedirect.com/science/article/pii/S0021999106001665 to compute the spectra of linear operators.

**Inputs:**

FFHM(period,numFourierModes,numFloqVals,P,opDim,opPrint,op1,op2=0)

period: Period of the coefficients in the operator

numFourierModes: The number of Fourier modes

numFloquetVals: The number of Floquet values

P: Every p-th diagonal will be filled in the Hill's method matrix (P=2 is the optimal value)

opDim: Dimension of operator (Ex: opDim=1 is a scalar operator, opDim=2 is a 2x2 operator)

opPrint: Boolean option to print out how the code interprets the user input operator

op1: Linear operator

op2 (optional): Optional second linear operator if you have a generalized eigenvalue problem

**Outputs:** 

(evalsArr,evecsArr) in addition to an interactive plot of the eigenvalues in the complex plane

evalsArr: Array of eigenvalues

evecsArr: Array of corresponding eigenvectors

**How the operator should be input:**

1) Each scalar operator should be input as a string
2) All operations should be included as needed (+,-,*,/,^)
3) Operators must be input as D#, where # is the order of the operator
4) Must input in decreasing order of operators
5) The independent variable should always be x
6) Vector problems should be input as a matrix of string operators (Ex: [[S1,S2],[S3,S4]], where each entry is a scalar operator of type string

Example: "-D4-sin(x)\*D2+a\*D1+(1-a^2)sin(x)cos(x)"
