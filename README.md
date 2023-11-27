# FloquetFourierHill.jl
Computational implementation of the Floquet-Fourier-Hill method described in https://www.sciencedirect.com/science/article/pii/S0021999106001665 to compute the spectra of linear operators.

Inputs: 

FFHM(period,numFourierModes,numFloqVals,P,opDim,opPrint,op1,op2=0)

period: Period of the coefficients in the operator

numFourierModes: The number of Fourier modes

numFloquetVals: The number of Floquet values

P: Every p-th diagonal will be filled in the Hill's method matrix (P=2 is the optimal value)

opDim: Dimension of operator (Ex: opDim=1 is a scalar operator, opDim=2 is a 2x2 operator)

opPrint: Boolean option to print out how the code interprets the user input operator

op1: Linear operator

op2 (optional): Optional second linear operator if you have a generalized eigenvalue problem
