!module containing global variables
!compile with:
!gfortran -c -fdefault-real-8 globalVars.f90

module globalVars
implicit none
real, parameter :: gr = 9.807, pi = 3.141592653589793238462643383279502884197
real, parameter :: m = 1.0, l = 1.0
real Energy
end module
