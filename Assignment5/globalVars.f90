!module containing global variables
!compile with:
!gfortran -c -fdefault-real-8 globalVars.f90

module globalVars
implicit none
real, parameter :: pi = 3.1415927, hBar = 1.0, m = 1.0, w = 1.0
real, parameter :: h = 1e-4
end module
