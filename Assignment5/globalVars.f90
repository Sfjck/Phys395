!module containing global variables
!compile with:
!gfortran -c -fdefault-real-8 globalVars.f90

module globalVars
implicit none
real, parameter :: pi = 3.1415927, hBar = 1.0, m = 1.0, w = 1.0
real, parameter :: dx = 1.0e-4, modx = 1e-2, eps = 1.0e-12
real, parameter :: initOdd(2) = [0.0, 1.0] !initial conditions [Psi, dPsi] to produce purely odd/even Psi
real, parameter :: initEven(2) = [1.0, 0.0]
real :: f(2), norm0
end module
