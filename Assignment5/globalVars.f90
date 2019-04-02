!module containing global variables
!compile with:
!gfortran -c -fdefault-real-8 globalVars.f90

module globalVars
implicit none
real, parameter :: pi = 3.14159259, eu = 2.71828, hBar = 1.0, m = 1.0, w = 1.0, lambda = 1.0, ell = 1.0
real, parameter :: dx = 1.0e-4, modx = 1e-2, eps = 1.0e-12
!initial conditions [Psi, dPsi] to produce purely odd/even Psi
real, parameter :: initOdd(2) = [0.0, 1.0], initEven(2) = [1.0, 0.0]
real :: f(2), norm0
integer, parameter :: n = 150 !spectral order
logical :: eigenFound, harmonic

real, dimension(n) :: x, theta
real, dimension(n,n) :: L, H

contains

!V of harmonic / anharmonic oscillator
pure function V(xx);	intent(in) xx
	real :: V, xx
	if (harmonic .eqv. .true.) then
		V = 0.5*m*(w*xx)**2.0
	else
		V = 0.25*lambda*(xx**4.0)
	end if
end function

end module
