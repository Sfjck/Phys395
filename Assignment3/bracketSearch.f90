!module containing code for bracketed minimum search by golden sections
!compile with:
!gfortran -c -fdefault-real-8 bracketSearch.f90 -llapack

module bracketSearch
implicit none

! Interface for calling the f&df functions
interface
    function real2real(x)
        real, intent(in) :: x
        real :: real2real
    end function real2real
end interface

contains

function bracketMin(f, a0, b0, eps)
	procedure(real2real) :: f
	real a0, b0, eps, bracketMin, a, b, c, d; intent(in) :: a0, b0, eps
	real, parameter :: phi = 0.5 + sqrt(5.0)/2.0
!	write(*,*) "phi=", phi
	
	!starting points
	a = a0
	b = b0
	c = b - (b-a)/phi
	d = a + (b-a)/phi
	
	do while (eps < abs(c-d))
		if (f(c) < f(d)) then
			b = d
		else
			a = c
		end if
		
		c = b - (b-a)/phi
		d = a + (b-a)/phi
	end do
	
	bracketMin = 0.5 * (b+a)
end function

pure function f2(x)
	real f2, x; intent(in) x
	
	f2 = (x**2.0 - 1.0)**2.0 + x
end function
end module
