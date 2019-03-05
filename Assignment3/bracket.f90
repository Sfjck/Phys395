!Copied from Andrei

module bracket
implicit none

! Interface for calling the f&df functions
interface
    function real2real(x)
        real, intent(in) :: x
        real :: real2real
    end function real2real
end interface

contains

	function bracket(f, xMin, xMax, eps)
		procedure(real2real) :: f
		real xMin, xMax, eps, bracket, a, b, c, d; intent(in) :: xMin, xMax, eps
		real, parameter :: phi = 1.61803398874989484820458683436563811772030917980576286
	
		!starting points
		a = xMin
		b = xMax
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
		
		bracket = 0.5 * (b+a)
	end function

elemental function f2(x)
	real f2, x; intent(in) x
	
	f2 = (x**2.0 - 1.0)**2.0 + x
end function
end module