!module containing code for Newton's method of approximation 
!compile with:
!gfortran -c -fdefault-real-8 newtonMethod.f90 -llapack

module newtonMethod
implicit none

! Interface for calling the f&df functions
interface
    function real2real(x)
        real, intent(in) :: x
        real :: real2real
    end function real2real
end interface

contains

! Newton's method of approximation, find roots of f(x)
function newton(f, df, x0, eps) result(x)
    procedure(real2real) :: f, df
    real x0, eps, x, y; intent(in) x0, eps
    
	!starting points
    x = x0
    y = f(x)

	!until y is close enough to 0 (within eps)
    do while(eps < abs(y))
        x = x - y / df(x)
        y = f(x)
    end do
end function


! f1(x) = x^3 -x +1/4
pure function f1(x)
    real f1, x; intent(in) x
    
    f1 = x**3.0 -x +0.25
end function

! df1(x) = f1'(x) = 3x^2 -1
pure function df1(x)
    real df1, x; intent(in) x
    
    df1 = 3.0 * x**2.0 -1.0
end function
end module
