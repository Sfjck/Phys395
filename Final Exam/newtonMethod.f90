!module containing code for Newton's method of approximation 
!compile with:
!gfortran -c -fdefault-real-8 newtonMethod.f90 -llapack
!Last edit: 3/5 12:24pm

module newtonMethod
implicit none

! Interface for calling the f&df functions
interface
    function f1Int(x)
        real, intent(in) :: x
        real :: f1Int
    end function f1Int
end interface

contains

! Newton's method of approximation, find roots of f(x)
function newton(f, df, x0, eps) result(x)
    procedure(f1Int) :: f, df
    real x0, eps, x, y; intent(in) x0, eps
    
	!starting points
    x = x0
    y = f(x)

	!repeat Newton until y is close enough to 0
    do while(eps < abs(y))
        x = x - y / df(x)
        y = f(x)
    end do
end function

! f1(x) = cos(x)-x/5
pure function f1(x)
    real f1, x; intent(in) x
    
    f1 = cos(x)-x/5.0
end function

! df1(x) = f1'(x) = -sin(x)-1/5
pure function df1(x)
    real df1, x; intent(in) x
    
    df1 = -1.0*sin(x)-0.2
end function
end module
