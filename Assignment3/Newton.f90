!module containing code for Newton's method of approximation 
!compile with:
!gfortran -c -fdefault-real-8 newtonMethod.f90 -llapack

module newtonMethod
implicit none

! Interface for calling the f&g functions
interface
    function real2real(x)
        real, intent(in) :: x
        real :: real2real
    end function real2real
end interface

contains

! Newton's method of approximation, find roots of f(x)
function newton (f, g, x0, tol)
    procedure(real2real) :: f, g
    real, intent(in) :: x0, tol
    real :: x, y
    
    x = x0
    y = f(x)

    do while(tol < abs(y))
        x = x - y / g(x)
        y = f(x)
    end do
end function


! f1(x) = x^3 -x +1/4
pure function f1(x)
    real, intent(in) :: x
    real :: f1
    
    f1 = x**3.0 -x +0.25
end function

! g1(x) = f1'(x) = 3x^2 -1
pure function g1(x)
    real, intent(in) :: x
    real :: f1
    
    g1 = 3.0 * x**2.0 -1.0
end function
end module
