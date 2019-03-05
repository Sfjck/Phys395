!compile with:
!gfortran -c -fdefault-real-8 tester.f90
!link and run with:
!gfortran tester.o newtonMethod.o -o tester && ./tester > output && cat output

program tester
use newtonMethod
implicit none

!Variable declarations



!Write format labels
1 format(

!Problem 1:
write (*,*) "Problem 1: Use Newton's method to solve x^3 -x +1/4 = 0 to double precision (eps = 1.0e-16)"
write (*,*) newton(f1, df1, x0=-1.0, tol=1.0e-16)

end program

