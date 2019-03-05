!compile with:
!gfortran -c -fdefault-real-8 tester.f90
!link and run with:
!gfortran tester.o newtonMethod.o -o tester && ./tester > output && cat output

program tester
use newtonMethod
implicit none

!Variable declarations



!Write format labels, x = space, a = characters, nfw.d = n (f)loats (w)idth (d)ecimals
1 format(x, a, 1f20.16)

!Problem 1:
write (*,*) "Problem 1: Use Newton's method to solve x^3 -x +1/4 = 0 to double precision (eps = 1.0e-16)"
write (*,1) "Root x1: ", newton(f1, df1, x0=-1.0, eps=1.0e-16)
write (*,1) "Root x2: ", newton(f1, df1, x0=0.0, eps=1.0e-16)
write (*,1) "Root x3: ", newton(f1, df1, x0=1.0, eps=1.0e-16)

!Problem 3:

end program

