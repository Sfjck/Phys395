!compile with:
!gfortran -c -fdefault-real-8 tester.f90
!link and run with:
!gfortran tester.o newtonMethod.o bracketSearch.o && ./tester > output && cat output

program tester
use newtonMethod
use bracketSearch
implicit none

!Variable declarations
Real min1X, min2X

!Write format labels, x = space, a = characters, nfw.d = n (f)loats (w)idth (d)ecimals
1 format(x, a, 1f20.16)
2 format(x, a, 1f8.5, a, 1f8.5, a)

!Problem 1:
write (*,*) "Problem 1: Use Newton's method to solve x^3 -x +1/4 = 0 to double precision (eps = 1.0e-16)"
write (*,1) "Root x1 = ", newton(f1, df1, x0=-1.0, eps=1.0e-16)
write (*,1) "Root x2 = ", newton(f1, df1, x0=0.0, eps=1.0e-16)
write (*,1) "Root x3 = ", newton(f1, df1, x0=1.0, eps=1.0e-16)

!Problem 2&3:
write (*,*) ""
write (*,*) "Problem 2&3: Use bracketed search with phi to minimize (x^2-1)^2 -x"
min1X = bracketMin(f2, a0=-1.5, b0=-0.5, eps=1.0e-6)
min2X = bracketMin(f2, a0=0.5, b0=1.5, eps=1.0e-6)
write (*,2) "Minimum (x1, y1) = (", min1X, ", ", f2(min1X), ")"
write (*,2) "Minimum (x2, y2) = (", min2X, ", ", f2(min2X), ")"

!Problem 4&5:

end program

