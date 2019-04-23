!compile with:
!gfortran -c -fdefault-real-8 Q1_2.f90
!link and run with:
!gfortran Q1_2.o newtonMethod.o bracketSearch.o && ./Q1_2 > output1_2 && cat output1_2

program tester
use newtonMethod
use bracketSearch
implicit none

!Variable declarations
real minX1, minX2, maxX1

!Format labels, x = space, a = characters, nfw.d = n (f)loats (w)idth (d)ecimals
1 format(x, a, 1f16.12)
2 format(/,x, a)
3 format(x, a, 1f16.12, a, 1f16.12, a)
4 format(a5, i1, a, 1f11.8)
5 format(a9, 1f11.8)
6 format(a9, 1f11.5)

!Problem 1:
write (*,*) "Problem 1: Use Newton's method to solve cos(x)-x/5 = 0 to 1e-12 precision"
write (*,1) "Root x1 = ", newton(f1, df1, x0=-4.0, eps=1.0e-12)
write (*,1) "Root x2 = ", newton(f1, df1, x0=-2.0, eps=1.0e-12)
write (*,1) "Root x3 = ", newton(f1, df1, x0=1.0, eps=1.0e-12)

!Problem 2:
write (*,2) "Problem 2: Use bracketed search with phi to minimize/maximize x^4+3x^3-4x^2-3x+4"
minX1 = bracketMin(f2, a0=-3.0, b0=-2.5, eps=1.0e-12)
maxX1 = bracketMin(f2neg, a0=-0.5, b0=0.5, eps=1.0e-12)
minX2 = bracketMin(f2, a0=0.5, b0=1.0, eps=1.0e-12)
write (*,3) "Minimum (x1, y1) = (", minX1, ", ", f2(minX1), ")"
write (*,3) "Maximum (x2, y2) = (", maxX1, ", ", f2(maxX1), ")"
write (*,3) "Minimum (x3, y3) = (", minX2, ", ", f2(minX2), ")"

end program

