!compile with:
!gfortran -c -fdefault-real-8 tester.f90
!link and run with:
!gfortran tester.o gaussLeg.o && ./tester > output && cat output

program tester
use gaussLeg
implicit none

!Variable declarations


!Format labels, x = space, a = chars, nfw.d = n (f)loats (w)idth (d)ecimals


!Problem 1:


!Problem 2:


!Problem 3:



end program

