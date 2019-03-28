!compile with:
!gfortran -c -fdefault-real-8 tester.f90
!link modules with:
!gfortran -g tester.o globalVars.o integrate.o -o tester
!run with:
!./tester > output && cat output

program tester
use globalVars
use integrate
implicit none

!Variable declarations
real :: init(2), init2(2)
integer :: i

!Format labels, x = space, a = chars, nfw.d = n (f)loats (w)idth (d)ecimals
1 format (a,f3.1, a, f2.0, a, f2.0, a)

!***************Problem 1***************
write(*,*) "Problem 1: Calculating typical solutions of Psi_+ and Psi_-"
write(*,*) "Refer to prob1.pdf for plots"

init = [0.0, 1.0]
init2 = [1.0, 0.0]
do i=0,9
	call integrate(0.2*i, init)
	call integrate(0.2*i, init2)
end do

!***************Problem 2***************
!***************Problem 3***************
!***************Problem 4***************
!***************Problem 5***************
	

end program

