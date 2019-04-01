!compile with:
!gfortran -c -fdefault-real-8 tester.f90
!link modules with:
!gfortran -g tester.o globalVars.o integrate.o -o tester
!run with:
!./tester > output && cat output

program tester
use globalVars
use integration
implicit none

!Variable declarations
integer :: i
real :: E
character(len=40) :: fileOdd, fileEven

!Format labels, x = space, a = chars, nfw.d = n (f)loats (w)idth (d)ecimals
1 format (a,f3.1, a)


!***************Problem 1***************
write(*,*) "Problem 1: Integrating schrodinger's equation from x=0 to inf"
write(*,*) "Using V(x) = 1/2 * m *(wx)^2 with m = w = h = 1"
write(*,*) "Eigenvalues are (n+1/2)hw = (n+1/2) with n = 0, 1, 2..."
write(*,*) "prob1.pdf shows typical solutions of Psi where E is NOT an eigenvalue"
write(*,*) "Boundary values of Psi are "

do i=0,9
	E = 0.2*i
	!open file to write, using odd initial condition (declared as globalVar)
	write(fileOdd,1) "Q1Out/E=", E, "_Init=(0,1)_odd.dat"
	write(fileEven,1) "Q1Out/E=", E, "_Init=(0,1)_even.dat"
	open(unit=1, file=fileOdd, status="replace", action="write")
	open(unit=2, file=fileEven, status="replace", action="write")
	call integrate(0.2*i, initOdd, verbose = .true.)
	close(1); close(2)
	
	!now even initial conditions
	write(fileOdd,1) "Q1Out/E=", E, "_Init=(1,0)_odd.dat"
	write(fileEven,1) "Q1Out/E=", E, "_Init=(1,0)_even.dat"
	open(unit=1, file=fileOdd, status="replace", action="write")
	open(unit=2, file=fileEven, status="replace", action="write")
	call integrate(0.2*i, initEven, verbose = .true.) 
	close(1); close(2)
end do

!***************Problem 2***************
!***************Problem 3***************
!***************Problem 4***************
!***************Problem 5***************
	

end program

