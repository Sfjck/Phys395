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
integer :: i, j
integer, parameter :: n=5
real :: E, EMinOdd, EMinEven, ERange
real :: y(2), EigenNormEven(n,2), EigenNormOdd(n,2)
!NOTE: delete y(2)?
character(len=80) :: fileOdd, fileEven, fileNorm

!Format labels for write
1 format (a,f3.1, a)
2 format (/, a)
3 format (x, a6, 5f6.3)

!***************Problem 1***************
write(*,*) "Problem 1: Integrating schrodinger's equation from x=0 towards inf"
write(*,*) "Using V(x) = 1/2 * m *(wx)^2 with m = w = h = 1"
write(*,*) "Eigenvalues are (n+1/2)hw = (n+1/2) with n = 0, 1, 2..."
write(*,*) "prob1.pdf shows typical solutions of Psi where E is NOT an eigenvalue"

do i=0,9
	E = 0.2*i
	!open file to write, using odd initial condition (declared as globalVar)
	write(fileOdd,1) "Q1Out/E=", E, "_Init=(0,1)_odd.dat"
	write(fileEven,1) "Q1Out/E=", E, "_Init=(0,1)_even.dat"
	open(unit=1, file=fileOdd, status="replace", action="write")
	open(unit=2, file=fileEven, status="replace", action="write")
	call integrate(0.2*i, initOdd, problem = 1)
	close(1); close(2)
	
	!now even initial conditions
	write(fileOdd,1) "Q1Out/E=", E, "_Init=(1,0)_odd.dat"
	write(fileEven,1) "Q1Out/E=", E, "_Init=(1,0)_even.dat"
	open(unit=1, file=fileOdd, status="replace", action="write")
	open(unit=2, file=fileEven, status="replace", action="write")
	call integrate(0.2*i, initEven, problem = 1) 
	close(1); close(2)
end do

!***************Problem 2***************
write(*,2) "Problem 2: Integrating for Psi and |Psi|^2 with Harmonic Potential"
!determine the exact eigenvalues (n+1/2) and nomalization factor
do i = 1,n
	EigenNormOdd(i,:) = bisect(2.0*(i-1), 2.0*i, initOdd)
	EigenNormEven(i,:) = bisect(2.0*(i-1), 2.0*i, initEven)
end do
write(*,*) "Calculated eigenvalues: (expecting n+1/2, n=0,1,2...)"
write(*,3) "Odd: ", EigenNormOdd(:,1)
write(*,3) "Even: ", EigenNormEven(:,1)

!calculate Psi and PDF of Psi (ie |Psi|^2/Norm0)
do i = 1,n
	!odd initial conditions
	write(fileNorm, 1) "Q2Out/PsiNorm, E=", EigenNormOdd(i,1), ".dat"
	open(unit=1, file=fileNorm, status="replace", action="write")
	norm0 = EigenNormOdd(i,2)
	call integrate(EigenNormOdd(i,1), initOdd, problem = 2)
	close(1)

	!even initial condition
	write(fileNorm, 1) "Q2Out/PsiNorm, E=", EigenNormEven(i,1), ".dat"
	open(unit=1, file=fileNorm, status="replace", action="write")
	norm0 = EigenNormEven(i,2)
	call integrate(EigenNormEven(i,1), initEven, problem = 2)
	close(1)
end do

!***************Problem 3***************
write(*,2) "Problem 3: Integrating for Psi and |Psi|^2 with Anharmonic Potential"
!determine the exact eigenvalues and nomalization factor
i=1; j=1
EMinOdd=0.0; EMinEven=0.0; ERange=1
do while (i+j .le. 2*n) !NOTE: can use only i? from 1 to 5?
	EigenNormOdd(i, :) = bisect(EMinOdd, EMinOdd + ERange, initOdd)
	if (eigenFound .eqv. .false.) then
		EMinOdd = EMinOdd + ERange
	else
		EMinOdd = EigenNormOdd(i, 1) + ERange !NOTE: might need to be much less? (epsilon)
		i = i+1
	end if
	
	EigenNormEven(j, :) = bisect(EMinEven, EMinEven + ERange, initEven)
	if (eigenFound .eqv. .false.) then
		EMinEven = EMinEven + ERange
	else
		EMinEven = EigenNormEven(j, 1) + ERange
	end if
end do

!calculate Psi and PDF of Psi (ie |Psi|^2/Norm0)

	
end do
!***************Problem 4***************
!***************Problem 5***************
	

end program

