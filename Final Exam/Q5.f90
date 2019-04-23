!compile with:
!gfortran -c -fdefault-real-8 tester.f90
!link modules with:
!gfortran -g tester.o globalVars.o integrate.o -o tester
!run with:
!./tester > output && cat output

program Q5
use globalVars5
use rayleigh
implicit none

!Variable declarations
integer :: i, j, k
integer, parameter :: nn=5
real :: E, EMinOdd, EMinEven, ERange
real, dimension(nn,2) :: EigenNormEven, EigenNormOdd
real, dimension(n) :: gaussian, psiV, psiVN !V is for vector, N is norm
character(len=100) :: fileOdd, fileEven, filePsi

!Format labels for write
1 format (a,f3.1, a)
2 format (/, a)
3 format (x, a6, 5f6.3)
4 format (a,f4.1, a)
5 format (f20.16)

!***************Problem 5***************
write(*,2) "Problem 5: Using Spectral method to calculate first 15 eigenvalues of quantum oscillator"

! initialize spectral operators
call initg()
call initl()
gaussian = eu**(-x**2/2.0)

! hamiltonian for harmonic oscillator
harmonic = .false.
H = -L/2.0; forall (i=1:n) H(i,i) = -L(i,i)/2.0 + V(x(i))

!determine eigenvalues
do i = 0,15!2*nn-1
	! initial guess of the wavefunction
	psiV = hermite(i) * gaussian

	! get Psi by relaxing with Rayleigh's iteration
	do k = 1,50
		E = dot_product(psiV,matmul(H,psiV))/dot_product(psiV,psiV)
		psiV = lsolve(psiV,E) 
		psiV = psiV/sqrt(dot_product(psiV,psiV))
	end do
	write (*,5) E
end do
end program

