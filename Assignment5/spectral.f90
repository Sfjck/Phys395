! spectral.f90  -  relaxation solver for BVP using spectral basis
! compile with: gfortran -O3 -fdefault-real-8 spectral.f90 -llapack

program spectral
use globalVars
use integration
implicit none

real, dimension(nS) :: x, psi, potential, temp1, temp2, temp3
real, dimension(nS,nS) :: L, hamiltonian, identity, matinv
real :: lamo, h
integer :: i

! initialize spectral operators
call initg(); call initl(); call inith()

identity = 0
do i = 1, nS
	identity(i,i) = 1
end do

! initial guess
call evalb(2, x, psi)

open(unit=1, file="RConv.dat")

do i = 1,20
	lamo = dot_product(psi, (matmul(hamiltonian, psi))) / dot_product(psi,psi)
	write(1,*) i, log10(abs(lamo-0.5))
	matinv = hamiltonian - lamo*identity
	psi = lsolve(matinv, psi)
	psi = psi/(sqrt(sum(psi**2)*h))
end do

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! spectral grid, basis functions and derivative operators
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! initialize the collocation grid
subroutine initg()
	integer i
	forall(i=1:nS) x(i) = atanh(cos(pi*(nS-i+0.5)/nS))
	forall(i=1:nS) potential(i) = V(x(i))
	h = (x(nS)-x(1))/nS
end subroutine

! evaluate rational Chebyshev basis on a grid theta
subroutine evalb(nS, x, Tn, Tnxx)
        integer nS, k
		real, dimension(nS) :: theta, x, Tn, Tnxx
        optional :: Tnxx
		
		theta = nS*acos(tanh(x))
		
		! Chebyshev basis and its derivatives
        Tn = cos(theta)
		if (nS > 1) then
			if (mod(nS,2) == 0) then
				Tn = Tn -1
			else if (mod(nS, 2) == 1) then
				Tn = Tn -x
			end if
		end if
        if (present(Tnxx)) Tnxx = -n * (sinh(x)*sin(theta) + n*cos(theta)) / (cosh(x)**2.0)
end subroutine evalb

! initialize linear spectral derivative operator
subroutine initl()
	integer i, pivot(nS), status; real A(nS,nS), B(nS,nS)
	
	! evaluate basis and differential operator values on collocation grid
	do i = 1,nS
		call evalb(i+1, x, A(i,:), B(i,:))
	end do
	
	! find linear operator matrix
	status = 0; select case (kind(A))
		case(4); call sgesv(nS, nS, A, nS, pivot, B, nS, status)
		case(8); call dgesv(nS, nS, A, nS, pivot, B, nS, status)
		case default; call abort
        end select
        
        ! bail at first sign of trouble
        if (status /= 0) call abort
        
        L = transpose(B)
end subroutine

subroutine inith()
	hamiltonian = -0.5 * L
	do i=i,nS
		hamiltonian(i,i) = Hamiltonian(i,i) + potential(i)
	end do
end subroutine

!solver for the classic Ax = B
function lsolve(A, B)
	real, dimension(nS) :: lsolve, psi, rhs, B
	real, dimension(ns,ns) :: A
	integer i, pivot(nS), status
	
	! find linear operator matrix
	status = 0; select case (kind(A))
		case(4); call sgesv(nS, 1, A, nS, pivot, B, nS, status)
		case(8); call dgesv(nS, 1, A, nS, pivot, B, nS, status)
		case default; call abort
	end select
        
	! bail at first sign of trouble
	if (status /= 0) call abort
        
	lsolve = B
end function
end program