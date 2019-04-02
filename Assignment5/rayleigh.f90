! rayleigh.f90  -  Rayleigh iteration solver for eigenvalue problem
! compile with: gfortran -O3 -fdefault-real-8 rayleigh.f90 -llapack

module rayleigh
use globalVars
implicit none

contains

! hermite polynomials, calculated using recursive loop
 function hermite(i); intent(in) i
	real :: hermite(n), hermiteNext(n), hermitePrev(n)
	integer i, k

	if (i == 0) then
		hermite = 1
	else if (i == 1) then
		hermite = 2*x
	else
		hermitePrev = 1
		hermite = 2*x
		do k = 2,i
			hermiteNext = 2*x*hermite - 2*(k-1)*hermitePrev
			hermitePrev = hermite
			hermite = hermiteNext
		end do
	end if
end function

! initialize the collocation grid
subroutine initg()
	integer i
	forall (i=1:n) theta(i) = pi*(n-i+0.5)/n; x = ell/tan(theta)
end subroutine

! evaluate rational Chebyshev basis on a grid theta
subroutine evalb(n, pts, theta, Tn, Tnx, Tnxx)
	integer n, pts; real, dimension(pts), intent(in) :: theta
	real, dimension(pts), intent(out), optional :: Tn, Tnx, Tnxx
        
	! Chebyshev basis and its derivatives
	if (present(Tn))   Tn = cos(n*theta)
	if (present(Tnx))  Tnx = n * sin(n*theta) * sin(theta)**2/ell
	if (present(Tnxx)) Tnxx = -n * (n*cos(n*theta)*sin(theta) + 2.0*sin(n*theta)*cos(theta)) * sin(theta)**3/ell**2

end subroutine

! initialize linear spectral derivative operator
subroutine initl()
	integer i, pivot(n), status; real A(n,n), B(n,n)
	
	! evaluate basis and differential operator values on collocation grid
	do i = 1,n
		call evalb(i-1, n, theta, Tn=A(i,:), Tnxx=B(i,:))
	end do
	
	! find linear operator matrix
	status = 0
	select case (kind(A))
		case(4); call sgesv(n, n, A, n, pivot, B, n, status)
		case(8); call dgesv(n, n, A, n, pivot, B, n, status)
		case default; call abort
	end select
        
	! bail at first sign of trouble
	if (status /= 0) call abort
        
	L = transpose(B)
end subroutine

! Rayleigh's iteration solving eigenvalue problem:
function lsolve(psiV, eigenValue)
	real eigenValue, psiV(n), lsolve(n), A(n,n), B(n)
	integer i, pivot(n), status
	
	! linear improvement inflating eigenvector
	A = H; forall (i=1:n) A(i,i) = H(i,i) - eigenValue
	B = psiV
	
	! find linear operator matrix
        status = 0; select case (kind(A))
                case(4); call sgesv(n, 1, A, n, pivot, B, n, status)
                case(8); call dgesv(n, 1, A, n, pivot, B, n, status)
                case default; call abort
        end select
        
        ! bail at first sign of trouble
        if (status /= 0) call abort
        
        lsolve = B
end function
end module
