!module containing code for fitting with SVD
!compile with:
!gfortran -c -fdefault-real-8 polysvd.f90 -llapack

module polysvd
implicit none
contains

! fit function for single value decomp
function fitsvd(x, y, numPoints, n) results(yFit)
	integer :: numPoints, n, a, i, j, k
	real, dimension(:) :: x, y
	real, dimension(numPoints) :: yFit
	real, dimension(numPoints, n+1) :: b_a
	real, dimension(n+1, n+1) :: B, U, VT
	real :: chi2Act, chi2Exp, condNum

	! minimize |y - b_a x|^2
	! solve Bc = p
	b_a = basis(numPoints, x, n)
	
	forall (j=1:n+1), k=1:n+1) B(j,k) = sum(b_a(:, j) * b_a(:, k)
	p = matmul(transpose(b_a), y)

	call svd(n+1, B, U, sig, VT)

	coeffs = svdSolver(p, U, sig, VT, eps = 1.0e-6)

	yFit = matmul(bxa, coeffs)
	chi2Act = sum((y-yFit)**2)
	chi2Exp = size(x) - (n+1)
	condNum = s(1) / s(n+1)
end function

! solve (U*s*Vh)x = b
function svdSolver(b, U, sig, VT, eps)! result(x)
	real :: b(:), s(:), U(:, :), x(size(VT, 2)), eps

	x = matmul(transpose(U), b)
	where (sig > eps * sig(1)); x = x / sig
		elsewhere; x = 0.0
	end where
	x = matmul(transpose(VT), x)
end function

! use dgesvd to calculate U, sig, VT from MxN A. VT = V transposed, sig=diagonal of sigma
subroutine svd(n, A, U, sig, VT)
	integer :: n, LWork, info
	real :: A(n,n), ACopy(n,n), U(n, n), sig(n), VT(n, n)
	real, dimension(:), allocatable :: work

	LWork = 6*n !must be at least 5*n, manually tested for optimal
	allocate(work(LWork))
	ACopy = A !dgesvd destroys ACopy so A  is untouched
	!DGESVD(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
	call dgesvd('A', 'A', n, n, ACopy, n, sig, U, n, VT, n, work, LWork, info)
	if (info /= 0) call abort
	if (Work(1) /= LWork) write(*,*) "LWork was not optimal, optimal Lwork=", LWork
end subroutine

! basis is b_a(x)=cos(2pi*a*x)
pure function basis(numPoints, x, n)
	integer numPoints, n, a, i
	real x(numPoints), basis(numPoints, n+1)
	intent(in) numPoints, x, n
	real, parameter :: pi = 3.14159265358979323846264338327950288419716939937510

	forall (a=0:n, i=1:numPoints) basis(i, a+1) = cos(2*pi*a*x(i))
end function

end module
