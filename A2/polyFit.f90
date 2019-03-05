!module containing code for polynomial fitting with SVD and LSS 
!compile with:
!gfortran -c -fdefault-real-8 polyFit.f90 -llapack

module polyFit
implicit none
contains

! fit function, uses SVD(method=0) or LSS(method=1)
function yFit(x, y, numPoints, n, fitResults, method)
	integer :: numPoints, n, method, a, i, j, k
	real, dimension(:) :: x, y, fitResults
	real, dimension(numPoints+1, n+1) :: b_a
	real, dimension(n+1, n+1) :: b_ab, U, VT
	real, dimension(n+1) :: fitParams, b, sig
	real, dimension(numPoints+1) :: yFit
	real :: chi2, DoF, condNum

	! set up basis matrices
	b_a = basis(numPoints, x, n) !want |y - b_a(x)|2 --> 0
	b_ab = basisSOP(numPoints, b_a, n)
	b = matmul(transpose(b_a), y) !solve b_ab*c = b

	!solving using chosen method
	if (method == 0) then !SVD
		call svd(n+1, b_ab, U, sig, VT)
		fitParams = svdFit(b, U, sig, VT) !=c_a
	else !LSS
		call lss(n+1, numPoints+1, b_a, y, -1.0, fitParams, sig)
	end if

	!calculate fit results
	yFit = matmul(b_a, fitParams)
	condNum = sig(1) / sig(n+1)
	DoF = numPoints - (n+1)
	chi2 = sum((y-yFit)**2) !std=1

	!populate fit results to return
	fitResults(1:n+1) = fitParams
	fitResults(n+2) = condNum
	fitResults(n+3) = chi2
	fitResults(n+4) = DoF 
end function

!pretty much same as subroutine svd below but with lss
subroutine lss(n, numPoints, A, B, rcond, fitParams, sig)
	integer :: n, numPoints, LWork, rank, info
	real :: A(numPoints,n), B(:), rcond, fitParams(n), sig(n)
	real :: ACopy(numPoints, n), BCopy(numPoints)
	real, dimension(:), allocatable :: work

	LWork=8*n !must be at least 5*n
	allocate(work(LWork))
	ACopy=A
	BCopy=B
	!DGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, LWORK, INFO )
	call dgelss(numPoints, n, 1, ACopy, numPoints, BCopy, numPoints, sig, rcond, rank, work, LWork, info)
	fitParams = BCopy(1:n)
	if (info /= 0) call abort
!	if (Work(1) /= LWork) write(*,*) "Optimal Lwork was", Work(1)
end subroutine

! solve (U * sig * VT)x = b for x, using SVD
function svdFit(b, U, sig, VT)
	real :: b(:), U(:, :), sig(:), VT(:,:), svdFit(size(VT, 2))

	svdFit = matmul(transpose(U), b) !UT*U = I; --> svdFit = UT*b = (sig*VT)x
	where (sig > 10e-6) ! prevent div by 0
		svdFit = svdFit / sig	!--> svdFit = UT*b/sig = VT*x
	elsewhere
		svdFit = 0.0
	end where
	svdFit = matmul(transpose(VT), svdFit) !VT*V = I; --> svdFit = V*UT*b/sig = x
end function

! use dgesvd to calculate U, sig, VT from MxN A. VT = V transposed, sig=diagonal of sigma
subroutine svd(n, A, U, sig, VT)
	integer :: n, LWork, info
	real :: A(n,n), ACopy(n,n), U(n, n), sig(n), VT(n, n)
	real, dimension(:), allocatable :: work

	LWork = 8*n !must be at least 5*n
	allocate(work(LWork))
	ACopy = A !dgesvd destroys ACopy so A  is untouched
	!DGESVD(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
	call dgesvd('A', 'A', n, n, ACopy, n, sig, U, n, VT, n, work, LWork, info)
	if (info /= 0) call abort
!	if (Work(1) /= LWork) write(*,*) "Optimal Lwork was", Work(1)
end subroutine

! basis is b_a(x)=cos(2pi*a*x)
pure function basis(numPoints, x, n)
	integer numPoints, n, a, i
	real x(numPoints+1), basis(numPoints+1, n+1)
	intent(in) numPoints, x, n
	real, parameter :: pi = 3.14159265358979323846264338327950288419716939937510

	forall (a=0:n, i=1:numPoints+1) basis(i, a+1) = cos(2*pi*a*x(i))
end function

! basisSOP is sum of product of b_a, b_b from eq3
pure function basisSOP(numPoints, b_a, n)
	integer numPoints, n, a, b
	real b_a(numPoints+1, n+1), basisSOP(n+1, n+1)
	intent(in) numPoints, b_a, n

	forall (a=1:n+1, b=1:n+1) basisSOP(a,b) = sum(b_a(:, a)*b_a(:, b))
end function

end module
