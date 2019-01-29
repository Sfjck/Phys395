!module containing code for polynomial approximation
!compile with:
!gfortran -c -fdefault-real-8 polyapprox.f90

module polyapprox
use gaussjordan
implicit none

contains

! uniformly spaced grid of n points covering interval [-1,1]
pure function uniform(n)
	integer i, n; real uniform(n); intent(in) n
	
	! uniform grid
	forall (i=1:n) uniform(i) = (2*i-n-1.0)/(n-1.0)
end function


! function to be approximated
elemental function f(x)
	real f, x; intent(in) x
	
	f = 1.0/(1.0 + 10.0*x*x)
end function

! approximate f(x) using chebyshev basis
function fapprox(x, y, xapprox, n, m)
	real, intent(in) :: x(n), y(n), xapprox(m)
	integer, intent(in) :: n, m
	integer i, j
	real fapprox(m), A(n,n), B(1,n), AT(n,n)

	do i=1,n
		do j=1,n
			A(j,i)=basis(j-1, x(i))
		end do
		B(1,i) = y(i)
	end do
	
	
!	AT = transpose(A)
	call gaussj(n,A,B)
	fapprox(1:m) = basisSums(B(1,:), xapprox, n, m)
end function

elemental function basis(n, x)
	integer n; real x
	intent(in) n, x
	real basis

	! Chebyshev polynomials
	basis = cos(n*acos(x))
end function

function basisSums(b, xapprox, n, m)
	real, intent(in) :: b(n), xapprox(m)
	integer, intent(in) :: n, m
	integer :: i,j
	real::basisSums(m), sums

	do i=1,m
	sums = 0
		do j=1,n
			sums = sums + b(j)*cos((j-1)*acos(xapprox(i)))
			basisSums(i)= sums
		end do
	end do
end function

end module
