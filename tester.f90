!compile with: 
!gfortran -c -fdefault-real-8 tester.f90
!link and run with:
!gfortran tester.o gaussjordan.o polysvd.o -o tester && ./tester > output && cat output

program tester
use polysvd
implicit none

!All variables
integer i, n, m
real A(3,3), B(3)
real, dimension(:), allocatable :: x, y, xapprox, yapprox, error
real :: basix(100)
integer a

n = 3
forall (a=0:n, i=1:size(bxa, 1)) bxa(i, a + 1) = basis(a, x(i))






end program
