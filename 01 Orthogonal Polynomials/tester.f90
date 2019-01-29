!compile with: 
!gfortran -c -fdefault-real-8 tester.f90
!link and run with:
!gfortran tester.o gaussjordan.o polyapprox.o -o tester && ./tester > output && cat output

program tester
use gaussjordan
use polyapprox
implicit none

!All variables
integer i, n, m
real A(3,3), B(3)
real, dimension(:), allocatable :: x, y, xapprox, yapprox


!Problem 1:
n = 3
A(:,1) = [0.0, 2.0, 1.0]
A(:,2) = [1.0, -2.0, -3.0]
A(:,3) = [-1.0, 1.0, 2.0]
B = [-8.0, 0.0, 3.0]

write (*,*) "Problem 1: Solve the linear system A|B:"
do i=1,3
	write (*,*) A(:,i), "|", B(i)
end do
call gaussj(n, A, B)
write (*,*) "System solution:"
write (*,*) B

!Problem 2 (Sample size 10):
n = 10
m = 1000
allocate (x(n), y(n), xapprox(m), yapprox(m))
open(1, file = "freal_10.csv", status = "unknown")
x = uniform(n)
do i=1,n
	y(i) = f(x(i))
	write(1,*) x(i), y(i)
end do
close(1)

open(1, file = "fapprox_10.csv", status = "unknown")
xapprox = uniform(m)
yapprox = fapprox(x,y,xapprox,n,m,0)
do i=i,m
	write(1,*) xapprox(i), yapprox(i)
end do
close(1)
deallocate (x, y, xapprox, yapprox)

!Problem 2 (Sample size 100):
n = 100
m = 1000
allocate (x(n), y(n), xapprox(m), yapprox(m))
open(1, file = "freal_100.csv", status = "unknown")
x = uniform(n)
do i=1,n
	y(i) = f(x(i))
	write(1,*) x(i), y(i)
end do
close(1)

open(1, file = "fapprox_100.csv", status = "unknown")
xapprox = uniform(m)
yapprox = fapprox(x,y,xapprox,n,m,0)
do i=i,m
	write(1,*) xapprox(i), yapprox(i)
end do
close(1)
deallocate (x, y, xapprox, yapprox)

write (*,*) "Problem 2: Approximations to f(x) = 1/(1+10*x^2) printed to csv files"
write (*,*) "Plots will be shown after"

!Problem 3 (Sample size 100):
n = 100
m = 1000
allocate (x(n), y(n), xapprox(m), yapprox(m))
open(1, file = "greal_100.csv", status = "unknown")
x = uniform(n)
do i=1,n
	y(i) = g(x(i))
	write(1,*) x(i), y(i)
end do
close(1)

open(1, file = "gapprox_100.csv", status = "unknown")
xapprox = uniform(m)
yapprox = fapprox(x,y,xapprox,n,m,1)
do i=i,m
	write(1,*) xapprox(i), yapprox(i)
end do
close(1)
deallocate (x, y, xapprox, yapprox)

write (*,*) "Problem 3: Approximations to g(x) = f'(x) printed to csv files"
write (*,*) "Plots will be shown after"



end program
