!compile with: 
!gfortran -c -fdefault-real-8 tester.f90
!link and run with:
!gfortran tester.o gaussjordan.o -o tester && ./tester > output && cat output

program tester
use gaussjordan
implicit none

integer, parameter :: n = 3
real A(n,n), B(n)
integer i
A(1,:) = [1.0, 3.0, 4.0]
A(2,:) = [2.0, 7.0, 3.0]
A(3,:) = [2.0, 8.0, 6.0]
B = [3.0, -7.0, -4.0]

write (*,*) "Problem 1: Solve the linear system A|B:"
do i=1,3
	write (*,*) A(i,:), "|", B(i)
end do
call gaussj(n, A, B)
write (*,*) "System solution:"
write (*,*) B



end program
