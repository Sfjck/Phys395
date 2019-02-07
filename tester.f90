!compile with: 
!gfortran -c -fdefault-real-8 tester.f90
!link and run with:
!gfortran tester.o polysvd.o -o tester && ./tester > output && cat output

program tester
use polysvd
implicit none

integer, parameter :: maxPoints = 1000000 !arbitrarily large
real xFile(maxPoints), yFile(maxPoints) !for reading the file
integer i, iost, numPoints
real, dimension(:), allocatable :: x, y !actual data arrays
real, dimension(:), allocatable :: yFit3svd, yFit7svd, yFit3lss, yFit7lss

! read data from file "data.dat", at most n points
	write(*,*) "123"

open(unit=1, file="data.dat", status="old", action="read")
do i = 1, maxPoints
	read (1, *, iostat=iost) xFile(i), yFile(i)
	if (iost < 0) exit
end do
close(1)

! check actual data extent
	write(*,*) "a"
numPoints = i-1
if (numPoints == maxPoints) write (0,*) "Read data extent was truncated, recompile with larger maxPoints"

! allocate data arrays
allocate(x(numPoints), y(numPoints))
x = xFile(1:i)
y = yFile(1:i)

!Problem 2
open(unit=2, file="svd3.csv", status="replace", action="write")
	write(*,*) "b"
yFit3svd = fitsvd(x, y, numPoints, n=3)
write(2, *) "x, y, yFit3svd(x)"
do i = 1, numPoints
	write (2, *) x(i), ",", y(i), ",", yFit3svd
end do
close(2)

deallocate (x, y)
end program
