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

! read data from file "data.dat", at most n points
open(unit=1, file="data.dat", status="old", action="read")
do i = 1, maxPoints
	read (1, *, iostat=iost) xFile(i), yFile(i)
	if (iost < 0) exit
end do
close(1)

! check actual data extent
numPoints = i-1
if (numPoints == maxPoints) write (0,*) "Read data extent was truncated, recompile with larger maxPoints"

! allocate data arrays
allocate(x(numPoints), y(numPoints))
x = xFile(1:i)
y = yFile(1:i)

! write data back to stdout
!do i = 1,m
!	write (*,*) x(i), y(i)
!end do

deallocate (x, y)
end program
