! read.f90  -  read data from a file demo
! compile with: gfortran -O3 -fdefault-real-8 read.f90 

program read; implicit none

integer, parameter :: maxPoints = 1000000 !arbitrarily large
integer i, iost, numPoints
real xFile(n), yFile(n) !for reading the file
real, dimension(:), allocatable :: x, y !actual data arrays 

! read data from file "data.dat", at most n points
open(unit=1, file="data.dat", status="old", action="read")
do i = 1, maxPoints
	read (1, *, iostat=iost) xFile(i), yFile(i)
	if (status < 0) exit
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
