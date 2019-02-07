!compile with: 
!gfortran -c -fdefault-real-8 tester.f90
!link and run with:
!gfortran tester.o polyFit.o -o tester && ./tester > output && cat output

program tester
use polyFit
implicit none

real :: SVDTestA(3,3), SVDTestU(3,3), SVDTestSig(3), SVDTestVT(3,3)
integer i, iost, numPoints
integer, parameter :: maxPoints = 1000000 !arbitrarily large
real xFile(maxPoints), yFile(maxPoints) !for reading the file
real, dimension(:), allocatable :: x, y !actual data arrays
real, dimension(:), allocatable :: yFit3svd, yFit7svd, yFitlss, fitResults
!Write format labels:
1 format(3f8.3)
2 format(3f8.3,a,1f9.3,a,3f7.3)
3 format(x,a,1f9.3,3f7.3)
4 format(x,a,1f15.3)
5 format(x,a,1f20.3)
6 format(x,a,1f13.3)
7 format(x,a,1f9.3,7f7.3)

!Problem 1:
write(*,*) "Problem 1: Testing single value decomposition of matrix A(3x3):"
SVDTestA = reshape((/1,2,3, 4,5,6, 7,8,9/), shape(SVDTestA))
do i = 1,3
	write(*,1) SVDTestA(i,:)
end do
call svd(3, SVDTestA, SVDTestU, SVDTestSig, SVDTestVT)
write(*,*) "SVD results:     U(3x3)| Sigma(3)| VT(3x3):"
do i = 1,3
	write(*,2) SVDTestU(i,:), "|", SVDTestSig(i), "|", SVDTestVT(i,:)
end do

!Problem 2/3:
write(*,*) ""
write(*,*) "Problem 2&3: Fitting noisy data with SVD, n=3. Then computing fit parameters, condition number, and X2."

! read data from file "data.dat", with limit maxPoints
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

! open file for writing
open(unit=2, file="svd3.csv", status="replace", action="write")
allocate(fitResults(3+1+3)) !3+1 for fit param; 3 for cond#/x2/dof
yFit3svd = yFit(x, y, numPoints, n=3, fitResults=fitResults, method=0)
write(2, *) "x, y, yFit3svd"
do i = 1, numPoints
	write (2, *) x(i), ",", y(i), ",", yFit3svd(i)
end do
close(2)
write(*,*) "Fit complete, results printed csv to be plotted later"
write(*,3) "SVD n=3 fit parameters:", fitResults(1:4)
write(*,4) "Condition number:", fitResults(5) 
write(*,5) "Chi-squared:", fitResults(6)
write(*,6) "Degrees of freedom:", fitResults(7)
write(*,*) "With the above results, the p-value of the chi-squared statistic is 0.991, making the fit very good."
deallocate(fitResults)

!Problem 4:
write(*,*) ""
write(*,*) "Problem 4: Fitting noisy data with SVD, n=7. Then computing fit parameters, condition number, and X2."
open(unit=2, file="svd7.csv", status="replace", action="write")
allocate(fitResults(7+1+3)) !7+1 for fit param; 3 for cond#/x2/dof
yFit7svd = yFit(x, y, numPoints, n=7, fitResults=fitResults, method=0)
write(2, *) "x, y, yFit7svd"
do i = 1, numPoints
	write (2, *) x(i), ",", y(i), ",", yFit7svd(i)
end do
close(2)
write(*,*) "Fit complete, results printed csv to be plotted later"
write(*,7) "SVD n=7 fit parameters:", fitResults(1:8)
write(*,4) "Condition number:", fitResults(9) 
write(*,5) "Chi-squared:", fitResults(10)
write(*,6) "Degrees of freedom:", fitResults(11)
write(*,*) "With the above results, the p-value of the chi-squared statistic is 0.250, making the fit poor."
deallocate(fitResults)

!Problem 5:
write(*,*) ""
write(*,*) "Problem 5: Fitting noisy data with LSS, n=3. Then computing fit parameters, condition number, and X2."
open(unit=2, file="lss.csv", status="replace", action="write")
allocate(fitResults(3+1+3)) !3+1 for fit param; 3 for cond#/x2/dof
yFitlss = yFit(x, y, numPoints, n=3, fitResults=fitResults, method=1)
write(2, *) "x, y, yFitlss"
do i = 1, numPoints
	write (2, *) x(i), ",", y(i), ",", yFitlss(i)
end do
close(2)
write(*,*) "Fit complete, results printed csv to be plotted later"
write(*,3) "LSS n=3 fit parameters:", fitResults(1:4)
write(*,4) "Condition number:", fitResults(5) 
write(*,5) "Chi-squared:", fitResults(6)
write(*,6) "Degrees of freedom:", fitResults(7)
write(*,*) "Compared to SVD fit, the LSS fit has a lower condition number, so it's  less sensitive to input errors."
write(*,*) "Fit parameters, chi-squared, and goodness of fit were all identical." 

deallocate (x, y, fitResults)
end program
