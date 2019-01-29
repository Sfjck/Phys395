! ass1- prob1
! gaussj.f90  - tester and subroutine for Gauss-Jordan elimination
! compile with: gfortran -fdefault-real-8 gaussj.f90 -o P1 && ./P1 > OutP1 && cat OutP1

program testP1
implicit none

integer, parameter :: n = 3
real A(n,n), B(n)
integer i

A(1,:) = [1.0, 3.0, 4.0]
A(2,:) = [2.0, 7.0, 3.0]
A(3,:) = [2.0, 8.0, 6.0]
B = [3.0, -7.0, -4.0]

write (*,*) "Test linear equation A|B:"
do i=1,3
	write (*,*) A(i,:), "|", B(i)
end do
call gaussj(n, A, B)
write (*,*) "Equation solution:"
write (*,*) B

contains
	
! solve A.x = B using Gauss-Jordan elimination
! A gets destroyed, answer is returned in B
subroutine gaussj(n, A, B)
    integer n; real A(n,n), B(n)
    integer i, j
    real tempA(n), tempB, maxA !for pivots ^ row swapping
    integer maxALoc
    
    ! get the matrix into upper-diagonal form first
    do i = 1,n
        ! find pivot element
        maxALoc = i
        maxA = A(i,i)
        do j = i+1,n
            if (abs(A(j,i)) > abs(maxA)) then
                maxA = A(j,i)
                maxALoc = j
            end if
        end do
        ! swap pivot to i-th row
        if (maxALoc /= i) then
            tempA = A(i,:)
            A(i,:) = A(maxALoc,:)
            A(maxALoc,:) = tempA
            tempB = B(i)
            B(i) = B(maxALoc)
            B(maxALoc) = tempB
        end if
        ! eliminate lower diagonal elements
        do j = i+1,n
            B(j) = B(j) - A(j,i)/A(i,i) * B(i)
            A(j,:) = A(j,:) - A(j,i)/A(i,i) * A(i,:)
        end do
    end do

    ! solve the equations by back-substitution
    do j = n,1,-1
        tempB = B(j)
        do i = j+1, n
            tempB = tempB - A(j,i)*B(i)
        end do
        B(j) = tempB / A(j,j)
    end do
end subroutine

subroutine writeMatrix(M, jMax, iMax)
	real M(jMax,iMax); integer jMax,iMax
	integer j,i

	do j=1,jMax
		write (*,*) M(:,i)
	end do
end subroutine

end program
