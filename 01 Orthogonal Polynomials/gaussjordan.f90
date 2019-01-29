!module containing code for Gauss-Jordan elimination
!compile with:
!gfortran -c -fdefault-real-8 gaussjordan.f90

module gaussjordan
implicit none
contains

! solve A.x = B
! A gets destroyed, answer is returned in B
subroutine gaussj(n, A, B)
    integer n; real A(n,n), B(n)
    integer i, j
    real tempA(n), tempB, maxA; integer maxALoc !for pivots & row swapping
    
    ! get matrix into upper-triangular form
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
        ! eliminate lower triangle elements
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
end module
