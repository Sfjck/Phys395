!compile with:
!gfortran -c -fdefault-real-8 ICSweep.f90 -lcfitsio
!link with
!gfortran -g ICSweep.o gaussLeg.o globalVars.o -o ICSweep

program ICSweep
use gaussLeg
use globalVars
implicit none

!Variable declarations
!Time step bigger here since we're running a LOT more simulations
character(len=8000) :: arg
real, parameter :: dt = 0.05
real :: y(4), t, tMax, theta1, theta2
integer i, j, N, nx, ny, status
real :: xrange(2) = [-pi, pi], yrange(2) = [-pi, pi]
real(4), allocatable :: image(:,:,:)

!Problem 3: Scan initial conditions of double pendulum for flip-time
!Get input N for sweep division
N = get_command_argument(1, arg)
write (*,*) trim(arg)
read (arg, *, iostat=status) N
if (status == 0 .and. (N .ge. 1) .and. N .le. 1000) then
	write (*,*) "Parses as integer: ", N
	nx = 2*N
	ny = 2*N
	allocate(image(1,nx,ny))
else 
	write (*,*) "Input argument invalid, enter integer between 1 and 1000"
	exit
end if
	
!Do sweep, no zoom
tMax = 1000.0*sqrt(l/gr)
do i=-N+1, N
	theta1 = pi * i/N
	do j=-N+1, N
		t = 0.0
		theta2 = pi * j/N
		write(*,*) theta1, theta2
		y = [theta1, theta2, 0.0, 0.0]
		do while ((abs(y(2)) .le. pi) .and. (t .le. tMax))
			call gl8(y, dt)
			t = t+dt
		end do
		image(1, i+N, j+N) = t
	end do
end do
call write2fits('flipTime1.fits', image, xrange, yrange, ['N'])
write(*,*)

contains

!write array data into FIT file as sequence of image extensions, from class
subroutine write2fits(file, array, xx, yy, vars, coords)
	character(len=*) file, vars(:), coords
	real(4) array(:,:,:); real(8) xx(2), yy(2)
	optional xx, yy, vars, coords
        
	integer i, j, status, unit
	integer :: hdus, naxis = 2, n(2), npix
	integer :: bitpix = -32, group = 1, blocksize = -1
        
	! data dimensions
	hdus = size(array,1)
	n(1) = size(array,2)
	n(2) = size(array,3)
	npix = n(1)*n(2)
        
	! delete file if it already exists
	open(unit=1234, iostat=status, file=file, status='old')
	if (status == 0) close(1234, status='delete'); status = 0
        
	! initialize FITS file
	call ftgiou(unit, status)
	call ftinit(unit, file, blocksize, status)
        
	! write image extensions
	do i = 1,hdus
		call ftiimg(unit, bitpix, naxis, n, status)
		call ftppre(unit, group, 1, npix, array(i,:,:), status)
                
		if (present(vars)) then
			if (present(coords)) then
				call ftpkys(unit, 'EXTNAME', vars(i)//coords, 'variable stored in extension', status)
			else
				call ftpkys(unit, 'EXTNAME', vars(i), 'variable stored in extension', status)
			end if
		end if
		if (present(xx)) then
 			call ftpkyj(unit, 'CRPIX1', 1, 'x-axis origin pixel', status)
			call ftpkyd(unit, 'CRVAL1', xx(1), 14, 'x-axis origin coordinate', status)
			call ftpkyd(unit, 'CDELT1', (xx(2)-xx(1))/n(1), 14, 'x-axis increment', status)
		end if
		if (present(yy)) then
			call ftpkyj(unit, 'CRPIX2', 1, 'y-axis origin pixel', status)
			call ftpkyd(unit, 'CRVAL2', yy(1), 14, 'y-axis origin coordinate', status)
			call ftpkyd(unit, 'CDELT2', (yy(2)-yy(1))/n(2), 14, 'y-axis increment', status)
		end if
	end do
        
	! clean up
	call ftclos(unit, status)
	call ftfiou(unit, status)
end subroutine
end
