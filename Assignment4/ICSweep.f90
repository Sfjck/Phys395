!compile with:
!gfortran -c -fdefault-real-8 ICSweep.f90 -lcfitsio
!link with
!gfortran -g ICSweep.o gaussLeg.o globalVars.o -o ICSweep

program ICSweep
use gaussLeg
use globalVars
implicit none

!Variable declarations
!Bigger dt for many simulations, [ZoomX, ZoomY] are focal points of zooms
real, parameter :: dt = 0.05, ZoomX = -0.99, ZoomY = 1.63
real :: y(4), t, tMax, theta1, theta2
integer i, j, N, status, numArgs, zoomPower, zoom
character(len=8000) :: arg, filename
real :: xrange(2) = [-pi, pi], yrange(2) = [-pi, pi]
real(4), allocatable :: image(:,:,:)

!Problem 3: Scan initial conditions of double pendulum for flip-time
!Get input N for sweep division
numArgs = command_argument_count()
if (numArgs == 1) then
	call get_command_argument(1, arg)
	write (*,*) trim(arg)
	read (arg, *, iostat=status) N
	if (status == 0 .and. (N .ge. 1) .and. N .le. 1000) then
		write (*,*) "Sucessful input: ", N
		allocate(image(1,2*N+1,2*N+1))
	else 
		write (*,*) "Input argument invalid, enter an integer between 1 and 1000"
		call exit
	end if
else
	write (*,*) "Incorrect number of input arguments, refer to readme.md"
	call exit
end if
	
tMax = 1000.0*sqrt(l/gr)
!Do sweep, no zoom, can use symmetry
!write (*,*) "Sweeping initial conditions for flip time, zoom =  0"
!do i=0, N !1-sided sweep, use symmetry to halve run time
!	theta1 = pi * i/N
!	do j=-N, N
!		t = 0.0
!		theta2 = pi * j/N
!		y = [theta1, theta2, 0.0, 0.0]
!		do while ((abs(y(2)) .le. pi) .and. (t .le. tMax))
!			call gl8(y, dt)
!			t = t+dt
!		end do
!		image(1, N+1+i, N+1-j) = t
!		image(1, N+1-i, N+1+j) = t !abusing symmetry
!	end do
!end do
!call write2fits('flipTime0.fits', image, xrange, yrange, ['N'])

!Do sweep, with zoom, can't use symmetry
do zoomPower=1,3
	zoom=3**zoomPower
	write (*,'(x,a,i3)') "Sweeping initial conditions for flip time, zoom =", zoom
	do i=-N, N
		theta1 = (pi * i/N)/zoom + ZoomX
		do j=-N, N
			t = 0.0
			theta2 = (pi * j/N)/zoom + ZoomY
			y = [theta1, theta2, 0.0, 0.0]
			do while ((abs(y(2)) .le. pi) .and. (t .le. tMax))
				call gl8(y, dt)
				t = t+dt
			end do
			image(1, N+1+i, N+1-j) = t
		end do
	end do
	write(filename, fmt='(a,i1,a)') 'flipTime', zoomPower, '.fits'
	call write2fits(filename, image, xrange, yrange, ['N'])
end do

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
