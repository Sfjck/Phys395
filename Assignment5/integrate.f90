!disc
!compile with:
!gfortran -c -fdefault-real-8 integrate.f90

module integrate
use globalVars
implicit none

contains

subroutine integrate(E, init)
	real :: E, init(2)
	real :: x, psi(2), psiMinus(2) !psiMinus = psi(-x), to get odd/even
	character(len=40) :: fileOdd, fileEven
	
	!determine file names, usually I do this in the tester for more modularity 
	!but I'd have to pass the filenames as reference and there's no big difference.
	write(fileOdd, 21) "E=", E, "_Init=(", init(1), "," init(2), ")_odd.dat"
	write(fileEven, 21) "E=", E, "_Init=(", init(1), "," init(2), ")_even.dat"
	open(unit=1, file=fileOdd)
	open(unit=2, file=fileEven)
	
	!using gl8 to integrate
	x = 0.0
	psi = init
	psiMinus = init
	do while ((abs(psi(1)) < 4.0) .and. ((abs(psi(1)) > 0.1) .or. (abs(psi(2)) > 0.1)))
		call gl8(psi, h, x, E)
		call gl8(psiMinus, -1.0*h, -1.0*x, E)
		
		write(1,*), x, (psi-psiMinus)/2.0 !odd parts
		write(1,*), -1.0*x, -1.0*(psi-psiMinus)/2.0
		write(2,*), x, (psi+psiMinus)/2.0 !even parts
		write(2,*), -1.0*x, (psi+psiMinus)/2.0
		x = x+h
	end do
	
	close (1)
	close (2)
	
end subroutine

function V(x)
	real :: V, x
	V = 0.5*m*(w*x)**2.0
end function	
	
! for psi derivative
subroutine evalf(psi, psiP, x, E)
	real :: psi(2), psiP(2), x, E
	psiP(1) = psi(2)
	psiP(2) = psi(1) * (V(x)-E)*2.0*m/(hBar**2.0)
end subroutine evalf

!8th order implicity Gauss-Legendre integrator, modified
subroutine gl8(psi, h, x, E)
	integer, parameter :: s = 4, n = 2
	real psi(n), g(n,s), h, x, E
	integer i, k
	
	! Butcher tableau for 8th order Gauss-Legendre method
	real, parameter :: a(s,s) = reshape((/ &
		0.869637112843634643432659873054998518Q-1, -0.266041800849987933133851304769531093Q-1, &
		0.126274626894047245150568805746180936Q-1, -0.355514968579568315691098184956958860Q-2, &
		0.188118117499868071650685545087171160Q0,   0.163036288715636535656734012694500148Q0,  &
		-0.278804286024708952241511064189974107Q-1, 0.673550059453815551539866908570375889Q-2, &
		0.167191921974188773171133305525295945Q0,   0.353953006033743966537619131807997707Q0,  &
		0.163036288715636535656734012694500148Q0,  -0.141906949311411429641535704761714564Q-1, &
		0.177482572254522611843442956460569292Q0,   0.313445114741868346798411144814382203Q0,  &
		0.352676757516271864626853155865953406Q0,   0.869637112843634643432659873054998518Q-1 /), (/s,s/))
	real, parameter ::   b(s) = (/ &
		0.173927422568726928686531974610999704Q0,   0.326072577431273071313468025389000296Q0,  &
		0.326072577431273071313468025389000296Q0,   0.173927422568726928686531974610999704Q0  /)

	! iterate trial steps
	g = 0.0 
	do k = 1,16
		g = matmul(g,a)
		do i = 1,s
			call evalf(psi + g(:,i)*h, g(:,i), x, E)
		end do
	end do
        
	! update the solution
	psi = psi + matmul(g,b)*h
end subroutine
end module