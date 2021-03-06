!module containing code to integrate equations of motion using 8th order Gauss-Legendre method
!compile with:
!gfortran -c -fdefault-real-8 gaussLeg.f90

module gaussLeg
use globalVars3_4
implicit none

contains

! evaluate derivatives & calculate energies for particle using equations of motion
!y(1) = x, y(2) = v
subroutine evalf(y, dydt)
	real y(2), dydt(2); intent(in) y
	
	dydt(1) = y(2)
	dydt(2) = -2.0*V0/(m*a) * sinh(y(1)/a)/(cosh(y(1)/a))**3.0
	
	Energy = m/2.0*y(2)**2.0 - V0/(cosh(y(1)/a))**2.0
end subroutine


!8th order implicity Gauss-Legendre integrator, from class
subroutine gl8(y, dt)
	integer, parameter :: s = 4, n = 2
	real y(n), g(n,s), dt
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
			call evalf(y + g(:,i)*dt, g(:,i))
		end do
	end do
        
	! update the solution
	y = y + matmul(g,b)*dt
end subroutine
end module
