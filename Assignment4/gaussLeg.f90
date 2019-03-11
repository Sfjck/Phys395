!module containing code to integrate equations of motion using 8th order Gauss-Legendre method
!compile with:
!gfortran -c -fdefault-real-8 gaussLeg.f90 -llapack


! Equations of motion:
! th1,2 = theta1,2; p1,2 = p_theta1,2;	m = mass;	l = length; g = gravity
! d = first derivative wrt time
! dth1 = 6/(ml^2) * (2p1 - 3cos(th1-th2)*p2) / (16-9cos^2(th1-th2))
! dth2 = 6/(ml^2) * (8p2 - 3cos(th1-th2)*p1) / (16-9cos^2(th1-th2))
! dp1 = -0.5(ml^2)* (dth1*dth2*sin(th1-th2) + 3(g/l)sin(th1))
! dp2 = -0.5(ml^2)*(-dth1*dth2*sin(th1-th2) + 1(g/l)sin(th2))

module gaussLeg
implicit none
	! Time step such that abs(E/E0-1) < 10^-12 for any initial cond
	real, parameter :: dt = 6.2831853071795864769252867665590057684Q0/2**7
	integer l
	real :: y(4) = (/ 1.0, 0.0, 1.0e-3, 0.0 /)

	do l = 1,2**12
		! output time, position, velocity, error in energy
		write (*,'(6g24.16)') (l-1)*dt, y, y(1)**4/4.0 + y(2)**2/2.0 - 0.25
	
		! step forward with your choice of numerical scheme
		call gl8(y, dt)
	end do
	

contains

! evaluate derivatives, from Andrei
subroutine evalf(y, dydx)
	real y(4), dydx(4)
	real, parameter :: lambda = 10.0
	
	dydx(1) = y(2)
	dydx(2) = -(1.0 + lambda*y(3)**2)*y(1)
	dydx(3) = y(4)
	dydx(4) = -(0.0 + lambda*y(1)**2)*y(3)
end subroutine


!8th order implicity Gauss-Legendre integrator, from Andrei
subroutine gaussLeg8(y, dt)
	integer, parameter :: s = 4, n = 4
	real y(n), g(n,s), dt; integer i, k
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