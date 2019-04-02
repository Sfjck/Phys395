!disc
!compile with:
!gfortran -c -fdefault-real-8 integration.f90

module integration
use globalVars
implicit none

contains

!simple bisection
function bisect(Emin, Emax, init)
	real :: Emin, Emax, init(2)
	real :: a, b, c, fa(2), fb(2), fc(2), bisect(2)
	
	a = Emin; fa = fFunc(a, init)
	b = Emax; fb = fFunc(b, init)
	
	! no eigenvalues found
	if (fa(1)*fb(1) > 0) then
		eigenFound = .false.
		call abort
	end if
	eigenFound = .true.
	
	do while (abs(b-a) > eps)
		c = (a+b)/2.0; fc = fFunc(c, init)
		if (abs(fc(1)) < 0.001) exit
		if (fa(1)*fc(1) < 0.0) then
			b = c
			fb = fc
		end if
		if (fb(1)*fc(1) < 0.0) then
			a = c
			fa = fc
		end if
	end do
	
	bisect = [c, fc(2)]
end function 

! basically converts integrate into a function that returns f
! also supresses file outputs
function fFunc(E, init)
	real :: E, init(2), fFunc(2)
	call integrate(E, init, problem = 0); fFunc = f
end function

subroutine integrate(E, init, problem)
	real :: E, init(2)
	integer :: problem
	real :: x, norm
	real :: Psi(2), PsiMinus(2) !PsiMinus = Psi(-x), to get odd/even

	!initialization
	x = 0.0
	Psi = init
	PsiMinus = init
	norm = Psi(1)**2

	!Psi is within 4 and either Psi or dPsi is above 0.1
	do while ((abs(Psi(1)) < 4.0) .and. ((abs(Psi(1)) > 0.1) .or. (abs(Psi(2)) > 0.1)))
		call gl8(Psi, x, dx, E)

		!problem specific code
		!modulo used to skip writing for speed, but still need to compute for accuracy
		select case (problem)
			case (1)
				call gl8(PsiMinus, -1.0*x, -1.0*dx, E) !integrate other side to get Psi(-x)
				if (modulo(x,modx) < dx) then !modulo to skip some points for speed
					write(1,*) x, (Psi-PsiMinus)/2.0 !odd part
					write(2,*) x, (Psi+PsiMinus)/2.0 !even part
				end if
			case (2)
			case (3)
				if (modulo(x,modx) < dx) then
					write(1,*) x, Psi(1), Psi(1)**2/(2.0*norm0)
				end if				
		end select

		x = x + dx
		norm = norm + dx*(Psi(1)**2) !summation for normalization factor
	end do
	f = [Psi(1), norm] !for bisection
end subroutine

!V of harmonic / anharmonic oscillator
function V(x)
	real :: V, x
	if (harmonic .eqv. .true.) then
		V = 0.5*m*(w*x)**2.0
	else
		V = 0.25*lambda*(x**4.0)
	end if
end function	
	
! evaluate derivatives
subroutine evalf(Psi, dPsi, x, E)
	real :: Psi(2), dpsi(2), x, E
	dPsi(1) = Psi(2)
	dPsi(2) = Psi(1) * (V(x)-E)*2.0*m/(hBar**2.0) !2nd derivative of Psi
end subroutine evalf

!8th order implicity Gauss-Legendre integrator, modified
!8th is faster than 10th so...
subroutine gl8(Psi, x, dx, E)
	integer, parameter :: s = 4, n = 2
	real Psi(n), g(n,s), x, dx, E
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
			call evalf(Psi + g(:,i)*dx, g(:,i), x, E)
		end do
	end do
        
	! update the solution
	Psi = Psi + matmul(g,b)*dx
end subroutine
end module
