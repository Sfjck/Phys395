!compile with:
!gfortran -c -fdefault-real-8 tester.f90
!link and run with:
!gfortran -g tester.o gaussLeg.o globalVars.o -o tester && ./tester > output && cat output

program Q3_4
use gaussLeg
use globalVars3_4
implicit none

!Variable declarations
!Time step such that abs(E/E0-1) < 10^-12 for any initial cond
real, parameter :: dt = 0.003, eps = 1.0e-4
real y(2), dydt(2), E0, t, tMax
integer i

!Format labels, x = space, a = chars, nfw.d = n (f)loats (w)idth (d)ecimals
1 format(x, a, ES9.2)
2 format(/, /, I5, a, I5)
3 format(/, x, a)

!Problem 1&2: Integrating equations of motion of double pendulum
write (*,*) "Problem 3: Integrate particle equations of motion"
write (*,*) "Initial conditions: x=a=1, v=0, m=1"
tMax = 20.0
y = [a, 0.0]
dydt = [0.0, 0.0]
call evalf(y, dydt)
E0 = Energy
t = 0

open(unit=1, file="energyViolation.csv", status="replace", action="write")
write(1,*) "t", ",", "|E/E0 -1|"
open(unit=2, file="y.csv", status="replace", action="write")
write(2,*) "t", ",", "x", ",", "v"

do while(t .le. tMax)
	!write data for energy conservation violation
	write(1,*) t, ",", abs(Energy/E0 - 1.0)
	!write data for position/velocity of particle
	write(2,*) t, ",", y(1), ",", y(2)
	
	call gl8(y, dt)
	t = t+dt
end do
close(1)
close(2)

write (*,3) "Problem 4: Find period of partcile with different initial conditions"
write (*,*) "Initial conditions: xi = i/10*a, v=0, m=1"

open(unit=1, file="periods.csv", status="replace", action="write")
write(1,*) "E0", ",", "T"
do i = 1,50
	y = [0.1*i*a, 0.0]
	dydt = [0.0, 0.0]
	call evalf(y, dydt)
	E0 = Energy
	t = 0

	!integrate until partcile reaches "other side", which is half period
	do while(abs(y(1)+0.1*i*a) .ge. eps)
		call gl8(y, dt)
		t = t+dt
	end do
	write(1,*) E0, ",", 2.0*t, ",", y(1)
end do
close(1)

write (*,*) "Data generation for problems 3 and 4 complete."

end program

