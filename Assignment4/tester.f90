!compile with:
!gfortran -c -fdefault-real-8 tester.f90
!link and run with:
!gfortran tester.o gaussLeg.o globalVars.o && ./tester > output && cat output

program tester
use gaussLeg
use globalVars
implicit none

!Variable declarations
!Time step such that abs(E/E0-1) < 10^-12 for any initial cond
real, parameter :: dt = 0.0006, eps = 1.0e-12
real y(4), dydt(4), E0, t, tMax

!Format labels, x = space, a = chars, nfw.d = n (f)loats (w)idth (d)ecimals
1 format(x, a, ES9.2)

!Problem 1&2: Integrating equations of motion of double pendulum
write (*,*) "Problem 1&2: Integrate double pendulum equations of motion"
write (*,*) "Initial conditions: theta1=pi/3, theta2=-pi/3, m=1, l=1"
tMax = 100.0 * sqrt(l/gr)
y = [pi/3.0, -pi/3.0, 0.0, 0.0]
dydt = [0.0, 0.0, 0.0, 0.0]
call evalf(y, dydt)
E0 = Energy
t = 0

open(unit=1, file="energyViolation.csv", status="replace", action="write")
write(1,*) "t", ",", "E/E0 -1"
open(unit=2, file="trajectory.csv", status="replace", action="write")
write(2,*) "x", ",", "y"
open(unit=3, file="animation.csv", status="replace", action="write")


do while(t .le. tMax)
	!write data for energy conservation violation
	write(1,*) t, ",", abs(Energy/E0 - 1.0)
	!write data for trajectory of end of pendlulum
	write(2,*) l * (sin(y(1))+sin(y(2))), ",", -l * (cos(y(1))+cos(y(2)))
	!write data for animation
	write(3,*) 0, 0
 	write(3,*) l*sin(y(1)), -l*cos(y(1))
 	write(3,*) l*(sin(y(1))+sin(y(2))), -l*(cos(y(1))+cos(y(2)))
 	write(3,*) " "
 	write(3,*) " "
	
	call gl8(y, dt)
	t = t+dt
end do

write (*,*) "See plots.pdf for plots of energy conservation violation vs time and trajectory of pendulum"

!Problem 3:



end program

