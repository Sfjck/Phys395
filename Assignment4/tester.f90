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
real, parameter :: dt = 0.0005, eps = 1.0e-12
real y(4), dydt(4), E0, t, tMax
logical EViolation

!Format labels, x = space, a = chars, nfw.d = n (f)loats (w)idth (d)ecimals


!Problem 1&2: Integrating equations of motion of double pendulum
write (*,*) "Problem 1&2: Integrate double pendulum equations of motion"
write (*,*) "Initial conditions: theta1 = pi/3, theta2 = -pi/3"
tMax = 100.0 * sqrt(l/gr)
y = [pi/3.0, -pi/3.0, 0.0, 0.0]
dydt = [0.0, 0.0, 0.0, 0.0]
call evalf(y, dydt)
E0 = Energy
EViolation = .false.
t = 0

open(unit=1, file="energyViolation.csv", status="replace", action="write")
write(1,*) "time", ",", "energy conservation violation"
open(unit=2, file="trajectory.csv", status="replace", action="write")
write(2,*) "time", ",", "x", ",", "y"

do while(t .le. tMax)
	!write data for energy conservation violation
	write(1,*) t, ",", abs(Energy/E0 - 1)
!	if (abs(Energy/E0 -1) .ge. eps) then
!		write (*,*) "Energy conservation violation of ", Energy/E0 - 1, " at t = ", t
!		EViolation = .true.
!	end if
	!write data for trajectory of end of pendlulum
!	write(2,*)  l*(sin(y(1))+sin(y(2))), -l*(cos(y(1))+cos(y(2)))

	write(2,*) t, ",", l * (sin(y(1))+sin(y(2))), ",", -l * (cos(y(1))+cos(y(2)))	
	
	call gl8(y, dt)
	t = t+dt
end do

!if (EViolation .eqv. .false.) then
!	write (*,*) "Energy conservation violations stayed below required 1e-12"
!end if
write (*,*) "Integration complete, see plots for energy conservation violation vs time and trajectory of end of pendulum"


!Problem 3:



end program

