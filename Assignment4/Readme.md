To compile and run (problems 1 & 2):
bash runme.sh
view results in plots.pdf

To animate double pendulum (problem 2, optional):
gnuplot makeAnimation.gpl
to end animation early, spam ctrl+c in the terminal (may vary based on setup)

For problem 3, output of time to flip data is in flipTimeX.fit, where 3^X is level of zoom (0 is no zoom)
To reproduce fit files, use (IC = Initial Condition):
./ICSweep N
Where 2N+1 is the dimension of the .fit image
(eg N=10 will produce an 11x11 image)
The full sweep may take several minutes to hours to complete, depending on N
Use N = 10 for low res, fast sweep (~5min)


**To Do
loops for zooms
