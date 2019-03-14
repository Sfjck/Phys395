To compile and run (problems 1 & 2):
bash runme.sh
view results in plots.pdf

To animate double pendulum (problem 2, optional):
gnuplot makeAnimation.gpl
to end animation early, spam ctrl+c in the terminal (may vary based on setup)

For problem 3, output of time to flip data is in flipTimeX.fit, where X is level of zoom (1 is no zoom)
To reproduce fit files, use (IC = Initial Condition):
./ICSweep N
Where N is an integer equal to HALF the .fit  image data dimension (eg. N=5 will produce 10x10 image)
Roughly, the time taken is 5 minutes for N=10


**To Do
calc dt such that max error ~1e-6 with theta1=theta2=pi/2
loops for zooms