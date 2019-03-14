**Problems 1 & 2, compile, run, and plot:
bash runme.sh
View results in plots.pdf

**Problem 2 optional, make animation:
gnuplot makeAnimation.gpl
To end animation early, spam ctrl+c in the terminal (may vary based on setup)

**Problem 3, output of time to flip data is in flipTimeX.fits, where 3^X is level of zoom (0 is no zoom). To reproduce fits files:
./ICSweep N
Replace N with 1 integer between 1 and 500 (inclusive)
Where 2N+1 is the dimension of the .fits data
(eg N=10 will produce an 21x21 image)
The full sweep may take several minutes to hours to complete, depending on N
Use N = 10 for a low res, fast sweep (~5min)
