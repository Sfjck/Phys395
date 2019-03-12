echo "compiling..."
gfortran -c -fdefault-real-8 globalVars.f90
gfortran -c -fdefault-real-8 gaussLeg.f90
gfortran -c -fdefault-real-8 tester.f90
gfortran -c -fdefault-real-8 fractal.f90 -lcfitsio
gfortran -g tester.o gaussLeg.o globalVars.o -o tester 
gfortran -g fractal.o gaussLeg.o globalVars.o -o fractal -lcfitsio

echo "running..."
./tester > output && cat output

echo "plotting..."
gnuplot plot.gpl

echo "plots complete."
