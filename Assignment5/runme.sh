echo "compiling..."
gfortran -c -fdefault-real-8 globalVars.f90
gfortran -c -fdefault-real-8 integration.f90
gfortran -c -fdefault-real-8 rayleigh.f90 -llapack
gfortran -c -fdefault-real-8 tester.f90

gfortran -g tester.o globalVars.o integration.o rayleigh.o -llapack -o tester 

echo "running... (This will take 2-10minutes)"
./tester > output && cat output

echo "plotting..."
#gnuplot plot1.gpl
#gnuplot plot2.gpl
#gnuplot plot3.gpl
#gnuplot plot4.gpl
gnuplot plot5.gpl

echo "All tasks complete."
