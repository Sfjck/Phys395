echo "compiling..."
gfortran -c -g -fdefault-real-8 polyFit.f90 -llapack
gfortran -c -g -fdefault-real-8 tester.f90 -llapack

echo "running..."
gfortran -g -llapack tester.o polyFit.o -o tester && ./tester > output && cat output

echo "testing code has completed."
echo "plotting results..."
gnuplot Plots.gpl

echo "all actions have completed."
echo "outputs can be found at 'Output' and 'Plots.pdf'"
