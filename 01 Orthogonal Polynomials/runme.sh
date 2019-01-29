echo "compiling..."
gfortran -c -fdefault-real-8 gaussjordan.f90
gfortran -c -fdefault-real-8 polyapprox.f90
gfortran -c -fdefault-real-8 tester.f90

echo "running..."
gfortran tester.o gaussjordan.o polyapprox.o -o tester && ./tester > output && cat output
