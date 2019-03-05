echo "compiling..."
gfortran -c -fdefault-real-8 newtonMethod.f90 -llapack
gfortran -c -fdefault-real-8 tester.f90

echo "running..."
gfortran -g --llapack test.o newtonMethod.o -o tester &&./tester > output && cat output
