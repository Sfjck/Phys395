echo "compiling..."
gfortran -c -fdefault-real-8 globalVars.f90
gfortran -c -fdefault-real-8 gaussLeg.f90
gfortran -c -fdefault-real-8 tester.f90

echo "running..."
gfortran -g tester.o gaussLeg.o globalVars.o -o tester && ./tester > output && cat output
