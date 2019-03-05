echo "compiling..."
gfortran -c -fdefault-real-8 newtonMethod.f90 -llapack
gfortran -c -fdefault-real-8 bracketSearch.f90 -llapack
gfortran -c -fdefault-real-8 tester.f90 -llapack

echo "running..."
gfortran -g -llapack tester.o newtonMethod.o bracketSearch.o -o tester &&./tester > output && cat output
