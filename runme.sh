echo "compiling..."
gfortran -c -g -fdefault-real-8 polysvd.f90 -llapack
gfortran -c -g -fdefault-real-8 tester.f90 -llapack
gfortran -g -llapack tester.o polysvd.o -o tester && ./tester > output && cat output
