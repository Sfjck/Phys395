echo "compiling..."
gfortran -c -fdefault-real-8 polysvd.f90 -llapack
gfortran -c -fdefault-real-8 tester.f90 -llapack
gfortran -llapack tester.o polysvd.o -o tester && ./tester > output && cat output
