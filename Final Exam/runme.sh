echo "compiling..."
#problems 1 & 2
gfortran -c -fdefault-real-8 newtonMethod.f90
gfortran -c -fdefault-real-8 bracketSearch.f90
gfortran -c -fdefault-real-8 Q1_2.f90
gfortran -g Q1_2.o newtonMethod.o bracketSearch.o -o Q1_2
#problems 3 & 4
gfortran -c -fdefault-real-8 globalVars3_4.f90
gfortran -c -fdefault-real-8 gaussLeg.f90
gfortran -c -fdefault-real-8 Q3_4.f90
gfortran -g Q3_4.o gaussLeg.o globalVars3_4.o -o Q3_4 
#problem 5
gfortran -c -fdefault-real-8 globalVars5.f90
gfortran -c -fdefault-real-8 rayleigh.f90 -llapack
gfortran -c -fdefault-real-8 Q5.f90
gfortran -g Q5.o globalVars5.o rayleigh.o -llapack -o Q5 

echo "running..."
#problems 1 & 2
./Q1_2 > output1_2 && cat output1_2
#problems 3 & 4
./Q3_4 > output3_4 && cat output3_4
#problem 5
./Q5 > output5 && cat output5

echo "plotting..."
#problems 3 & 4
gnuplot plot3.gpl
gnuplot plot4.gpl

echo "all actions complete"
