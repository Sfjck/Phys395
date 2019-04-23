echo "compiling..."
#problem 1 & 2
gfortran -c -fdefault-real-8 newtonMethod.f90
gfortran -c -fdefault-real-8 bracketSearch.f90
gfortran -c -fdefault-real-8 Q1_2.f90

echo "running..."
#problem 1 & 2
gfortran -g Q1_2.o newtonMethod.o bracketSearch.o -o Q1_2 &&./Q1_2 > output1_2 && cat output1_2

#echo "testing code complete"
#echo "plotting..."
#gnuplot Plot.gpl

#echo "all actions complete"
#echo "outputs can be found at 'Output' and 'Plot.pdf'"
