echo "compile and link"
gfortran -c -fdefault-real-8 ICSweep.f90 -lcfitsio
gfortran -g ICSweep.o gaussLeg.o globalVars.o -o ICSweep -lcfitsio
