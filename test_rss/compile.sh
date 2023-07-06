gfortran -c -g MOD_Soilsurface_Resistance.F90
gfortran -c -g test_rss.F90
gfortran -o main test_rss.o MOD_Soilsurface_Resistance.o MOD_Precision.o
./main > 1.txt

