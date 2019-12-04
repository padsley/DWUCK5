gfortran -c *.FOR

cd culib4

gfortran -c *.FOR

cd ../culib8

gfortran -c *.FOR

cd ..

gfortran *.o culib8/*.o -o DWUCK5.EXE

./DWUCK5.EXE < DW5TST.DAT
