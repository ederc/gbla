#!/bin/bash

./configure CFLAGS="-march=native -O3" LDFLAGS="-static" ;
make clean ; 
make -j8 ;
gcc -std=gnu99 -fopenmp -march=native -O3 -o gbla -static cli/gbla.o cli/io.o  ./.libs/libgbla.a -lm -lhwloc -fopenmp
mv gbla gbla.intelavx.static ;

./configure CFLAGS="-O2" LDFLAGS="-static" ;
make clean ; 
make -j8 ;
gcc -std=gnu99 -fopenmp -O2 -o gbla -static cli/gbla.o cli/io.o  ./.libs/libgbla.a -lm -lhwloc -fopenmp ;
mv gbla gbla.intel.static ;
# cd tools ;
