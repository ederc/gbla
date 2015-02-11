#!/bin/sh
make clean;
make converter   OPT="-O2" STA="-static" 
make dump_matrix OPT="-O2" STA="-static" 
mv converter   converter.intel.static
mv dump_matrix dump_matrix.intel.static
make converter   OPT="-O3 -march=native" STA="-static" 
make dump_matrix OPT="-O3 -march=native" STA="-static" 
mv converter   converter.intelavx.static
mv dump_matrix dump_matrix.intelavx.static


