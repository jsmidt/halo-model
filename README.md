halo-model
==========

Halo Model with CAMB and CosmoMC (not yet implemented)

This code implements the halo model described in Cooray and Sheth (2002) astro-ph/0206508 and
calculates the linear matter power spectrum. This code is very preliminary and
still does not work so don't download until completed. 

To run the program check that src/Makefile matches your compiler. The intel
and Gfortran options are already implemented. You also need to have the CAMB folder in the same
level as this project and have compiled it using the option "make all." Make sure you do not use the
compiler option "-fast" as that seems to break the linking. There are some Fortran 2003
features in this code so a modern Fortran compiler is needed (especially for
using CosmoMC). The code has been tested and works on the Gfortran 4.6 and intel 11.1 compilers. 

Then to run type:

make

./halo

This will calculate the dark matter power spectrum with output in the output
directory. If you have python and matplotlib, a
quick plot of the results are given py typing:

python tplot.py

* This is considered alpha-quality code use at own risk *

Any patches, bug reports and useful comments are welcome so thanks in advance.
