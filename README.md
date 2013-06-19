halo-model
==========

Halo Model with CAMB and CosmoMC (not yet implemented)

Based on Joseph Smidt's code, this code implements the halo model assuming to have access to CAMB to
calculate the linear matter power spectrum. This code is very preliminary and
still does not work so don't download until completed.  For more details of the equations refer to Cooray and Sheth (2002) astro-ph/0206508.

To run the program check that src/Makefile matches your compiler.  The intel
and Gfortran options are already implemented. There are some Fortran 2003
features in this code so a modern Fortran compiler is needed (especially for
using CosmoMC). The code has been tested and works on the Gfortran 4.6 and intel 12.1 compilers. 

Then to run type:

make

./halo

to calculate the dark matter power spectrum.  The output will be in the output
directory which you may need to create. If you have python and matplotlib, a
quick plot of the results are given py typing:

python tplot.py

* This is considered alpha-quality code use at own risk *

Any patches, bug reports and useful comments are welcome so thanks in advance.

Cameron Thacker
