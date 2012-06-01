halo-model
==========

Halo Model with CAMB

This code implements the halo model assuming to have access to CAMB to
calculate the linear matter power spectrum. This code is very preliminary and
still does not work so don't download until completed.  For more details of the equations refer to Cooray and Sheth (2002) astro-ph/0206508.

To run the program check that src/Makefile matches your compiler.  The intel
and Gfortran options are already implemented. There are some Fortran 2003
features in this code so a modern Fortran compiler is needed. The code has been
tested and works on the Gfortran 4.6 and intel 12.1 compilers.  

Then to run type:

make
./halo

to calculate the dark matter power spectrum.  The output will be in the output
directory which you may need to create. If you have python and matplotlib a
quick plot of the results are given py typing:

python tplot.py

* This is considered alpha-quality code use at own risk *

The release schedule is expected to be a few betas followed by a release
candidate followed by a version 1.0 release.  The goals for each release 
are:

  -v0.6: (Beta) Working dark matter power spectrum. (eq. 95)
  -v0.7: (Beta) Working dark matter bispectrum.     (eq. 96)
  -v0.8: (Beta) Working dark matter trispectrum.    (eq. 97)
  -v0.9: (Release candidate) De-bugged & well-structured version of the above.
  -v1.0: Bug free release.  

The hope is to have the version 1.0 out sometime in July.  A version 2.0 will
quickley be in the works implementing the correlations between matter and
galaxies as well as various angular spectra which are projections on the 2D sky.

