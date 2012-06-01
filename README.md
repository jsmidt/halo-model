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

  -v0.1: Working dark matter power spectrum. (eq. 95)
  -v0.2: Working dark matter bispectrum.     (eq. 96)
  -v0.3: Working dark matter trispectrum.    (eq. 97)
  -v0.4: Galaxy-Galaxy and Galaxy-Matter versions of the above.
  -v0.5: General restructuring of code to facilitate next calculations.
  -v0.6: Weak lensing and secondary halo calculations.
  -v0.7: Thermal and Kinetic SZ halo calculations.
  -v0.8: Reionization halo calculations.
  -v0.9: Any more physics I decide to add before v1.0 and final restructuring.
  -v1.0: A final debugged version of the code for the general public.

The hope is to have through versions 0.5 by July and a final verison 1.0 by the
end of the summer with with a few weekly beta releases between version 0.9 and
1.0 to address bugs and restructure to facilitate use for the general user. I
hope to have test programs where results may be easily verified with the
literature. 

Any patches, bug reports and useful comments are welcome so thanks in advance.

Joseph Smidt

