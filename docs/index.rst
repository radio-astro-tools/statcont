STATCONT documentation
----------------------

To be implemented.  For now, just see the readme.


The web page (https://hera.ph1.uni-koeln.de/~sanchez/statcont) includes
more complete instructions.

Test data are available `here
<https://hera.ph1.uni-koeln.de/~sanchez/software/STATCONT/test_cases.tar.gz>`_


Interactive Example
-------------------

To reproduce something like figure 13 in `the paper
<http://adsabs.harvard.edu/abs/2018A%26A...609A.101S>`_, you can do the
following (which is based off of the ``test_cases`` data linked above):

.. code:: python

   from astropy.io import fits
   from statcont import cont_finding

   data = fits.getdata('SYNTHETIC_cube.fits')
   sigmaclip_flux_prev, sigmaclip_flux, sigmaclip_noise, filtered_data = cont_finding.c_sigmaclip(data, 1.5, 0)

   import pylab as pl
   pl.plot(data[:,28,43], 'k')
   pl.plot(filtered_data[:,28,43], 'r')

which will show the "filtered data", i.e., the data used for the continuum
measurement, in red overlaid on the black spectrum.
