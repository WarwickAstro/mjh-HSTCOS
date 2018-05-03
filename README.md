# mjh-HSTCOS
Tools for analysis of HST COS data, based on ongoing work.

### Programs:

* lsf-convoluter.py

- - - -

## lsf-convoluter ##

lsf-convoluter convolves a model spectrum (such as a .dk synthetic spectrum) to observed resolution,
using the provided observation wavelengths and Line Spread Function File for the appropraite grating.

The LSF files can be obtained from the COS STScI webpages, [here](http://www.stsci.edu/hst/cos/performance/spectral_resolution/).

Currently tested for FUV gratings in lifetime positions 2 and 3.

#### Dependancies: #####
NumPy, matplotlib, astropy, scipy

#### File Requirements: ####
* LSF file
* Model spectrum
* Observed wavelengths (file can be to the observed data, will read only first column)


