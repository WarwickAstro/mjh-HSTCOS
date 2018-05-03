# mjh-HSTCOS
Tools for analysis of HST COS data, based on ongoing work.

### Programs:

* lsf-convoluter.py

- - - -

## lsf-convoluter ##

___I intend for this program to become a module once it is more complete and I have the time. Please provide input as you see fit.___

lsf-convoluter convolves a model spectrum (such as a .dk synthetic spectrum) to observed resolution,
using the provided observation wavelengths and Line Spread Function File for the appropriate grating.
With command line argument `p` a plot will be produced upon completing the convolution.

The LSF files can be obtained from the COS STScI webpages, [here](http://www.stsci.edu/hst/cos/performance/spectral_resolution/).

Currently tested for FUV G130M grating (although this should generalise) in lifetime position 3, and briefly position 2.
Support for lifetime position 1 files to be added.

#### Dependancies: #####
numpy, matplotlib, astropy, scipy

#### File Requirements: ####
* LSF file
* Model spectrum
* Observed wavelengths (file can be to the observed data, will read only first column)

Program currently reads input variables from a file in the pwd called `lsf_convoluter_inputs.txt`.
The file should be formatted as lines of `variable_name = variable`.

The following inputs are recognised:
* `lsf_name = path/to/lsf.extension` (REQUIRED)
* `lifetime_pos = 3)` (REQUIRED)
* `spec_name = path/to/spec.extension` (REQUIRED)
* `obs_name = path/to/obs.extension` (REQUIRED)
* `beta = v/c` (optional, if spectrum is to be redshifted before convolution)
* `output_name = path/to/output.extension` (Needed to output convolved spectrum to file)

#### Arguments: ###
`p` --> show plot of original vs convolved spectra
