import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy import interpolate

"""
Convolve a HST/COS Line Spread Function with a chosen input spectrum

This program takes a spectrum (typically a .dk file),
observed wavelengths,
and a COS Line Spread Function file, available at
http://www.stsci.edu/hst/cos/performance/spectral_resolution/

This intended for synthetic spectra to be
reproduced at observing resolution.

The spectrum is linearly interpolated to observed wavelengths.

A convolution grid is produced from the LSF file,
and manually convolved with the interpolated spectrum,
resulting in a synthetic spectrum of typical observing resolution.

Inputs:
    Required:
        LSF file name
        Spectrum data file name
        Observed data file name (for wavelengths)
        Lifetime position
    Optional:
        Radial velocity
        Output file destination

Outputs:
Convolved spectrum

"""

#c = 2.99792458e5 # speed of light in km/s


def apply_redshift(beta, wavs, plot=0):
    """Redshift spec by beta = v/c"""
    if plot:
        plt.figure()
        plt.plot(wavs, c='red')
    wavs += wavs * (np.sqrt(1+beta)/np.sqrt(1-beta) - 1)
    if plot:
        plt.plot(wavs, c='green')
        plt.show()
    return wavs


def fold_tails(flux_toolong, lsf_size):
    """Fold flux_toolong, to reflect tails into edge data
    
    flux_toolong = unfolded array
    lsf_size = size of line spread function ie number of pixels
    (should be the same for all lsfs if multiple used)
    
    Returns:
    flux_justright = tail-folded 1d flux array"""
    fl = lsf_size//2
    # separate tails of conv procedure
    lead_flux = flux_toolong[:fl]
    follow_flux = flux_toolong[-1*fl:]
    #follow_flux[:] = 0
    flux_justright = flux_toolong[fl:-1*fl]
    # add flipped tails
    flux_justright[:fl] += np.flip(lead_flux, axis=0)
    flux_justright[-1*fl:] += np.flip(follow_flux, axis=0)
    
    return flux_justright


def model_x_lsf(specflx, lsf):
    """Perform convolution of model with lsf profiles
    
    mflx = model flux array (1d)
    lsf = line spread functions (2d array; 'stack' of lsfs interpolated from file)"""
    
    if lsf.shape[0] != specflx.shape[0]:
        if lsf.shape[1] == specflx.shape[0]:
            lsf = np.transpose(lsf)
        else:
            print("LSF array length not equal to flux array length! Exiting...")
            sys.exit(1)
    
    # define lengths required, and lsf holder array
    lwavs, lpix = lsf.shape
    if lpix%2 == 0:
        print("Detected pixel width of LSF is not odd!\nExiting...")
        sys.exit(1)
    holder = np.zeros((lwavs, lwavs+lpix-1))
    
    # multiply each lsf by relevant flux
    fluxed_lsf = specflx * lsf
    # write into holder, along lead diagonal
    for l in range(lwavs):
        holder[l,l:l+lpix] = fluxed_lsf[l]

    # and sum holder along pix dimension, giving a 1d flux array with tails
    wide_flux = np.sum(holder, axis=0)
    
    #return wide_flux, lpix
    # fold tails in, to account for flux loss
    convolved_flux = fold_tails(wide_flux, lpix)
    
    return convolved_flux


def getinput_file():
    """Read parameters in input file into control dictionary
    """
    dic = {}
    try:
        with open('lsf_convoluter_inputs.txt') as inp_file:
            inputs = inp_file.readlines()
        for line in inputs:
            # discard blank lines
            if not line.strip():
                continue
            param, val = line.split('=')
            dic[param.strip().lower()] = val.strip()
        return dic
    except OSError:
        return OSError


def getlsf(lsf_filename, lp):
    """Read lsf file (lsfname), which is in format of lp2 onwards
    
    Returns:
    wavelengths  (shape n,)
    spread_functions (shape n,pix)
    """
    
    if lp == 2 or lp == 3:
        try:
            lsf_file = np.transpose(np.genfromtxt(lsf_filename, delimiter=''))
        except OSError:
            print('File {} not found! Exiting...'.format(lsf_filename))
            sys.exit(1)
        wavls = np.copy(lsf_file[:,0])
        spread_funcs = np.copy(lsf_file[:,1:])
        return wavls, spread_funcs
    elif lp == 1:
        print('No code for LPPOS=1 yet!')
        sys.exit(1)
    else:
        print('LPPOS not defined or invalid ({})'.format(lifetime_pos))
        sys.exit(1)


def getspec_file(sname, dlim=''):
    """Retrieve data from spectrum file
    
    Does not accomodate errors on synthetic spectra"""
    try:
        spec_data = np.genfromtxt(sname, delimiter=dlim, usecols=(0,1))
    except OSError:
        print('File {} not found!\nExiting...'.format(sname))
        return OSError
    return spec_data


def getobs_wavelengths(obsname):
    """Retrieve column 0 of observations file
    Requires wavelengths are first column of input file"""
    try:
        obs_wavelengths = np.genfromtxt(obsname, usecols=(0))
    except OSError:
        print('File {} not found!\nExiting...'.format(obsname))
        return NameError
    return obs_wavelengths


def main():
    # arguements
    try:
        cmdline = sys.argv[1]
    except IndexError:
        cmdline = ""
    
    try:
        INPUT_DICT = getinput_file()
    except OSError:
        print('lsf_conv_inputs.txt file deos not exist!\nExiting...')
        return
    
    # open LSF file
    try:
        INPUT_DICT['lsf_name']
    except KeyError:
        print("No input LSF file given!\nExiting...")
        return
    try:
        INPUT_DICT['lifetime_pos']
    except KeyError:
        print("No input lifetime position given!\nExiting...")
        return
    lsf_wavs, lsf_spreads = getlsf(INPUT_DICT['lsf_name'], int(INPUT_DICT['lifetime_pos']))
    
    # open spectrum file
    try:
        spec = getspec_file(INPUT_DICT['spec_name'])
    except KeyError:
        print("No input spectrum file given!\nExiting...")
        return
    # if beta given, apply redshift
    try:
        beta = float(INPUT_DICT['beta'])
        # if beta > 1 or < -1
        spec = apply_redshift(beta, spec)
    except KeyError:
        pass
    
    # open wavelengths file
    try:
        wavelengths = getobs_wavelengths(INPUT_DICT['obs_name'])
    except:
        print("No input wavelengths file given\nExiting...")
        return INPUT_DICT
    # interpolate model to observed wavelengths (and make 2d for conv)
    fluxes = np.interp(wavelengths, spec[:,0], spec[:,1])[np.newaxis].T
    
    # set up pixel dimension for convolution
    pix = np.arange(lsf_spreads.shape[1])
    # create 2d interpolation grid for LSF
    lsf_interp = interpolate.interp2d(pix, lsf_wavs, lsf_spreads)
    # interpolate onto grid of observed wavelength * pixel
    lsf_vals = lsf_interp(pix, wavelengths)
    # integrate each lsf and normalise
    lsf_i = np.trapz(lsf_vals, pix, axis=1).reshape(len(lsf_vals), 1)
    lsf_vals /= lsf_i
    
    #return fluxes, lsf_vals
    # convolve, then stack result with wavelengths
    conv_spec = model_x_lsf(fluxes, lsf_vals)
    
    output_spec = np.stack((wavelengths, conv_spec), axis=1)
    
    # plot and save

    if 'p' in cmdline.lower():
        f = plt.figure()
        a = plt.gca()
        a.plot(spec[:,0], spec[:,1], lw=2, c='black', label='Original Spectrum')
        a.plot(output_spec[:,0], output_spec[:,1], lw=1, c='red', label='Convolved Spectrum')
        a.set_xlim(output_spec[:,0].min()*0.95, output_spec[:,0].max()*1.05)
        a.set_ylim(output_spec[:,1].min()*0.95, output_spec[:,1].max()*1.05)
        f.legend()
        f.show()
    
    try:
        if INPUT_DICT['output_name']:
            np.savetxt(INPUT_DICT['output_name'], output_spec, delimiter=' ')
            print("Convolved spectrum saved at {}".format(INPUT_DICT['output_name']))
    except NameError:
        print("No output name given - File not saved.")
    


if __name__ == "__main__":
    print("\nBeginning convolution program.\n")
    main()
    print("\nConvolution program complete!")
    
