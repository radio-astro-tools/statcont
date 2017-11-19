#!/usr/bin/env python

# See http://docs.astropy.org/en/stable/development/scripts.html

""" Description:

-----------------------------------------------------------------
 STATCONT: A statistical continuum level determination method for
           line-rich sources"

 Sanchez-Monge et al (2017, A&A, submitted)
 Details in https://hera.ph1.uni-koeln.de/~sanchez/statcont
 version 1.2

"""

from __future__ import print_function
import argparse

from .statcont_main import process_files

def main(args=None):

    # Creating the list of options
    pars = argparse.ArgumentParser(description="+++ ----------------------------------------------------------------------- \
                                                +++ STATCONT : A statistical continuum level determination method (v 1.2) ")
    grou = pars.add_mutually_exclusive_group()
    grou.add_argument('-i', '--iname', nargs='*',
                      help='NECESSARY: unless parameters -f or -s are considered. \
                            One or more FITS files, stored in the directory --path \
                            Name of the files without the extension [ .fits ]')
    grou.add_argument('-f', '--ifile', nargs=1,
                      help='NECESSARY: unless parameters -i or -s are considered. \
                            File containing a 1-column list with the file names. \
                            Name of the files without the extension [ .fits ]')
    grou.add_argument('-s', '--ispec', nargs='*',
                      help='NECESSARY: unless parameters -i or -f are considered. \
                            One single ASCII-format file with two columns: \
                            frequency (c.1) and intensity (c.2), and no header. \
                            Name of the files without the extension [ .dat ]')
    pars.add_argument('-m', '--imerge', nargs=1,
                      help='OPTIONAL: All the files in --iname are merged \
                            The argument is the output name for the merged file')
    pars.add_argument('-p', '--ipath', nargs=1,
                      help='OPTIONAL: Specify the PATH where there original files \
                            are saved, within the DATA directory')
    pars.add_argument('-n', '--noise', nargs=1, type=float,
                      help='NECESSARY: RMS noise level of the observataions')
    pars.add_argument('--continuum', action='store_true',
                      help='OPTIONAL: Determination of the continuum (and noise) \
                            Method corrected SIGMACLIP (fast) is used by default\
                            Subtraction of continuum to line data (--cfree)')
    pars.add_argument('--cmax', action='store_true',
                      help='OPTIONAL: Continuum from the MAXIMUM of the histogram')
    pars.add_argument('--cmean', action='store_true',
                      help='OPTIOANL: Continuum from the MEAN of the intensties')
    pars.add_argument('--cmedian', action='store_true',
                      help='OPTIONAL: Continuum from the MEDIAN of the intensities')
    pars.add_argument('--cpercent', action='store_true',
                      help='OPTIONAL: Continuum from the 25th and 75th PERCENTILE')
    pars.add_argument('--cKDEmax', action='store_true',
                      help='OPTIONAL: Continuum from the maximum of the KDE')
    pars.add_argument('--cGaussian', action='store_true',
                      help='OPTIONAL: Continuum from a Gaussian fit to the histogram')
    pars.add_argument('--csigmaclip', action='store_true',
                      help='OPTIONAL: Continuum from the SIGMA-CLIPPING algorithm')
    pars.add_argument('--call', action='store_true',
                      help='OPTIONAL: Continuum using ALL the methods')
    pars.add_argument('--cfree', action='store_true',
                      help='OPTIONAL: Remove the continuum to the original datacube')
    pars.add_argument('--nooffset', action='store_true',
                      help='OPTIONAL: Search for and remove a general continuum offset')
    pars.add_argument('--spindex', action='store_true',
                      help='OPTIONAL: Determine the spectral index (ALPHA), as \
                            a linear function: Flux = FACTOR * frequency ^ (ALPHA) \
                            Provide a list of files to be processed, for which \
                            the continuum needs to be created (--continuum)')
    pars.add_argument('--model', action='store_true',
                      help='OPTIONAL: Create a continuum model using information \
                            from --spindex, and a line+continuum(model) cube')
    pars.add_argument('--cutout', nargs=3, type=int,
                      help='OPTIONAL: Create a cutout image of the original file. \
                            Three integer numbers are required for this option: \
                            xcen: x coordinate of the central pixel of the cutout \
                            ycen: y coordinate of the central pixel of the cutout \
                            size: size (in pixels) of the cutout')
    pars.add_argument('--plots', action='store_true',
                      help='OPTIONAL: Create plots on a pixel-by-pixel basis \
                            Spectrum with continuum (and noise) levels indicated \
                            (computing time increases considerably)')
    pars.add_argument('--verbose', nargs='?', type=int, const=1, default=2,
                      help='OPTIONAL: Modify level of output verbosity (default 2) \
                            0: Switch off from terminal all message communications \
                            1: Only important messages are displayed in terminal \
                            2: Standard level of output verbosity (default value) \
                            3: Increase level of output verbosity (more details) \
                            4: Switch on all messages (for developers)')
    pars.add_argument('--betaversion', action='store_true',
                      help='DEVELOPERS: For developers and testing')
    op = pars.parse_args(args)
    
    # Modify level of output verbosity
    verbose = op.verbose

    # Activate level 2 of verbosity for ASCII spectra
    # Can be overwritten by specifying the verbosity output level
    if op.ispec and verbose == 2:
        verbose = 3
    
    # For developers, activate maximum verbosity
    # Can be overwritten by specifying the verbosity output level
    if op.betaversion and verbose == 2:
        verbose = 4

    # Select all continuum determination methods when --call is used
    if op.call:
        op.cmax = True
        op.cmean = True
        op.cmedian = True
        op.cpercent = True
        op.cKDEmax = True
        op.cGaussian = True
        op.csigmaclip = True

    # Activate the main code for continuum determination when using any method
    if op.cmax or op.cmean or op.cmedian or op.cpercent or op.cKDEmax or op.cGaussian or op.csigmaclip:
        op.continuum = True

    # Use the SIGMA-CLIPPING method and subtract continuum to line, if --continuum is used
    if op.continuum:
        op.cfree = True

    # Activate the subtraction of the continuum --cfree, if the general offset wants to be removed
    if op.nooffset:
        op.continuum = True
        op.cfree = True

    # Activate the determination of the spectral index --spindex, if the model wants to be created
    if op.model:
        op.spindex = True
            
    # Noise level of your data cubes (in units of the FITS file)
    rms_noise = op.noise[0]
    
    process_files(iname=op.iname,
                  ifile=op.ifile,
                  ispec=op.ispec,
                  imerge=op.imerge,
                  ipath=op.ipath,
                  rms_noise=rms_noise,
                  continuum=op.continuum,
                  cmax=op.cmax,
                  cmean=op.cmean,
                  cmedian=op.cmedian,
                  cpercent=op.cpercent,
                  cKDEmax=op.cKDEmax,
                  cGaussian=op.cGaussian,
                  csigmaclip=op.csigmaclip,
                  cfree=op.cfree,
                  nooffset=op.nooffset,
                  spindex=op.spindex,
                  model=op.model,
                  plots=op.plots,
                  cutout=op.cutout,
                  verbose=verbose)
