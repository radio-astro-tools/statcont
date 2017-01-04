#!/usr/bin/env python

from __future__ import print_function
import argparse

from .subcont import process_files


def main(args=None):

    # see http://docs.astropy.org/en/stable/development/scripts.html

    ################################################################################
    #
    # Creating the list of options that can be used in this script:
    #

    pars = argparse.ArgumentParser(description="Continuum-subtraction script")

    grou = pars.add_mutually_exclusive_group()
    grou.add_argument('-i', '--iname', nargs='*', help='NECESSARY: unless parameters -f or -s are considered. \
                                                        One or more fits files, stored in the directory --path \
                                                        Name of the files without the extension [ .fits ]')
    pars.add_argument('-m', '--imerge', nargs=1,  help='OPTIONAL: Files in --iname are merged \
                                                        The argument is the output file name')
    grou.add_argument('-f', '--ifile', nargs=1,   help='NECESSARY: unless parameters -i or -s are considered. \
                                                        File with the names, in one column, of the fits files. \
                                                        Name of the files without the extension [ .fits ]')
    grou.add_argument('-s', '--ispec', nargs=1,   help='NECESSARY: unless parameters -i or -f are considered. \
                                                        One single ascii-format file with two columns: \
                                                        frequency (c.1) and intensity (c.2), and no header. \
                                                        Name of the files without the extension [ .dat ]')
    pars.add_argument('-p', '--ipath', nargs=1,   help='OPTIONAL: Specify the path where there original files \
                                                        are saved, within the data directory')
    pars.add_argument('-n', '--noise', nargs=1, type=float,   help='NECESSARY: RMS noise level of the observataions')
    pars.add_argument('-v', '--verbose', action='store_true', help='OPTIONAL: Increase output verbosity')
    pars.add_argument('-c', '--cutout', nargs=3, type=int,    help='OPTIONAL: Create a cutout image of the original file. \
                                                                    Three integer numbers are required for this option: \
                                                                    xcen: x coordinate of the central pixel of the cutout \
                                                                    ycen: y coordinate of the central pixel of the cutout \
                                                                    size: size (in pixels) of the cutout')
    pars.add_argument('--plots', action='store_true',         help='OPTIONAL: Create plots on a pixel-by-pixel basis \
                                                                    (computing time increases considerably)')
    pars.add_argument('--continuum', action='store_true',     help='Determination of the continuum level (SIGMACLIP)')
    pars.add_argument('--cmax', action='store_true',          help='Continuum maps using method MAX')
    pars.add_argument('--cmean', action='store_true',         help='Continuum maps using method MEAN and MEANSEL')
    pars.add_argument('--cmedian', action='store_true',       help='Continuum maps using method MEDIAN and MEDIANSEL')
    pars.add_argument('--cpercent', action='store_true',      help='Continuum maps using method PERCENTILE (25, 50 or 75)')
    pars.add_argument('--cKDEmax', action='store_true',       help='Continuum maps using method KDEMAX')
    pars.add_argument('--cGaussian', action='store_true',     help='Continuum maps using method GAUSSIAN and GAUSSNW')
    pars.add_argument('--csigmaclip', action='store_true',    help='Continuum maps using method SIGMACLIP')
    pars.add_argument('--call', action='store_true',          help='Continuum maps using ALL the methods')
    pars.add_argument('--cfree', action='store_true',         help='Remove the continuum to the original datacube')
    pars.add_argument('--spindex', action='store_true',       help='Fit a linear equation to the continuum of different files')
    op = pars.parse_args(args)

    if op.call:
        op.cmax = True
        op.cmean = True
        op.cmedian = True
        op.cpercent = True
        op.cKDEmax = True
        op.cGaussian = True
        op.csigmaclip = True
        
    if op.imerge:
        op.continuum = True

    if op.continuum:
        op.csigmaclip = True
        op.cfree = True

    if op.cmax or op.cmean or op.cmedian or op.cpercent or op.cKDEmax or op.cGaussian:
        op.csigmaclip = True

    if op.cmax or op.cmean or op.cmedian or op.cpercent or op.cKDEmax or op.cGaussian or op.csigmaclip:
        op.continuum = True



    # Noise level of your data cubes (in units of the fits files)
    rms_noise = op.noise[0]

    process_files(tmp_files=op.tmp_files, ispec=op.ispec, iname=op.iname,
                  ifile=op.ifile, rms_noise=rms_noise, continuum=op.continuum, cmax=op.cmax,
                  cmean=op.cmean, cmedian=op.cmedian, cpercent=op.cpercent,
                  cKDEmax=op.cKDEmax, cGaussian=op.cGaussian, cfree=op.cfree,
                  csigmaclip=op.csigmaclip, verbose=op.verbose, plots=op.plots,
                  imerge=op.imerge, cutout=op.cutout)
