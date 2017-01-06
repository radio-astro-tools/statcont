################################################################################
#
# SUBCONT: A statistical continuum level determination method
#
################################################################################

"""

----------------------------------------------------------------"
 SUBCONT: A statistical continuum level determination method for"
          line-rich sources"
 "
 Sanchez-Monge et al (2017, A&A, submitted)              v 1.0.0"
 "

"""


################################################################################
#
# Importing packages and general setup:
#

from __future__ import print_function
import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import astropy
import astropy.io.ascii as ascii
from astropy.io import fits
from astropy.stats import sigma_clip
from scipy import stats
from scipy.optimize import leastsq

from . import fits_cutout
# not used from . import kdestats



################################################################################
#
# Variables and functions used throughout the code:
#

# Definition of the Gaussian function
fitfunc = lambda p, x: p[0]*np.exp(-0.5*((x-p[1])/p[2])**2.)
errfunc = lambda p, x, y: (y - fitfunc(p, x))

# Definition of the Gaussian function (version 2 of the fit, currently not used)
def gaussian(x, a, mu, sigma):
    return a * np.exp(-(1./2.) * np.power((x - mu)/sigma, 2.))


################################################################################
#
# Main code: Reading data and header, and determining the continuum level:
#

def process_files(ispec=False, iname=False, ifile=False, ipath=False,
                  rms_noise=None, continuum=False, cmax=False, cmean=False,
                  cmedian=False, cpercent=False, cKDEmax=False,
                  cGaussian=False, cfree=False, csigmaclip=False, cutout=False,
                  spindex=False, verbose=False, plots=False, imerge=False):

    ################################################################################
    #
    # List of files/datacubes to be analyzed:
    #

    if iname:
        input_files = iname
        extension = '.fits'

    elif ifile:
        print(ifile[0])
        lines = [line.rstrip('\n') for line in open(ifile[0])]
        input_files = lines
        extension = '.fits'

    elif ispec:
        input_files = ispec
        extension = '.dat'
        verbose = True



    ################################################################################
    #
    # Definition of directories:
    #

    os.system('mkdir -p data/')
    os.system('mkdir -p products/')

    data_path = "data/"
    cont_path = "products/"
    line_path = cont_path

    if ipath:
        source = ipath[0]
        sourcedir = source + '/'
        data_path = data_path + sourcedir
        cutout_path = data_path + 'cutout/'
        os.system('mkdir -p ' + cutout_path)
        merged_path = data_path + 'merged/'
        os.system('mkdir -p ' + merged_path)
        cont_path = cont_path + sourcedir
        os.system('mkdir -p ' + cont_path)
        line_path = line_path + sourcedir
        os.system('mkdir -p ' + line_path)

    plots_path = cont_path + 'plots/'
    os.system('mkdir -p ' + plots_path)




    ################################################################################
    #
    # Setting the path and file names
    # ... and in case of fits files, using cutout to create smaller files:
    #

    if ispec:
        tmp_files = []
        for file_name in input_files:
            tmp_path = data_path
            tmp_file = file_name
            tmp_files.append(tmp_file)

    if iname or ifile:
        if cutout:
            print("+++ Producing cutout files ...")
            tmp_files = []
            for file_name in input_files:
                data_fitsfile = data_path + file_name + extension
                central_xpixel = cutout[0]
                central_ypixel = cutout[1]
                number_pixels = cutout[2]
                if verbose:
                    print("  . Cutout of %s \n    at central pixel %i, %i with size %i" % (data_fitsfile, central_xpixel, central_ypixel, number_pixels))
                cutout_fitsfile = cutout_path + file_name + '_cutout' + extension
                fits_cutout.cutout(data_fitsfile, central_xpixel, central_ypixel, number_pixels, cutout_fitsfile)
                tmp_path = cutout_path
                tmp_file = file_name + '_cutout'
                tmp_files.append(tmp_file)
        else:
            tmp_files = []
            for file_name in input_files:
                tmp_path = data_path
                tmp_file = file_name
                tmp_files.append(tmp_file)



    ################################################################################
    #
    # Merging the files, in case this is required:
    #

    if imerge:
        print("+++ Merging files ...")
        
        #
        # name of the output file
        merged_file_name = imerge[0]

        #
        # defining ranges for channels
        ninterval = len(tmp_files)
        bnchan = np.empty([ninterval])
        enchan = np.empty([ninterval])

        #
        # total number of channels of merged file
        icount = 0
        for tmp_file in tmp_files:
            
            header = fits.getheader(tmp_path + tmp_file + extension)
            ndim = header.get('NAXIS')
            nxpix = header.get('NAXIS1')
            nypix = header.get('NAXIS2')
            nchan = header.get('NAXIS3')
            npolz = header.get('NAXIS4')
            
            if icount == 0:
                bnchan[icount] = 0
                enchan[icount] = nchan
            else:
                jcount = icount-1
                bnchan[icount] = enchan[jcount]
                enchan[icount] = nchan+enchan[jcount]
                
            icount = icount + 1

        #
        # creating output merged file
        merged_data = np.empty([npolz,enchan[icount-1],nypix,nxpix])
        
        #
        # writing the data on the new merged file
        icount = 0
        for tmp_file in tmp_files:
            
            data = fits.getdata(tmp_path + tmp_file + extension)
            merged_data[0,bnchan[icount]:enchan[icount],:,:] = data
            icount = icount + 1

        merged_file = merged_path + merged_file_name + extension
        os.system('rm -rf ' + merged_file)
        fits.writeto(merged_file, np.float32(merged_data), header=header, clobber=True)
        
        unmerged_files = tmp_files
        unmerged_path = tmp_path
        tmp_files = []
        tmp_files.append(merged_file_name)
        tmp_path = merged_path


    #
    # loop through all the files that will be processed
    for tmp_file in tmp_files:

        print("")
        print("+++ PROCESSING " + tmp_path + tmp_file + extension)

        #
        # reading data and header of the fits file
        if iname or ifile:

            header = fits.getheader(tmp_path + tmp_file + extension)
            data = fits.getdata(tmp_path + tmp_file + extension)

            ndim = header.get('NAXIS')
            nxpix = header.get('NAXIS1')
            nypix = header.get('NAXIS2')
            nchan = header.get('NAXIS3')
            # not used npolz = header.get('NAXIS4')
            bunit = header.get('BUNIT')

        #
        # reading the data and header of the ascii/spectrum file
        if ispec:

            specfile = open(tmp_path + tmp_file + extension, 'r')

            nxpix = 1
            nypix = 1
            nchan = 1
            # not used npolz = 1

        #
        # loop through the y and x pixels to determine the continuum level
        if continuum:

            #
            # setting up the variables that will contain
            # the continuum level and the uncertainty (dispersion)
            if cmax:
                continuum_flux_maximum = []
            if cmean:
                continuum_flux_mean = []
                continuum_flux_meansel = []
            if cmedian:
                continuum_flux_median = []
                continuum_flux_mediansel = []
            if cpercent:
                continuum_flux_percent25 = []
                continuum_flux_percent50 = []
                continuum_flux_percent75 = []
            if cKDEmax:
                continuum_flux_KDEmax = []
            if cGaussian:
                continuum_flux_GaussNw = []
                continuum_noise_GaussNw = []
                continuum_flux_Gaussian = []
                continuum_noise_Gaussian = []
            if csigmaclip:
                continuum_flux_sigmaclip = []
                continuum_noise_sigmaclip = []

            #
            # loop through y pixels
            for ypix in range(nypix):

                print("... analyzing column " + str(ypix+1) + " out of " + str(nypix))

                if cmax:
                    continuum_flux_maximum.append([])
                if cmean:
                    continuum_flux_mean.append([])
                    continuum_flux_meansel.append([])
                if cmedian:
                    continuum_flux_median.append([])
                    continuum_flux_mediansel.append([])
                if cpercent:
                    continuum_flux_percent25.append([])
                    continuum_flux_percent50.append([])
                    continuum_flux_percent75.append([])
                if cKDEmax:
                    continuum_flux_KDEmax.append([])
                if cGaussian:
                    continuum_flux_GaussNw.append([])
                    continuum_noise_GaussNw.append([])
                    continuum_flux_Gaussian.append([])
                    continuum_noise_Gaussian.append([])
                if csigmaclip:
                    continuum_flux_sigmaclip.append([])
                    continuum_noise_sigmaclip.append([])

                #
                # loop through x pixels
                for xpix in range(nxpix):

                    #
                    # for fits files
                    if iname or ifile:

                        if verbose:
                            print("  . corresponding to pixel " + str(xpix+1) + "," + str(ypix+1))

                        #
                        # writting the intensity of a given pixel for all the channels
                        # into the array flux
                        if ndim == 4:
                            flux = data[0, :, ypix, xpix]

                        if ndim == 3:
                            flux = data[:, ypix, xpix]

                        #
                        # writting the frequencies (in GHz) in the array freqs
                        chans = []
                        freqs = []
                        for channel in range(nchan):
                            chans.append(channel)
                            freqs.append((header.get('CRVAL3') + (channel - header.get('CRPIX3') - 1) * header.get('CDELT3')) / 1.e9)
                        freq = np.array(freqs)


                    #
                    # for ascii files
                    if ispec:

                        #
                        # writting the intensity of a given pixel for all the channels
                        # into the array flux
                        lflux = []
                        freqs = []
                        for line in specfile:
                            freqs.append(float(line.split()[0]))
                            lflux.append(float(line.split()[1]))
                        flux = np.array(lflux)
                        freq = np.array(freqs)
                        nchan = len(flux)
                        flux = flux

                    if cmax or cmean or cmedian or cpercent or cGaussian:
                        all_bins, all_hist, sel_bins, sel_hist, sel_flux = cont_histo(flux, rms_noise)

                    # CONTINUUM FLUX as the maximum of the histogram
                    if cmax:
                        maximum_flux = c_max(flux, rms_noise)

                        continuum_flux_maximum[ypix].append(maximum_flux)

                        if verbose:
                            print("    flux of maximum      = " + str(int(maximum_flux*1.e5)/1.e5))

                    # CONTINUUM FLUX as the mean of the histogram
                    if cmean:
                        mean_flux = np.mean(flux)
                        mean_variance = np.var(flux)
                        mean_sigma = np.std(flux)

                        continuum_flux_mean[ypix].append(mean_flux)

                        if verbose:
                            print("    flux of mean (all)   = " + str(int(mean_flux*1.e5)/1.e5))

                    # CONTINUUM FLUX as the median of the histogram
                    if cmedian:
                        median_flux = np.median(flux)
                        median_variance = np.var(flux)
                        median_sigma = np.std(flux)

                        continuum_flux_median[ypix].append(median_flux)

                        if verbose:
                            print("    flux of median (all) = " + str(int(median_flux*1.e5)/1.e5))

                    # CONTINUUM FLUX as the percentile(s) of the histogram
                    if cpercent:
                        percent25_flux = c_percent(flux, percentile=25)
                        percent75_flux = c_percent(flux, percentile=75)

                        continuum_flux_percent25[ypix].append(percent25_flux)
                        continuum_flux_percent75[ypix].append(percent75_flux)

                        if verbose:
                            print("    flux of percent 25   = " + str(int(percent25_flux*1.e5)/1.e5))
                            print("    flux of percent 75   = " + str(int(percent75_flux*1.e5)/1.e5))


                    # CONTINUUM FLUX as the mean of the "Gaussian" distribution
                    if cmean or cGaussian:
                        meansel_flux = np.mean(sel_flux)
                        meansel_variance = np.var(sel_flux)
                        meansel_sigma = np.std(sel_flux)

                        if cmean:
                            continuum_flux_meansel[ypix].append(meansel_flux)

                            if verbose:
                                print("    flux of mean (sel)   = " + str(int(meansel_flux*1.e5)/1.e5))

                    # CONTINUUM FLUX as the median of the "Gaussian" distribution
                    if cmedian or cGaussian:
                        mediansel_flux = np.median(sel_flux)
                        mediansel_variance = np.var(sel_flux)
                        mediansel_sigma = np.std(sel_flux)

                        if cmedian:
                            continuum_flux_mediansel[ypix].append(mediansel_flux)

                            if verbose:
                                print("    flux of median (sel) = " + str(int(mediansel_flux*1.e5)/1.e5))

                    # CONTINUUM FLUX as the maximum of a KDE distribution
                    if cKDEmax:
                        KDE_bandwidth = rms_noise/10.
                        scipy_kde = stats.gaussian_kde(flux, bw_method=KDE_bandwidth)
                        KDExmin, KDExmax = min(flux), max(flux)
                        KDEx = np.mgrid[KDExmin:KDExmax:100j]
                        positions = np.vstack([KDEx.ravel()])
                        KDEpos = scipy_kde(positions)
                        KDEmax_flux = positions.T[np.argmax(KDEpos)]

                        continuum_flux_KDEmax[ypix].append(KDEmax_flux)

                        if verbose:
                            print("    flux of KDEmax       = " + str(int(KDEmax_flux*1.e5)/1.e5))

                    # CONTINUUM FLUX from a global Gaussian fit
                    if cGaussian:
                        init = [all_hist.max(), meansel_flux, meansel_sigma]
                        out = leastsq(errfunc, init, args=(all_bins, all_hist))
                        c = out[0]
                        Gaussian_flux = c[1]
                        Gaussian_noise = c[2]

                        continuum_flux_Gaussian[ypix].append(Gaussian_flux)
                        continuum_noise_Gaussian[ypix].append(Gaussian_noise)

                        if verbose:
                            print("    flux of Gaussian     = " + str(int(Gaussian_flux*1.e5)/1.e5))
                            print("                         = " + str(int(Gaussian_noise*1.e5)/1.e5))

                    # CONTINUUM FLUX as the peak of the "Gaussian" fit
                    if cGaussian:
                        init = [all_hist.max(), meansel_flux, meansel_sigma]
                        out = leastsq(errfunc, init, args=(sel_bins, sel_hist))
                        d = out[0]
                        GaussNw_flux = d[1]
                        GaussNw_noise = d[2]

                        continuum_flux_GaussNw[ypix].append(GaussNw_flux)
                        continuum_noise_GaussNw[ypix].append(GaussNw_noise)

                        if verbose:
                            print("    flux of Gauss (sel)  = " + str(int(GaussNw_flux*1.e5)/1.e5))
                            print("                         = " + str(int(GaussNw_noise*1.e5)/1.e5))

                    # CONTINUUM FLUX from the method of sigma-clipping
                    if csigmaclip:

                        sigmaclip_flux, sigmaclip_noise = c_sigma_clip(flux, rms_noise)

                        # determination of the continuum level (mean)
                        # and noise/uncertainty (dispersion/standard deviation)

                        continuum_flux_sigmaclip[ypix].append(sigmaclip_flux)
                        continuum_noise_sigmaclip[ypix].append(sigmaclip_noise)

                        if verbose:
                            print("    flux of sigma-clip   = " + str(int(sigmaclip_flux*1.e5)/1.e5))
                            print("    error sigma-clip     = " + str(int(sigmaclip_noise*1.e5)/1.e5))

                    #
                    # creating some plots with spectra and continuum levels
                    if plots:
                        fig_file = plots_path + tmp_file + '_' + str(xpix+1) + '_' + str(ypix+1) + '.png'
                        fig1 = plt.figure()
                        gs = gridspec.GridSpec(3,1)
                        ax1 = fig1.add_subplot(gs[0,0])
                        ax1.axis('off')
                        ax2 = fig1.add_subplot(gs[1:3,0])

                        ax2.plot(freq,flux, 'k-')
                        plt.xlim(freq.min(), freq.max())
                        plt.ylim(flux.min(), flux.max())
                        ax1.set_title('Spectrum and continuum level at pixel (' + str(xpix+1) + ',' + str(ypix+1) + ')')
                        if cmax:
                            ax2.axhline(y=maximum_flux, linestyle='--', color='green', linewidth='1.5')
                            if iname or ifile:
                                ax1.text(0.0, 0.9, "Maximum = " + str(int(maximum_flux*1.e5)/1.e5) + " " + bunit)
                            if ispec:
                                ax1.text(0.0, 0.8, "Maximum = " + str(int(maximum_flux*1.e5)/1.e5))
                        if cmean:
                            ax2.axhline(y=mean_flux, linestyle='--', color='orange', linewidth='1.5')
                            ax2.axhline(y=meansel_flux, linestyle='--', color='yellow', linewidth='1.5')
                            if iname or ifile:
                                ax1.text(0.0, 0.6, "Mean = " + str(int(mean_flux*1.e5)/1.e5) + " " + bunit)
                                ax1.text(0.0, 0.4, "Mean (sel.) = " + str(int(meansel_flux*1.e5)/1.e5) + " " + bunit)
                            if ispec:
                                ax1.text(0.0, 0.6, "Mean = " + str(int(mean_flux*1.e5)/1.e5))
                                ax1.text(0.0, 0.4, "Mean (sel.) = " + str(int(meansel_flux*1.e5)/1.e5))
                        if cmedian:
                            ax2.axhline(y=median_flux, linestyle='--', color='orange', linewidth='1.5')
                            ax2.axhline(y=mediansel_flux, linestyle='--', color='yellow', linewidth='1.5')
                            if iname or ifile:
                                ax1.text(0.0, 0.2, "Median = " + str(int(median_flux*1.e5)/1.e5) + " " + bunit)
                                ax1.text(0.0, 0.0, "Median (sel.) = " + str(int(mediansel_flux*1.e5)/1.e5) + " " + bunit)
                            if ispec:
                                ax1.text(0.0, 0.2, "Median = " + str(int(median_flux*1.e5)/1.e5))
                                ax1.text(0.0, 0.0, "Median (sel.) = " + str(int(mediansel_flux*1.e5)/1.e5))
                        if cpercent:
                            ax2.axhline(y=percent25_flux, linestyle='--', color='red', linewidth='1.5')
                            ax2.axhline(y=percent50_flux, linestyle='--', color='red', linewidth='1.5')
                            ax2.axhline(y=percent75_flux, linestyle='--', color='red', linewidth='1.5')
                            if iname or ifile:
                                ax1.text(0.4, 0.8, "Percent 25th = " + str(int(percent25_flux*1.e5)/1.e5) + " " + bunit)
                            if ispec:
                                ax1.text(0.4, 0.8, "Percent 25th = " + str(int(percent25_flux*1.e5)/1.e5))
                        if cKDEmax:
                            ax2.axhline(y=KDEmax_flux, linestyle='-', color='black', linewidth='1.5')
                            if iname or ifile:
                                ax1.text(0.4, 0.6, "KDE max = " + str(int(KDEmax_flux*1.e5)/1.e5) + " " + bunit)
                            if ispec:
                                ax1.text(0.4, 0.6, "KDE max = " + str(int(KDEmax_flux*1.e5)/1.e5))
                        if cGaussian:
                            ax2.axhline(y=Gaussian_flux, linestyle='-', color='blue', linewidth='3.0', alpha=0.5)
                            ax2.axhline(y=GaussNw_flux, linestyle='-', color='cyan', linewidth='3.0', alpha=0.5)
                            if iname or ifile:
                                ax1.text(0.4, 0.4, "Gaussian = " + str(int(Gaussian_flux*1.e5)/1.e5) + " " + bunit + " (+/- " + str(int(Gaussian_noise*1.e5)/1.e5) + ")")
                                ax1.text(0.4, 0.2, "Gaussian (sel.) = " + str(int(GaussNw_flux*1.e5)/1.e5) + " " + bunit + " (+/- " + str(int(GaussNw_noise*1.e5)/1.e5) + ")")
                            if ispec:
                                ax1.text(0.4, 0.4, "Gaussian = " + str(int(Gaussian_flux*1.e5)/1.e5) + " (+/- " + str(int(Gaussian_noise*1.e5)/1.e5) + ")")
                                ax1.text(0.4, 0.2, "Gaussian (sel.) = " + str(int(GaussNw_flux*1.e5)/1.e5) + " (+/- " + str(int(GaussNw_noise*1.e5)/1.e5) + ")")
                        if csigmaclip:
                            ax2.axhline(y=sigmaclip_flux, linestyle='-', color='red', linewidth='1.5')
                            if iname or ifile:
                                ax1.text(0.4, 0.0, "corrSigma-clip = " + str(int(sigmaclip_flux*1.e5)/1.e5) + " " + bunit + " (+/- " + str(int(sigmaclip_noise*1.e5)/1.e5) + ")")
                            if ispec:
                                ax1.text(0.4, 0.0, "corrSigma-clip = " + str(int(sigmaclip_flux*1.e5)/1.e5) + " (+/- " + str(int(sigmaclip_noise*1.e5)/1.e5) + ")")
                        plt.xlabel('Frequency (GHz)')
                        plt.ylabel('Intensity')

                        fig1.savefig(fig_file)
                        plt.close(fig1)

            #
            # writing the output continuum fits file
            if cfree is False:
                print(" ")
                print("... CONTINUUM FILEs CREATED: ")

            # CONTINUUM FLUX for the sigma-clipping method
            #
            # for the continuum level
            continuum_file_sigmaclip = cont_path + tmp_file + '_continuum' + extension
            os.system('rm -rf ' + continuum_file_sigmaclip)
            if cfree is False:
                print("  . " + continuum_file_sigmaclip)
            if iname or ifile:
                fits.writeto(continuum_file_sigmaclip, np.float32(continuum_flux_sigmaclip), header=header, clobber=True)
            if ispec:
                np.savetxt(continuum_file_sigmaclip, continuum_flux_sigmaclip)
            #
            # for the error/noise on the continuum level
            noise_file_sigmaclip = cont_path + tmp_file + '_noise' + extension
            os.system('rm -rf ' + noise_file_sigmaclip)
            if cfree is False:
                print("  . " + noise_file_sigmaclip)
            if iname or ifile:
                fits.writeto(noise_file_sigmaclip, np.float32(continuum_noise_sigmaclip), header=header, clobber=True)
            if ispec:
                np.savetxt(noise_file_sigmaclip, continuum_noise_sigmaclip)

            # CONTINUUM FLUX for the maximum method
            if cmax:
                #
                continuum_file_maximum = cont_path + tmp_file + '_continuum_maximum' + extension
                os.system('rm -rf ' + continuum_file_maximum)
                if cfree is False:
                    print("  . " + continuum_file_maximum)
                if iname or ifile:
                    fits.writeto(continuum_file_maximum, np.float32(continuum_flux_maximum), header=header, clobber=True)
                if ispec:
                    np.savetxt(continuum_file_maximum, continuum_flux_maximum)

            # CONTINUUM FLUX for the mean method
            if cmean:
                #
                continuum_file_mean = cont_path + tmp_file + '_continuum_mean' + extension
                os.system('rm -rf ' + continuum_file_mean)
                if cfree is False:
                    print("  . " + continuum_file_mean)
                if iname or ifile:
                    fits.writeto(continuum_file_mean, np.float32(continuum_flux_mean), header=header, clobber=True)
                if ispec:
                    np.savetxt(continuum_file_mean, continuum_flux_mean)
                #
                continuum_file_meansel = cont_path + tmp_file + '_continuum_meansel' + extension
                os.system('rm -rf ' + continuum_file_meansel)
                if cfree is False:
                    print("  . " + continuum_file_meansel)
                if iname or ifile:
                    fits.writeto(continuum_file_meansel, np.float32(continuum_flux_meansel), header=header, clobber=True)
                if ispec:
                    np.savetxt(continuum_file_meansel, continuum_flux_meansel)

            # CONTINUUM FLUX for the median method
            if cmedian:
                #
                continuum_file_median = cont_path + tmp_file + '_continuum_median' + extension
                os.system('rm -rf ' + continuum_file_median)
                if cfree is False:
                    print("  . " + continuum_file_median)
                if iname or ifile:
                    fits.writeto(continuum_file_median, np.float32(continuum_flux_median), header=header, clobber=True)
                if ispec:
                    np.savetxt(continuum_file_median, continuum_flux_median)
                #
                continuum_file_mediansel = cont_path + tmp_file + '_continuum_mediansel' + extension
                os.system('rm -rf ' + continuum_file_mediansel)
                if cfree is False:
                    print("  . " + continuum_file_mediansel)
                if iname or ifile:
                    fits.writeto(continuum_file_mediansel, np.float32(continuum_flux_mediansel), header=header, clobber=True)
                if ispec:
                    np.savetxt(continuum_file_mediansel, continuum_flux_mediansel)

            # CONTINUUM FLUX for the percentile(s) method
            if cpercent:
                #
                continuum_file_percent25 = cont_path + tmp_file + '_continuum_percent25' + extension
                os.system('rm -rf ' + continuum_file_percent25)
                if cfree is False:
                    print("  . " + continuum_file_percent25)
                if iname or ifile:
                    fits.writeto(continuum_file_percent25, np.float32(continuum_flux_percent25), header=header, clobber=True)
                if ispec:
                    np.savetxt(continuum_file_percent25, continuum_flux_percent25)
                #
                continuum_file_percent50 = cont_path + tmp_file + '_continuum_percent50' + extension
                os.system('rm -rf ' + continuum_file_percent50)
                if cfree is False:
                    print("  . " + continuum_file_percent50)
                if iname or ifile:
                    fits.writeto(continuum_file_percent50, np.float32(continuum_flux_percent50), header=header, clobber=True)
                if ispec:
                    np.savetxt(continuum_file_percent50, continuum_flux_percent50)
                #
                continuum_file_percent75 = cont_path + tmp_file + '_continuum_percent75' + extension
                os.system('rm -rf ' + continuum_file_percent75)
                if cfree is False:
                    print("  . " + continuum_file_percent75)
                if iname or ifile:
                    fits.writeto(continuum_file_percent75, np.float32(continuum_flux_percent75), header=header, clobber=True)
                if ispec:
                    np.savetxt(continuum_file_percent75, continuum_flux_percent75)

            # CONTINUUM FLUX for the KDEmax method
            if cKDEmax:
                #
                continuum_file_KDEmax = cont_path + tmp_file + '_continuum_KDEmax' + extension
                os.system('rm -rf ' + continuum_file_KDEmax)
                if cfree is False:
                    print("  . " + continuum_file_KDEmax)
                if iname or ifile:
                    fits.writeto(continuum_file_KDEmax, np.float32(continuum_flux_KDEmax), header=header, clobber=True)
                if ispec:
                    np.savetxt(continuum_file_KDEmax, continuum_flux_KDEmax)

            # CONTINUUM FLUX for the Gaussian method
            if cGaussian:
                #
                continuum_file_Gaussian = cont_path + tmp_file + '_continuum_Gaussian' + extension
                os.system('rm -rf ' + continuum_file_Gaussian)
                if cfree is False:
                    print("  . " + continuum_file_Gaussian)
                if iname or ifile:
                    fits.writeto(continuum_file_Gaussian, np.float32(continuum_flux_Gaussian), header=header, clobber=True)
                if ispec:
                    np.savetxt(continuum_file_Gaussian, continuum_flux_Gaussian)
                #
                noise_file_Gaussian = cont_path + tmp_file + '_noise_Gaussian' + extension
                os.system('rm -rf ' + noise_file_Gaussian)
                if cfree is False:
                    print("  . " + noise_file_Gaussian)
                if iname or ifile:
                    fits.writeto(noise_file_Gaussian, np.float32(continuum_noise_Gaussian), header=header, clobber=True)
                if ispec:
                    np.savetxt(noise_file_Gaussian, continuum_noise_Gaussian)
                #
                continuum_file_GaussNw = cont_path + tmp_file + '_continuum_GaussNw' + extension
                os.system('rm -rf ' + continuum_file_GaussNw)
                if cfree is False:
                    print("  . " + continuum_file_GaussNw)
                if iname or ifile:
                    fits.writeto(continuum_file_GaussNw, np.float32(continuum_flux_GaussNw), header=header, clobber=True)
                if ispec:
                    np.savetxt(continuum_file_GaussNw, continuum_flux_GaussNw)
                #
                noise_file_GaussNw = cont_path + tmp_file + '_noise_GaussNw' + extension
                os.system('rm -rf ' + noise_file_GaussNw)
                if cfree is False:
                    print("  . " + noise_file_GaussNw)
                if iname or ifile:
                    fits.writeto(noise_file_GaussNw, np.float32(continuum_noise_GaussNw), header=header, clobber=True)
                if ispec:
                    np.savetxt(noise_file_GaussNw, continuum_noise_GaussNw)



    if imerge:

        #
        # re-setting the varibles to the unmerged files
        tmp_path = unmerged_path
        tmp_files = unmerged_files

    for tmp_file in tmp_files:

        if imerge:

            #
            # copy of the merged continuuum to single continuum files
            # for each one of the files used in the merging
            os.system('cp -rp ' + continuum_file_sigmaclip + ' ' + cont_path + tmp_file + '_continuum' + extension)

        #
        # subtraction of the continuum to the data file
        # it produces a line-only cube and a continuum-only file
        if cfree is True:

            print("")
            print("... REMOVING CONTINUUM FROM DATA ... " + tmp_path + tmp_file + extension)

            #
            # defining the cube file (original file)
            # and the continuum file
            cube_file = tmp_path + tmp_file + extension
            cont_files = []
            cont_files.append(tmp_file + '_continuum')

            print(" ")
            print("... FILEs CREATED: ")

            for cont_file in cont_files:

                #
                # for ascii files
                if ispec:

                    fdata_cont = open(cont_path + cont_file + extension, 'r')

                    for line in fdata_cont:

                        data_cont = float(line.strip())

                    line_outfile = line_path + cont_file + '.line' + extension
                    ascii.write((freqs, flux[:]-data_cont), output=line_outfile)

                    print("  . " + line_outfile)

                #
                # for fits files
                if iname or ifile:

                    data_cube = fits.getdata(cube_file)
                    header_cube = fits.getheader(cube_file)
                    data_cont = fits.getdata(cont_path + cont_file + extension)
                    header_cont = fits.getheader(cont_path + cont_file + extension)

                    #
                    # to determine a general noise to be considered
                    nxpix = header_cube.get('NAXIS1')
                    rmspix = int(nxpix / 8)
                    npixmin = 30

                    #
                    # to calculate the rms noise level in four different regions
                    # throughout the continuum fits file
                    if nxpix > npixmin:
                        rms = []
                        rms.append(np.mean(data_cont[rmspix*1:rmspix*2,rmspix*1:rmspix*2]))
                        rms.append(np.mean(data_cont[rmspix*1:rmspix*2,rmspix*2:rmspix*3]))
                        rms.append(np.mean(data_cont[rmspix*1:rmspix*2,rmspix*5:rmspix*6]))
                        rms.append(np.mean(data_cont[rmspix*1:rmspix*2,rmspix*6:rmspix*7]))
                        #
                        rms.append(np.mean(data_cont[rmspix*2:rmspix*3,rmspix*1:rmspix*2]))
                        rms.append(np.mean(data_cont[rmspix*2:rmspix*3,rmspix*2:rmspix*3]))
                        rms.append(np.mean(data_cont[rmspix*2:rmspix*3,rmspix*5:rmspix*6]))
                        rms.append(np.mean(data_cont[rmspix*2:rmspix*3,rmspix*6:rmspix*7]))
                        #
                        rms.append(np.mean(data_cont[rmspix*5:rmspix*6,rmspix*1:rmspix*2]))
                        rms.append(np.mean(data_cont[rmspix*5:rmspix*6,rmspix*2:rmspix*3]))
                        rms.append(np.mean(data_cont[rmspix*5:rmspix*6,rmspix*5:rmspix*6]))
                        rms.append(np.mean(data_cont[rmspix*5:rmspix*6,rmspix*6:rmspix*7]))
                        #
                        rms.append(np.mean(data_cont[rmspix*6:rmspix*7,rmspix*1:rmspix*2]))
                        rms.append(np.mean(data_cont[rmspix*6:rmspix*7,rmspix*2:rmspix*3]))
                        rms.append(np.mean(data_cont[rmspix*6:rmspix*7,rmspix*5:rmspix*6]))
                        rms.append(np.mean(data_cont[rmspix*6:rmspix*7,rmspix*6:rmspix*7]))

                        #print(np.median(rms))
                        data_finalcont = data_cont + np.absolute(np.median(rms))
                        data_line = data_cube - (data_cont + np.absolute(np.median(rms)))

                    #
                    # if the size of the map is too small (less than 100x100 pixels)
                    # no rms noise level is subtracted
                    else:

                        print("  . WARNING: The image has less than %i pixels" % (npixmin))
                        print("  .          No residual noise level subtracted for")
                        print("  .          %s " % (cube_file))

                        data_finalcont = data_cont
                        data_line = data_cube - data_cont

                    cont_outfile = cont_path + cont_file + '.cont' + extension
                    line_outfile = line_path + cont_file + '.line' + extension

                    os.system('rm -rf ' + cont_outfile)
                    fits.writeto(cont_outfile, np.float32(data_finalcont), header=header_cont, clobber=True)

                    os.system('rm -rf ' + line_outfile)
                    fits.writeto(line_outfile, np.float32(data_line), header=header_cube, clobber=True)

                    #
                    # creating a line-only cube
                    # replacing the old continuum-only image by the new one
                    os.system('mv ' + cont_path + cont_file + extension + ' ' + cont_path + cont_file + '_original' + extension)
                    os.system('mv ' + cont_outfile + ' ' + cont_path + cont_file + extension)
                    os.system('mv ' + line_outfile + ' ' + cont_path + tmp_file + '_line' + extension)

                    print("  . " + cont_path + tmp_file + '_continuum' + extension)
                    print("  . " + cont_path + tmp_file + '_noise' + extension)
                    print("  . " + cont_path + tmp_file + '_line' + extension)

    #
    # combintation of several continuum files to determine the spectral index
    # it also combines them and creates a continuum model, and line+continuum model
    if spindex:

        print(" ")
        print("+++ DETERMINING SPECTRAL INDEX (CONTINUUM MODEL) ...")

        #
        # to combine all the continuum images in one single cube
        my_frequency = []
        my_cube = np.empty([len(tmp_files),nypix,nxpix])
        icount = 0

        for tmp_file in tmp_files:

            contcube_data = fits.getdata(cont_path + tmp_file + '_continuum' + extension)
            contcube_header = fits.getheader(cont_path + tmp_file + '_continuum' + extension)
            my_frequency.append(contcube_header.get('RESTFRQ'))
            my_cube[icount,:,:] = contcube_data
            icount = icount+1

        #
        # applying a blanking value of 1.e-10
        my_cube[my_cube<1.e-10] = 1.e-10
        my_frequency = np.array(my_frequency)

        #
        # to fit all the continuum images and determine the spectral index (m)
        # following expression:  flux = A x frequency ^ spindex
        # and in logarithmic version:  log(flux) = log(A) + spindex * log(frequency)
        m = []
        n = []

        for ypix in range(nypix):

            print("... analyzing column " + str(ypix+1) + " out of " + str(nypix))

            m.append([])
            n.append([])

            for xpix in range(nxpix):

                y = np.log10(my_cube[:,ypix,xpix])
                x = np.log10(my_frequency)
                z = np.polyfit(x, y, 1)
                m[ypix].append(z[0])
                n[ypix].append(z[1])

        cont_m_file = cont_path + source + '_spindex.fits'
        cont_n_file = cont_path + source + '_intercept.fits'

        os.system('rm -rf ' + cont_m_file)
        fits.writeto(cont_m_file, np.float32(m), header=header, clobber=True)

        os.system('rm -rf ' + cont_n_file)
        fits.writeto(cont_n_file, np.float32(n), header=header, clobber=True)

        #
        # reading the spectral index (m) and the intercept (n)
        m_data = fits.getdata(cont_m_file)
        n_data = fits.getdata(cont_n_file)

        #
        # to create a continuum model (from the spectral index)
        # and a line+continuum cube, using the continuum model
        for tmp_file in tmp_files:

            #
            # to create a frequency array from the original data cube
            cube_data =  fits.getdata(tmp_path + tmp_file + extension)
            cube_header = fits.getheader(tmp_path + tmp_file + extension)
            frequency_array = []

            for i in range(cube_header['NAXIS3']):

                frequency_array.append(cube_header['CRVAL3'] + (i - cube_header['CRPIX3'])*cube_header['CDELT3'])

            frequency_array = np.array(frequency_array)

            #
            # to read the continuum real emission (determined in --continuum mode)
            cont_real = fits.getdata(cont_path + tmp_file + '_continuum.fits')

            #
            # to generate the continuum model (using the n, m parameters)
            print("... creating continuum model of " + tmp_file)
            cont_model = np.power(10, n_data[np.newaxis, :, :] + (np.log10(frequency_array[:,np.newaxis, np.newaxis]) * m_data[np.newaxis, :, :]))
            cont_model_file = cont_path + tmp_file + '_cont_model' + extension
            os.system('rm -rf ' + cont_model_file)
            fits.writeto(cont_model_file, np.float32(cont_model), header=cube_header, clobber=True)

            #
            # to compute the factor between real and model continuum
            factor = cont_model[:, :, :] / cont_real[np.newaxis, :, :]
            factor[cont_real < (rms_noise*1.0e-3)] = 1.0
            factor_file = cont_path + tmp_file + '_cont_factor' + extension
            os.system('rm -rf ' + factor_file)
            fits.writeto(factor_file, np.float32(factor), header=cube_header, clobber=True)

            #
            # to compute the line and continuum model
            print("... creating the line+continuum model of " + tmp_file)
            line_cont_model = cube_data[:, :, :, :] * factor[np.newaxis, :, :, :]
            line_cont_model_file = cont_path + tmp_file + '_line_cont_model' + extension
            os.system('rm -rf ' + line_cont_model_file)
            fits.writeto(line_cont_model_file, np.float32(line_cont_model), header=cube_header, clobber=True)

            #
            # line only model data
            line_model = line_cont_model - cont_model
            line_model_file = cont_path + tmp_file + '_line_model' + extension
            os.system('rm -rf ' + line_model_file)
            fits.writeto(line_model_file, np.float32(line_model), header=cube_header, clobber=True)

        #
        # to indicate where the created files can be found
        print("... FILEs CREATED are found in " + cont_path)
        print("  . search for spindex, cont_model and line_cont_model")
        print(" ")

def cont_histo(flux, rms_noise):
    """
    Create histogram distribution of the flux data
    and select a narrower range around the maximum
    of the histogram distribution
    
    Parameters
    ----------
    flux : np.ndarray
        One-dimension array of flux values
    rms_noise : float
        The estimated RMS noise level of the data
    
    Returns
    -------
    all_bins : np.ndarray
        One-dimension array with the value of bins of the histogram
    all_hist : np.ndarray
        One-dimension array with the value of the position of the bins
    sel_bins : np.ndarray
        One-dimension array with the value of bins of the histogram
        for the selected bins around the maximum
    sel_hist : np.ndarray
        One-dimension array with the value of position of the bins
        for the selected bins around the maximum
    sel_flux : np.ndarray
        One-dimension array of the flux values selected
        around the maximum of the histogram
    """
    
    #
    # creating a general histogram of the flux data
    # main variables are:
    #   all_hist     - counts in each bin of the histogram
    #   all_bins     - location of the bins (fluxes)
    #   all_number_* - index of the array
    number_bins = int((np.amax(flux)-np.amin(flux))/(2*rms_noise))
    all_hist, all_bin_edges = np.histogram(flux, number_bins)
    all_bins = all_bin_edges[0:len(all_bin_edges)-1]
    all_bins = [x + (all_bins[1]-all_bins[0])/2. for x in all_bins]
    all_number_max_array = (np.where(all_hist == all_hist.max())[0])
    all_number_max = all_number_max_array[0]
    all_bins_max = (all_bin_edges[all_number_max] + (all_bins[1]-all_bins[0])/2.)

    # Gaussian fit around the maximum of the distribution
    # determining the range to fit the Gaussian function
    all_number_left  = (np.where(((all_hist == 0) & (all_bins <= all_bins_max)) | (all_bins == all_bins[0]))[0]).max()
    all_number_right = (np.where(((all_hist == 0) & (all_bins >= all_bins_max)) | (all_bins == all_bins[number_bins-1]))[0]).min()
    all_number_total = abs(all_number_right-all_number_max)+abs(all_number_left-all_number_max)
    emission_absorption_ratio = abs(all_number_right-all_number_max)*1.0/(all_number_total*1.0)
    if (emission_absorption_ratio >= 0.66):
        lower_all_bins = all_bins_max - 8. * (all_bins[1]-all_bins[0])
        upper_all_bins = all_bins_max + 4. * (all_bins[1]-all_bins[0])
    if (emission_absorption_ratio <= 0.33):
        lower_all_bins = all_bins_max - 4. * (all_bins[1]-all_bins[0])
        upper_all_bins = all_bins_max + 8. * (all_bins[1]-all_bins[0])
    if ((emission_absorption_ratio > 0.33) and (emission_absorption_ratio < 0.66)):
        lower_all_bins = all_bins_max - 5. * (all_bins[1]-all_bins[0])
        upper_all_bins = all_bins_max + 5. * (all_bins[1]-all_bins[0])
    sel_bins_array = np.where((all_bins >= lower_all_bins) & (all_bins <= upper_all_bins))[0]
    if (len(sel_bins_array) < 3):
        sel_bins_array = [sel_bins_array[0]-2, sel_bins_array[0]-1, sel_bins_array[0], sel_bins_array[0]+1, sel_bins_array[0]+2]
        lower_all_bins = all_bins[sel_bins_array[0]]
        upper_all_bins = all_bins[sel_bins_array[len(sel_bins_array)-1]]
    sel_bins = all_bins[sel_bins_array[0]:sel_bins_array[len(sel_bins_array)-1]+1]
    sel_hist = all_hist[sel_bins_array[0]:sel_bins_array[len(sel_bins_array)-1]+1]
    sel_flux = flux[(flux >= lower_all_bins) & (flux <= upper_all_bins)]
    
    return all_bins, all_hist, sel_bins, sel_hist, sel_flux


def c_max(flux, rms_noise):
    """
    Perform histogram distribution of variable flux, and determine
    the flux level of the maximum of the histogram
    
    Parameters
    ----------
    flux : np.ndarray
        One-dimension array of flux values
    rms_noise : float
        The estimated RMS noise level of the data
    
    Returns
    -------
    maximum_flux : float
        The measured continuum flux as the maximum of the histogram
    """
    all_bins, all_hist, sel_bins, sel_hist, sel_flux = cont_histo(flux, rms_noise)
    
    maximum_flux = all_bins[(np.where(all_hist == all_hist.max())[0])[0]]
    
    return maximum_flux

def c_percent(flux, percentile):
    """
    Perform numpy percentile to determine the level of the selected
    percentile
    
    Parameters
    ----------
    flux : np.ndarray
        One-dimension array of flux values
    percentile : float
        The selected percentile
    
    Returns
    -------
    percent_flux : float
        The measured continuum flux at the selected percentile
    """
    
    percent_flux = np.percentile(flux, percentile)
    
    return percent_flux

def c_sigma_clip(flux, rms_noise, sigma_clip_threshold=1.8):
    """
    Perform sigma-clipping to determine the mean flux level, with different
    adaptations for emission- and absorption-dominated spectra

    Parameters
    ----------
    flux : np.ndarray
        One-dimension array of flux values
    rms_noise : float
        The estimated RMS noise level of the data
    sigma_clip_threshold : float
        The threshold in number of sigma above/below which to reject outlier
        data

    Returns
    -------
    sigmaclip_flux : float
    sigmaclip_noise : float
        The measured continuum flux and estimated 1-sigma per-channel noise
        around that measurement
    """

    # sigma-clipping method applied to the flux array
    if astropy.version.major >= 1:
        filtered_data = sigma_clip(flux, sigma=sigma_clip_threshold,
                                   iters=None)
    elif astropy.version.major < 1:
        filtered_data = sigma_clip(flux, sig=sigma_clip_threshold, iters=None)

    # for deep absorptions/emission like SgrB2
    sigmaclip_flux_prev = sigmaclip_flux = np.mean(filtered_data)
    sigmaclip_noise = sigmaclip_sigma = np.std(filtered_data)
    mean_flux = np.mean(flux)

    # for EMISSION-dominated spectra
    if (mean_flux-sigmaclip_flux_prev) > (+1.0*rms_noise):
        sigmaclip_flux = sigmaclip_flux_prev - sigmaclip_sigma
    # for ABSORPTION-dominated spectra
    elif (mean_flux-sigmaclip_flux_prev) < (-1.0*rms_noise):
        sigmaclip_flux = sigmaclip_flux_prev + sigmaclip_sigma

    return sigmaclip_flux, sigmaclip_noise
