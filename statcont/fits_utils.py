import os
import warnings
import astropy.io
from astropy.io import fits
import numpy as np

##======================================================================
def fits_cutout(filename, xc, yc, size, outfile):

    """
    Create a cutout of a larger FITS file
    
    Parameters
    ----------
    filename : string
        Name of the FITS file
    xc : float
    yc : float
        X and Y coordinates of the central pixel of the cutout image
    size : float
        Size of the cutout image in pixels
    outfile : string
        Name of the output FITS file
    
    Returns
    -------
    outfile : string
        Name of the output FITS file
    """

    if isinstance(filename, str):
        fitsfile = astropy.io.fits.open(filename)
        opened = True

    # Get header information
    head = fitsfile[0].header.copy()
    ndim = head.get('NAXIS')
    cdelt1 = head.get('CDELT1')
    cdelt2 = head.get('CDELT2')
    nmax = head.get('NAXIS3')

    if cdelt1 is None or cdelt2 is None:
        raise Exception("Missing CD or CDELT keywords in header")

    xmin, xmax = np.max([0, xc - size / 2]), np.min([head['NAXIS1'], xc + size / 2])
    ymin, ymax = np.max([0, yc - size / 2]), np.min([head['NAXIS2'], yc + size / 2])
    
    xmin = int(xmin)
    xmax = int(xmax)
    ymin = int(ymin)
    ymax = int(ymax)

    if xmax < 0 or ymax < 0:
        raise ValueError("Max Coordinate is outside of map: %f, %f." % (xmax, ymax))

    if ymin >= head.get('NAXIS2') or xmin >= head.get('NAXIS1'):
        raise ValueError("Min Coordinate is outside of map: %f, %f." % (xmin, ymin))

    head['CRPIX1'] -= xmin
    head['CRPIX2'] -= ymin
    head['NAXIS1'] = int(xmax - xmin)
    head['NAXIS2'] = int(ymax - ymin)
    head['NAXIS3'] = 1

    if head.get('NAXIS1') == 0 or head.get('NAXIS2') == 0:
        raise ValueError("Map has a 0 dimension: %i, %i." % (head.get('NAXIS1'), head.get('NAXIS2')))

    # Copy data
    if ndim == 4:
        img = fitsfile[0].data[0:1, 0:nmax, ymin:ymax, xmin:xmax]
    if ndim == 3:
        img = fitsfile[0].data[0:nmax, ymin:ymax, xmin:xmax]
    newfile = astropy.io.fits.PrimaryHDU(data=img, header=head)

    # Output
    if isinstance(outfile, str):

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')

            newfile.writeto(outfile, clobber=True)

    if opened:
        fitsfile.close()

    return newfile

##======================================================================
def fits_merge(tmp_files, tmp_path, extension, merged_file_name, merged_path):
    """
    Merge files to be processed into one single output merged file
    
    Parameters
    ----------
    tmp_files : np.ndarray
    tmp_path : string
        One-dimension array with the names of the files to be merged
        and path where the files are saved
    extension : string
        Extension of the input files, either .fits or .dat
    merged_file_name : string
    merged_path : string
        Strings with the name of the output merged file and the
        path to this file
    
    Returns
    -------
    tmp_files : np.ndarray
        One-dimension array with the name of the merged file
    tmp_path : np.ndarray
        One-dimension array with the path to the merged file
    unmerged_files : np.ndarray
        One-dimension array with the names of the individual files
        that were merged
    unmerged_path : np.ndarray
        One-dimension array with the path to the individual files
    """
    
    merged_file = merged_path + merged_file_name + extension
    print (merged_file)
    
    if extension=='.dat':
        filenames = []
        for tmp_file in tmp_files:
            filenames.append(tmp_path+tmp_file+extension)
        with open(merged_file, 'w') as outfile:
            for fname in filenames:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)
    
    if extension=='.fits':
        # Define ranges for channels, to combine all the files
        ninterval = len(tmp_files)
        bnchan = np.empty([ninterval])
        enchan = np.empty([ninterval])

        # Determine the total number of channels (icount) of the merged file
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

        # Create output merged file
        merged_data = np.empty([npolz,enchan[icount-1],nypix,nxpix])
        
        # Write data to the new merged file
        icount = 0
        for tmp_file in tmp_files:
            data = fits.getdata(tmp_path + tmp_file + extension)
            merged_data[0,bnchan[icount]:enchan[icount],:,:] = data
            icount = icount + 1

        os.system('rm -rf ' + merged_file)
        fits.writeto(merged_file, np.float32(merged_data), header=header, clobber=True)
        
    # Keep track of the individual file names, and the merged file
    unmerged_files = tmp_files
    unmerged_path = tmp_path
    tmp_files = []
    tmp_files.append(merged_file_name)
    tmp_path = merged_path
    
    return tmp_files, tmp_path, unmerged_files, unmerged_path
