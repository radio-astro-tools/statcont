import os
import numpy
import warnings
import astropy.io

##=====================================================================================================================
def cutout(filename, xc, yc, size, outfile):

    """
    Inputs:
        filename    - .fits filename
        xc,yc       - x and y center coordinates in the fits files' coordinate system (CTYPE)
        size        - size of the cutout in pixel
        outfile     - output file
    """

    if isinstance(filename, str):
        fitsfile = astropy.io.fits.open(filename)
        opened = True

    #get header information
    head = fitsfile[0].header.copy()
    ndim = head.get('NAXIS')
    cdelt1 = head.get('CDELT1')
    cdelt2 = head.get('CDELT2')
    nmax = head.get('NAXIS3')

    if cdelt1 is None or cdelt2 is None:
        raise Exception("Missing CD or CDELT keywords in header")

    xmin, xmax = numpy.max([0, xc - size / 2]), numpy.min([head['NAXIS1'], xc + size / 2])
    ymin, ymax = numpy.max([0, yc - size / 2]), numpy.min([head['NAXIS2'], yc + size / 2])

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

    # copy data

    if ndim == 4:
        img = fitsfile[0].data[0:1, 0:nmax, ymin:ymax, xmin:xmax]
    if ndim == 3:
        img = fitsfile[0].data[0:nmax, ymin:ymax, xmin:xmax]
    newfile = astropy.io.fits.PrimaryHDU(data=img, header=head)
###ASM    print 'Cut image %s with dims %s to %s.  xrange: %f:%f, yrange: %f:%f' % (filename, fitsfile[0].data.shape, img.shape, xmin, xmax, ymin, ymax)
###ASM    print 'New outputfile is %s' % outfile

    # output
    if isinstance(outfile, str):

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')

            newfile.writeto(outfile, clobber=True)

    if opened:
        fitsfile.close()

###ASM    # export outfile to MIRIAD format; remove any existing file
###ASM    os.system('rm -rf %s' % outfile.replace('.fits', ''))
###ASM    os.system('fits op=xyin in=%s out=%s' % (outfile, outfile.replace('.fits', '')))

    return newfile



