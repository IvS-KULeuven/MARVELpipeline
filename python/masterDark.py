import os
import numpy as np
import time
import yaml
from pathlib import Path
from astropy.io import fits
import tools




    
class MasterDark:
    """
    Class that creates the master Dark image. This image is obtained by taking the median of multiple dark images,
    and then subtracting the master bias file. Without this electronic offset, the dark current is the only thing
    that remains.
    """

    def __init__(self):
        pass




    def run(self, rawDarkImagePaths, masterBiasPath, outputFileName=None):
        """
        Create the Dark image

        Input:
            rawDarkImagePaths: list of strings containing the full path of the 2D raw dark CCD FITs files
            masterBiasPath:    string containing the full path of the 2D master bias FITS file
            outputFileName:    If None, nothing is saved. Otherwise, a string with the name of the outputfile,
                               incl. the extension ".fits".

        Output:
            masterDark:     master dark image [ADU]
        """

        # Get all the relevant CCD images

        masterBias  = tools.getImage(masterBiasPath)
        darkImages = tools.getImages(rawDarkImagePaths)

        # Compute the master dark as the median of the raw dark images, then subtract the electronic offset

        masterDark = np.median(darkImages, axis=0) - masterBias

        # We assume all darks have the same exposure time. Use the one of the
        # first raw dark to extract it, so that we can store it in the master dark as well.

        hdulist = fits.open(rawDarkImagePaths[0])
        hdu = hdulist[0]
        darkExposureTime = float(hdu.header['EXPTIME'])

        # If required, save to a FITS file

        if outputFileName is not None:
            num_row, num_col = masterDark.shape
            hdr = fits.Header()
            hdr["rows"] = num_row
            hdr["cols"] = num_col
            hdr["std_dark"] = np.std(masterDark)
            hdr["exptime"] = darkExposureTime
            outputParentPath = Path(outputFileName).parent.absolute()
            if not os.path.exists(outputParentPath):
                os.makedirs(outputParentPath)
            hdu = fits.PrimaryHDU(masterDark, header=hdr)
            hdu.writeto(outputFileName, overwrite=True)

        # That's it!

        return masterDark


