import os
import tools
import numpy as np
from pathlib import Path
from astropy.io import fits







class MasterFlat:
    """
    Class that creates the master flat image. This image is obtained by taking the median of multiple flat images,
    and then subtracting the master bias file and the rescaled master dark file.
    """

    def __init__(self):
        pass










    def run(self, rawFlatImagePaths, masterBiasPath, masterDarkPath, outputFilePath=None):
        """
        Create a master flat field

        Input:
            rawFlatImagePaths: list of strings containing the full path of the 2D raw flatfield CCD FITs files
            masterBiasPath:    string containing the full path of the 2D master bias FITS file
            masterDarkPath:    string containing the full path of the 2D master dark FITS file
            outputFilePath:    If None, nothing is saved. Otherwise, a string with the name of the outputfile,
                               incl. the extension ".fits".

        Output:
            masterFlat: 2D numpy array:  master flatfield image [ADU]
        """

        # Get all the relevant CCD images

        masterBias = tools.getImage(masterBiasPath)
        masterDark = tools.getImage(masterDarkPath)
        flatImages = tools.getImages(rawFlatImagePaths)

        # Take the median of all flat images

        masterFlat = np.median(flatImages, axis=0)

        # Find out the exposure time for the flat images 

        flatFileStem = Path(rawFlatImagePaths[0]).stem
        flatExposureTime = float(flatFileStem[-4:])

        # Find out the exposure time for the master dark image

        darkFileStem = Path(masterDarkPath).stem
        darkExposureTime = float(darkFileStem[-4:])

        # Correct the master flat for the bias and the dark current
        # Note that the masterDark was already corrected for the bias 

        masterFlat = masterFlat - masterBias - masterDark / darkExposureTime * flatExposureTime

        # Zero all negative values 

        masterFlat[masterFlat < 0] = 0

        # If required, save to a FITS file

        if outputFilePath is not None:
            num_row, num_col = masterFlat.shape
            hdr = fits.Header()
            hdr["rows"] = num_row
            hdr["cols"] = num_col
            hdr["std_dark"] = np.std(masterDark)
            outputParentPath = Path(outputFilePath).parent.absolute()
            if not os.path.exists(outputParentPath):
                os.makedirs(outputParentPath)
            hdu = fits.PrimaryHDU(masterFlat, header=hdr)
            hdu.writeto(outputFilePath, overwrite=True)

        # That's it!

        return masterFlat




