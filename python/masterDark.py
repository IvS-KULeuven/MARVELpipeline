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

        biasLevel = np.median(masterBias)
        masterDark = np.median(darkImages, axis=0) - biasLevel

        # If required, save to a FITS file

        if outputFileName is not None:
            num_row, num_col = masterDark.shape
            hdr = fits.Header()
            hdr["rows"] = num_row
            hdr["cols"] = num_col
            outputParentPath = Path(outputFileName).parent.absolute()
            if not os.path.exists(outputParentPath):
                os.makedirs(outputParentPath)
            hdu = fits.PrimaryHDU(masterDark, header=hdr)
            hdu.writeto(outputFileName, overwrite=True)

        # That's it!

        return masterDark









if __name__ == "__main__":

    t1 = time.time()

    params   = yaml.safe_load(open("params.yaml"))

    root     = (params["Configuration"])["rootFolder"]
    master_dark_params = params["MasterDarkImage"]
    master_bias_params = params["MasterBiasImage"]

    raw_dark_paths = [ root+path for path in master_dark_params["inputPath"] ]
    master_bias_path = root + master_bias_params["outputPath"]

    masterDark = MasterDark()
    masterDark.run(raw_dark_paths, master_bias_path, root + master_dark_params["outputPath"])

    t2 = time.time()

    print(f"[{t2-t1}]")


