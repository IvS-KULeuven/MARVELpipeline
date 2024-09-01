
import os
import yaml
import tools
import numpy as np
import time
from pathlib import Path
from astropy.io import fits







class MasterThAr:
    """
    Class that creates the master ThAr image. This image is obtained by taking the median of multiple ThAr images,
    and then subtracting the master bias file.
    """

    def __init__(self):
        pass




    def run(self, rawThArImagePaths, masterBiasPath, outputFileName=None):
        """
        We run through the algorithm to create the master ThAr images.

        Input:
            rawThArImagePaths: list of strings containing the full path of the 2D raw ThAr CCD FITs files
            masterBiasPath:    string containing the full path of the 2D master bias FITS file
            outputFileName:    If None, nothing is saved. Otherwise, a string with the name of the outputfile,
                               incl. the extension ".fits".

        Output:
            masterFlat:     master flat image [ADU]
        """

        # Get all the fits files corresponding to these hashes

        masterBias  = tools.getImage(masterBiasPath)
        ThArImages = tools.getImages(rawThArImagePaths)

        # Use the image in the fits files, and use mean_combining to obtain the the master image

        masterThAr = np.median(ThArImages, axis=0) - masterBias

        # Add offset so that all the values in the master ThAr image are positive

        if np.min(masterThAr) < 0:
                  masterThAr = masterThAr - np.min(masterThAr)

        # If required, save to a FITS file

        if outputFileName is not None:
            num_row, num_col = masterThAr.shape
            hdr = fits.Header()
            hdr["rows"] = num_row
            hdr["cols"] = num_col
            outputParentPath = Path(outputFileName).parent.absolute()
            if not os.path.exists(outputParentPath):
                os.makedirs(outputParentPath)
            hdu = fits.PrimaryHDU(masterThAr, header=hdr)
            hdu.writeto(outputFileName, overwrite=True)

        # That's it!

        return masterThAr









if __name__ == "__main__":

    t1 = time.time()

    params   = yaml.safe_load(open("params.yaml"))

    root     = (params["Configuration"])["rootFolder"]
    masterThAr_params = params["MasterThArImage"]
    masterbias_params = params["MasterBiasImage"]

    # Mater Flat Image
    raw_ThAr_paths = [ root+path for path in masterThAr_params["inputPath"] ]
    master_bias_path = root + masterbias_params["outputPath"]

    masterThAr = MasterThAr()
    masterThAr.run(raw_ThAr_paths, master_bias_path, root + masterThAr_params["outputPath"])

    t2 = time.time()

    print(f"[{t2-t1}]")

