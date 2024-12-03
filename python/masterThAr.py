
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




    def run(self, rawThArImagePaths, masterBiasPath, masterDarkPath, outputFilePath=None):
        """
        Create the ThAr master image

        Input:
            rawThArImagePaths: list of strings containing the full path of the 2D raw ThAr CCD FITs files
            masterBiasPath:    string containing the full path of the 2D master bias FITS file
            masterDarkPath:    string containing the full path of the 2D master dark FITS file
            outputFilePath:    If None, nothing is saved. Otherwise, a string with the name of the outputfile,
                               incl. the extension ".fits".

        Output:
            masterThAr: 2D numpy array:  master ThAr image [ADU]
        """

        # Get all the relevant CCD images

        masterBias = tools.getImage(masterBiasPath)
        masterDark = tools.getImage(masterDarkPath)
        ThArImages = tools.getImages(rawThArImagePaths)

        # Take the median of all ThAr images

        masterThAr = np.median(ThArImages, axis=0)

        # Find out the exposure time for the ThAr images 

        ThArFileStem = Path(rawThArImagePaths[0]).stem
        ThArExposureTime = float(ThArFileStem[-4:])

        # Find out the exposure time for the master dark image

        darkFileStem = Path(masterDarkPath).stem
        darkExposureTime = float(darkFileStem[-4:])

        # Correct the master ThAr for the bias and the dark current

        masterThAr = masterThAr - masterBias - masterDark / darkExposureTime * ThArExposureTime

        # Zero all negative values 

        masterThAr[masterThAr < 0] = 0
        
        # If required, save to a FITS file

        if outputFilePath is not None:
            num_row, num_col = masterThAr.shape
            hdr = fits.Header()
            hdr["rows"] = num_row
            hdr["cols"] = num_col
            outputParentPath = Path(outputFilePath).parent.absolute()
            if not os.path.exists(outputParentPath):
                os.makedirs(outputParentPath)
            hdu = fits.PrimaryHDU(masterThAr, header=hdr)
            hdu.writeto(outputFilePath, overwrite=True)

        # That's it!

        return masterThAr









if __name__ == "__main__":

    t1 = time.time()

    params   = yaml.safe_load(open("params.yaml"))

    rootFolderRawData  = params["Configuration"]["rootFolderRawData"]
    rootFolderProcessedData = params["Configuration"]["rootFolderProcessedData"]
    masterThAr_params = params["MasterThArImage"]
    masterbias_params = params["MasterBiasImage"]
    masterdark_params = params["MasterDarkImage"]

    raw_ThAr_paths = [ rootFolderRawData + path for path in masterThAr_params["inputPath"] ]
    master_bias_path = rootFolderProcessedData + masterbias_params["outputPath"]
    master_dark_path = rootFolderProcessedData + masterdark_params["outputPath"]
    ThAr_output_path = rootFolderProcessedData + masterThAr_params["outputPath"]

    masterThAr = MasterThAr()
    masterThAr.run(raw_ThAr_paths, master_bias_path, master_dark_path, ThAr_output_path)

    t2 = time.time()

    print(f"[{t2-t1}]")

