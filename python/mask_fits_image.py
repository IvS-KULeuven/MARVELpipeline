

import os
import yaml
import numpy as np
import time
from pathlib import Path
from astropy.io import fits







class MaskImageCreator:
    """
    Class that takes the order mask boundaries, and uses them to create a FITS image where the
    order masks are visible.
    """

    def __init__(self):
        pass




    def run(self, maskBoundaryPath, outputFilePath=None):
        """
        Create the FITS file containing the image with the order masks visible 

        Input:
            maskBoundaryPath: path to FITS file containing the order boundaries
            outputFilePath:   If None, nothing is saved. Otherwise, a string with the name of the outputfile,
                              incl. the extension ".fits".

        Output:
            maskImage: 2D numpy array:   [ADU]
        """

        # Get the mask boundaries

        maskBoundaries = fits.open(maskBoundaryPath)

        # The lower_column contains for each row the leftmost column and upper_column the rightmost column.

        lower_cols = maskBoundaries[1].data 
        upper_cols = maskBoundaries[2].data 

        maskImage = np.zeros((10560, 10560))      # FIXME: don't hardcode the dimensions of the CCD

        for order in range(len(lower_cols)):          # all orders, of every fiber
            lower_col = lower_cols[order]
            upper_col = upper_cols[order]
            for row in range(len(lower_col)):
                for col in range(lower_col[row], upper_col[row]):
                    maskImage[row, col] = 1
        
        # If required, save to a FITS file

        if outputFilePath is not None:
            num_row, num_col = maskImage.shape
            hdr = fits.Header()
            hdr["rows"] = num_row
            hdr["cols"] = num_col
            outputParentPath = Path(outputFilePath).parent.absolute()
            if not os.path.exists(outputParentPath):
                os.makedirs(outputParentPath)
            hdu = fits.PrimaryHDU(maskImage, header=hdr)
            hdu.writeto(outputFilePath, overwrite=True)

        # That's it!

        return maskImage









if __name__ == "__main__":

    t1 = time.time()

    params   = yaml.safe_load(open("params.yaml"))

    root     = (params["Configuration"])["rootFolder"]
    mask_visualization_params = params["TwoDimensionalOrderMaskVisualisation"]

    inputPath = root + mask_visualization_params["inputPath"]
    outputPath = root + mask_visualization_params["outputPath"]

    maskImageCreator = MaskImageCreator()
    maskImageCreator.run(inputPath, outputPath)

    t2 = time.time()

    print(f"[{t2-t1}]")

