import os
import yaml
import tools
import time
from pathlib import Path
import numpy as np
from astropy.io import fits
from pipeline import PipelineComponent






class BiasAndDarkCorrectedScienceFrames(PipelineComponent):
    """
    Class that subtracts the electronic offset (bias) and the dark current from the science image.
    """

    def __init__(self):
        pass



    def run(self, rawScienceImagePaths, masterBiasPath, masterDarkPath, outputPaths=None):
        """
        Create the bias and dark subtracted science image

        Input:
            rawSciencePaths: list of strings, full paths to the 2D raw science fits file
            masterBiasPath:  string containing the full path to the 2D master bias fits file
            masterDarkPath:  string containing the full path to the 2D master dark fits file
            outputPaths:     If None, nothing is saved. Otherwise, a list of the same size
                             as rawScienceImagePath, containing the output filenames, incl
                             the extension ".fits".

        Output:
            biasAndDarkSubtractedScienceImages: list of 2D numpy arrays: bias and dark corrected science images [ADU]
        """

        # Get all the relevant bias and dark images

        masterBias = tools.getImage(masterBiasPath)
        masterDark = tools.getImage(masterDarkPath)

        # Extract the exposure time for which the dark current image was measured. 

        darkFileStem = Path(masterDarkPath).stem
        darkExposureTime = float(darkFileStem[-4:])
        
        # Compute a median bias and dark level

        darkLevel = np.median(masterDark)
        biasLevel = np.median(masterBias)

        # Verify if the output file list has the same size as the input file list

        if outputPaths is not None:
            if len(rawScienceImagePaths) != len(outputPaths):
                print("Error: BiasAndDarkCorrectedScienceFrames: list of input paths has not the same size as list of output paths")
                return None

        # Process each of the raw 2D science images

        biasAndDarkSubtractedScienceImages = []

        for n in range(len(rawScienceImagePaths)):
            scienceImage = tools.getImage(rawScienceImagePaths[n])
            scienceFileStem = Path(rawScienceImagePaths[n]).stem
            scienceExposureTime = float(scienceFileStem[-4:])

            scienceImage = scienceImage - biasLevel - darkLevel / darkExposureTime * scienceExposureTime

            # Zero all negative values

            scienceImage[scienceImage < 0] = 0

            # Convert float64 values in science images to int32 

            biasAndDarkSubtractedScienceImages.append(scienceImage.astype(np.int32))

            if outputPaths is not None:

                # Make sure that the output directories exist, if not create them

                outputDirectory = Path(outputPaths[n]).parent.absolute()
                if not os.path.exists(outputDirectory):
                    os.makedirs(outputDirectory)

                # Save the shape of the image in a header 

                num_row, num_col = scienceImage.shape
                hdr = fits.Header()
                hdr["rows"] = num_row
                hdr["cols"] = num_col

                # Save the the image to the FITS file

                hdu = fits.PrimaryHDU(biasAndDarkSubtractedScienceImages[-1], header=hdr)
                hdu.writeto(outputPaths[n], overwrite=True)


        # That's it!

        return biasAndDarkSubtractedScienceImages



