from fluxExtraction import FluxExtraction
import numpy as np
import os
import tools
from numba import njit, jit, vectorize
import matplotlib.pyplot as plt
from astropy.io import fits
import hashlib

class ScienceOrderExtraction2(FluxExtraction):
    """
    DESCRIPTION: This modeule uses the mask obtained from the flat field extraction to extract
                 the flux from the Calibrated Science images for every fiber and order.

    INPUT: Input hashes corresponding to one Extracted Flat Flux image and one Calibrated Science Image.

    ALGORITHM:
    1. From the Extract Flat Flux image we are able to obtain a mask of pixels that should be extracted.
    2. We extract the pixels in within this mask for every fiber/order and save this in a fits file.
    """

    def __init__(self, input, debug=0):
        """
        Initialize the component
        """
        super().__init__(input, "Science", debug=debug)



    def getFileName(self):
        """
        Return the name of the output file
        """
        return "extracted_science_orders.fits"



    def checkInput(self, input):
        """
        This function checks that the input has the right format. The input should consist of
        a calibrated science image from which we want to extract the flux and an image with extracted flat orders.
        """
        if not super().checkInput(input):
            return False

        scienceFrames = self.db["ScienceImages"]
        amountOfScienceOrders = 0
        for hash in input:
            isScience = len([x for x in scienceFrames.find({"_id": hash})]) == 1
            if isScience:
                isCalibratedScience = np.all([ x["type"] == "Calibrated Science Image" for x in scienceFrames.find({"_id": hash})])
                if isCalibratedScience:
                    amountOfScienceOrders += 1
        return amountOfScienceOrders == 1





if __name__ == "__main__":
    hash = ["b0ef6a99bde7cdbc968a46fcd7a57e450a554c548d9cc89d7a9555e7236fe05f", "2133e1778dd6f8208763eeeed3a5ae6efd525fabb33a5fdf8809bd77bf89bb2b"]

    FExtra = ScienceOrderExtraction2(hash, debug=3)
    FExtra.run()

 
