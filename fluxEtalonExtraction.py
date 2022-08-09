from fluxExtraction import FluxExtraction
import numpy as np
import os
import tools
from numba import njit, jit, vectorize
import matplotlib.pyplot as plt
from astropy.io import fits
import hashlib

class EtalonOrderExtraction2(FluxExtraction):
    """
    DESCRIPTION: This modeule uses the mask obtained from the flat field extraction to extract
                 the flux from the Etalon calibration images for every fiber and order.

    INPUT: Input hashes corresponding to one Extracted Flat Flux image and one Master etalon image.

    ALGORITHM:
    1. From the Extract Flat Flux image we are able to obtain a mask of pixels that should be extracted.
    2. We extract the pixels in within this mask for every fiber/order and save this in a fits file.
    """

    def __init__(self, input, debug=0):
        """
        Initialize the component
        """
        super().__init__(input, "Etalon", debug=debug)


    def getFileName(self):
        """
        Retun the name of the output file
        """
        return "extracted_etalon_orders.fits"



    def checkInput(self, input):
        """
        This function checks that the input has the right format. The input should consist of
        an etalon image from which we want to extract the flux and an image with extracted flat orders.
        """
        if not super().checkInput(input):
            return False

        etalonFrames = self.db["EtalonImages"]
        amountOfEtalonOrders = 0
        for hash in input:
            isEtalon = len([x for x in etalonFrames.find({"_id": hash})]) == 1

            if isEtalon:
                isEtl = np.all([ x["type"] == "Raw Etalon Image" for x in etalonFrames.find({"__id": hash})])
                if isEtl:
                    amountOfEtalonOrders += 1

        return amountOfEtalonOrders == 1


if __name__ == "__main__":
    
    # Raw Etalon Image <-> Extracted Flat Orders
    hash = ["e0ac021d19ce5520d0ba92df1d5aadf6541f35f76a84121576828287937ca508", "2133e1778dd6f8208763eeeed3a5ae6efd525fabb33a5fdf8809bd77bf89bb2b"]
    FExtra = EtalonOrderExtraction2(hash, debug=3)

    FExtra.run()



    
                
            
