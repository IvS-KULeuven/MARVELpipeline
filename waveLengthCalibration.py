import os
import h5py
import pandas as pd
import numpy as np
import debug
from pipeline import PipelineComponent
from numba import njit, jit, vectorize
import matplotlib.pyplot as plt
from astroquery.nist import Nist
import astropy.units as u
import tools


class WaveLengthCalibration(PipelineComponent):
    def __init__(self, input, debug=0):
        """
        Initialize the component
        """
        self.inputDict = {}
        super().__init__(input)
        self.col = self.setCollection(input)
        self.debug = debug
        self.outputPath = os.getcwd() + "/Data/ProcessedData/WaveLengthCalibration/"

        if not os.path.isdir(self.outputPath):
            os.mkdir(self.outputPath)

        self.type = "WaveLength Calbrated Orders"

        # Load in the reference line list
        self.df = pd.read_csv("Data/ReferenceLineList.csv")
        self.lines = pd.read_csv("Data/ThAr_linelist_NIST_Nave2017.csv")
        self.lines = self.lines[self.lines["Ritz_wavelength"] != "-"]





    def setCollection(self, input):
        """
        Where in the database we should look for the input hashes
        """
        return self.db["OptimalExtracted"]



    def checkInput(self, input):
        """
        This function checks that the input has the right format. This input should consist of
        an image with Optimal Extracted Etalon.
        """

        # 1. Check that the input only consist of one hash
        inputAmount = len(input) == 1
        if not inputAmount:
            print("Input is not in the correct format, we expect only one hash, instead {} were given".format(len(input)))
            return False
                  

        collection = self.db["OptimalExtracted"]

        # 2. We check that one of the hashes corresponds to the Optimal Extracted Etalon
        amountOfEtalon = 0
        for hash in input:
            isOptimalExtracted = len([x for x in collection.find({"_id": hash})]) == 1

            if isOptimalExtracted:
                isOptimalExtractedEtalon = np.all([x["type"] == "Optimal Extracted Etalon" for x in collection.find({"_id": hash})])
                if isOptimalExtractedEtalon:
                    amountOfEtalon += 1
                    self.inputDict["etalon"] = hash

        return amountOfEtalon == 1



    def test1(self):
        etalonPath = [ x["path"] for x in self.db["OptimalExtracted"].find({"_id": self.inputDict["etalon"]})][0]
        fibers, orders = tools.getFibersAndOrders(etalonPath)


        # for i in np.arange(5, 667, 5):
        #     position, flux = (tools.getOptimalExtractedSpectrum(etalonPath)[3])[5]
        #     plt.plot(position, flux)
        #     plt.show()
        



if __name__ == "__main__":

    etalon_hash = ["05510755d260ae30f45b8a7742141595b327cb01c13b970c03a5780ea254a891"]
    
    calibrated = WaveLengthCalibration(etalon_hash)
    calibrated.test1()
