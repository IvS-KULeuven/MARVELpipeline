import os
import h5py
import pandas as pd
import numpy as np
import debug
from pipeline import PipelineComponent
from numba import njit, jit, vectorize
import matplotlib.pyplot as plt
import scipy.signal
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


        etalonSpectrum = tools.getOptimalExtractedSpectrum(etalonPath)
        position, flux = (etalonSpectrum[19])[5]
        self.doForFiberOrder(position, flux)



    def doForFiberOrder(self, position, flux):
        maxim, minim = self.getMaxAndMin(position, flux)


    def getMaxAndMin(self, position, flux):

        # Replace cosmics as nan values
        cleanmax = np.median(np.sort(flux)[-200:])
        flux[np.array([e>1.5*cleanmax for e in flux])] = np.nan

        # Replace the nan values by interpolated values
        dataClean = flux.copy()
        mask = np.isnan(flux)
        dataClean[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), flux[~mask])

        # Apply the savgol filter 
        d = scipy.signal.savgol_filter(dataClean, polyorder=1, window_length=3, mode="interp")

        maxima = scipy.signal.argrelextrema(d, np.greater_equal, order=2)[0]
        minima = scipy.signal.argrelextrema(d, np.less_equal, order=2)[0]

        def ranges(nums):
            gaps = [[s, e] for s, e in zip(nums, nums[1:]) if s + 1 < e]
            edges = iter(list(nums[:1]) + sum(gaps, []) + list(nums[-1:]))
            return list(zip(edges, edges))

        # clean possible double minima / maxima
        minima = np.array([ np.mean(r, dtype=int) for r in ranges(minima)]).flatten()
        maxima = np.array([ np.mean(r, dtype=int) for r in ranges(maxima)]).flatten()

        # drop incomplete peaks on either side of the spectrum
        left, right = minima[0], minima[-1]
        maxima = maxima[(maxima > left) & (maxima < right)]


        # make sure there  is a minimum inbetween each max
        k = np.asarray([np.any(np.logical_and(minima >= maxima[j], minima <= maxima[j+1])) for j in range(len(maxima) - 1)])
        missing = np.where(~k)
        for idx in missing:
            minima = np.append(minima, ((maxima[idx] + maxima[idx+1])/2).astype('int'))
        minima.sort()
        print(len(maxima), len(minima))


        # make sure there is a maxima between each minima
        k = np.asarray([(np.any(np.logical_and(maxima > minima[j], maxima <= minima[j+1]))) for j in range(len(minima) - 1)])
        missing = np.where(~k)
        for idx in missing:
            maxima = np.append(maxima, ((minima[idx] + minima[idx+1])/2).astype('int'))
        maxima.sort()
        print(len(maxima), len(minima))



        wrongMax = [ True if idx in minima else False for idx in maxima]
        for idx in maxima[wrongMax]:
            maxima = np.delete(maxima, idx)
            print("Wrong maxima found in pixel {} and removed".format(idx))


        if len(minima) != len(maxima) + 1 or np.any((minima[:-1] > maxima)|(maxima > minima[1:])):
            plt.figure()
            d_nonan = d.copy()
            plt.plot(d)
            plt.plot(maxima, d[maxima], 'g+')
            plt.plot(minima, d[minima], 'r+')
            #plt.hlines(0.25 * np.max(d), 0, len(d))
            plt.show()
            pass

        return maxima, minima



        # Make sure there is a max between each min
        # k = np.asarray([(np.any(np.logical_and(minima >= maxima[j], minima <= maxima[j + 1]))) for j in range(len(maxima) - 1)])
        # missing = np.where(~k)
        # for idx in missing:
        #     minima = np.append(minima,((maxima[idx]+maxima[idx+1])/2).astype('int'))
        # minima.sort()

        
        






if __name__ == "__main__":

    etalon_hash = ["05510755d260ae30f45b8a7742141595b327cb01c13b970c03a5780ea254a891"]
    
    calibrated = WaveLengthCalibration(etalon_hash)
    calibrated.test1()
