from pipeline import PipelineComponent
from pymongo import MongoClient
from numba import njit, jit, vectorize
from tqdm import tqdm
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import hashlib
import os
import tools
import debug


class OptimalExtraction(PipelineComponent):
    """
    Docstring
    """


    def __init__(self, input, debug=0):
        """
        Initializes the optimal extraction component.
        """

        super().__init__(input)
        self.col = self.setCollection(input)
        self.debug = debug
        self.outputPath = os.getcwd() + "/Data/ProcessedData/OptimalExtraction/"

        if not os.path.isdir(self.outputPath):
            os.mkdir(self.outputPath)

        self.type = "Optimal Extracted"


    def setCollection(self, input):
        """
        returns the collection in the database where the input hashes should be.
        """
        return self.db["ExtractedOrders"]


    def getFileName(self):
        return "optimal_extracted_flux.fits"



    def checkInput(self, input):
        """
        make sure that the input hashes correspond to correct image types.
        """

        # Check that the input only consist of two hashes
        inputAmount = len(input) == 2
        if not inputAmount:
            print("Input is not in the currect format, we expect a list with 2 elements, instead {} were fiven.".format(len(input)))
            return

        isCorrect = []
        collection = self.db["ExtractedOrders"]
        self.inputDict = {}
        
        for hash in input:
            instances = collection.find({"_id" : hash})
            isFromFlat = [( x["type"] == "Extracted Flat Orders") for x in instances]

            instances = collection.find({"_id" : hash})
            isFromScience = [( x["type"] == "Extracted Science Orders") for x in instances]

            if np.all(isFromFlat):
                if not len(isFromFlat) == 1:
                    return False
                self.inputDict["flat"] = hash
            elif np.all(isFromScience):
                if not len(isFromScience) == 1:
                    return False
                self.inputDict["science"] = hash
            else:
                return False
        return True



    def make(self):
        """
        run through the steps for optimal extraction
        """

        getPath = lambda x : ([x["path"] for x in self.col.find({"_id" : self.inputDict[x]})])[0]

        # Obtain and return the optimal extracted spectrum 
        return self.extractSpectrum(getPath("flat"), getPath("science"))



    def extractSpectrum(self, flatPath, sciencePath):

        fiberAndOrder   = []
        sSpectra = []
        fSpectra = []
        oSpectra = []

        fibers, orders = tools.getFibersAndOrders(flatPath)

        flatInfo = tools.getAllExtractedInfo(flatPath)
        scienceInfo = tools.getAllExtractedInfo(sciencePath)
        
        for o in tqdm(orders):
            for f in fibers:
                xFlat, yFlat, flat = (flatInfo[o])[f]
                xPosition, yPosition, science = (scienceInfo[o])[f]

                # We check that for the same order/fiber that the extracted positions for the science and flat image are the same.
                if not np.all(xFlat == xPosition) and np.all(yFlat == yPosition):
                    print("Order: {}, fiber: {} do not have matching coordinates for flat {} and science {}".format(o, f, flatPath, sciencePath))

                if f == 1:
                    continue
                flats, science, optim = getSpectrum(science, flat, xPosition, yPosition)

                fiberAndOrder.append((o, f))
                fSpectra.append(flats)
                sSpectra.append(science)
                oSpectra.append(optim)

        if self.debug > 2:
            debug.plotOrdersWithSlider(fSpectra, yMax=20000)
            debug.plotOrdersWithSlider(sSpectra, yMax=2000)
            debug.plotOrdersWithSlider(oSpectra, yMax=0.4)

        return oSpectra, fiberAndOrder



    def runComponent(self):
        """
        ...
        """
        spectrum, orders = self.make()
        self.saveImage(spectrum, orders)
        print("Block Generated!")



    def saveImage(self, spectrum, orders):
        """
        Save the image and add it to the database
        """
        hash = hashlib.sha256(bytes("".join(self.input), 'utf-8' )).hexdigest()
        path = self.outputPath + self.getFileName()
        orders, fibers = zip(*orders)

        # Save Optimal Extracted as FITS file
        primary_hdr = fits.Header()
        primary_hdr["hash"] = hash
        primary_hdr["path"] = path
        primary_hdr["type"] = self.type
        primary_hdr["orders"] = str(set(orders))
        primary_hdr["fibers"] = str(set(fibers))
        primary_hdr["input"] = str(self.input)

        hdu = fits.PrimaryHDU(header=primary_hdr)
        hdul = fits.HDUList([hdu])

        for i in np.arange(len(spectrum)):

            hdr1 = fits.Header()
            hdr1["order"] = orders[i]
            hdr1["fiber"] = orders[i]

            spect1 = np.array(spectrum[i], dtype=np.float64)
            col = fits.Column(name="Spectrum", format='D', array=spect1)

            cols = fits.ColDefs([col])
            hdu1 = fits.BinTableHDU.from_columns(cols, header=hdr1)

            hdul.append(hdu1)

        hdul.writeto(path, overwrite=True)

        # Add image to the database
        dict = {"_id" : hash, "path" : path, "type" : self.type}
        tools.addToDataBase(dict, overWrite = True)
        

        

        
     
@njit()
def getSpectrum(sFlux, fFlux, xPos, yPos, readout=2000):

    flats = np.zeros_like(np.unique(xPos), dtype=np.float32)
    scien = np.zeros_like(np.unique(xPos), dtype=np.float32)
    optim = np.zeros_like(np.unique(xPos), dtype=np.float32)

    for i, x in enumerate(np.unique(xPos)):
        mask = (xPos == x)
        signal = sFlux[mask]
        flat   = fFlux[mask]

        w = 1/(readout*np.ones_like(signal) + signal)

        flats[i] = np.sum(flat)/len(flat)
        scien[i] = np.sum(signal)/len(signal)
        if not np.sum(w * flat * flat)  == 0:
            optim[i] = np.sum(w * flat * signal) / np.sum(w * flat * flat)
        else:
            optim[i] = 1

    return flats, scien, optim



            


        

        




        

        



if __name__ == "__main__":

    # Extracted flat <-> Extracted science
    hash_list = ["2133e1778dd6f8208763eeeed3a5ae6efd525fabb33a5fdf8809bd77bf89bb2b", "626de973a22fe042b8355bfdb868260e1dc13cbbc411f4baaf6730b813e3a26d"]
    oExtracted = OptimalExtraction(hash_list, debug=2)
    oExtracted.runComponent()
