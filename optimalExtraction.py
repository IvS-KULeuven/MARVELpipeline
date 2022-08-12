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


    def __init__(self, input, extractedType, debug=0):
        """
        Initializes the optimal extraction component.
        """
        self.exType = extractedType
        self.inputDict = {}
        super().__init__(input)
        self.col = self.setCollection(input)
        self.debug = debug
        self.outputPath = os.getcwd() + "/Data/ProcessedData/OptimalExtraction/"

        if not os.path.isdir(self.outputPath):
            os.mkdir(self.outputPath)

        self.type = "Optimal Extracted {}".format(extractedType)




    def setCollection(self, input):
        """
        returns the collection in the database where the input hashes should be.
        """
        return self.db["ExtractedOrders"]



    def checkInput(self, input):
        """
        This function checks that the input has the right format. The input should consist of
        an image with Extracted Flat Orders and another image with Extracted Orders that we
        want to optimaly extract.        
        """

        # 1. Check that the input only consist of two hashes
        inputAmount = len(input) == 2
        if not inputAmount:
            print("Input is not in the correct format, we expect a list with 2 elements, instead {} were given.".format(len(input)))
            return False

        # 2. We check that one of the hashes corresponds to the Extracted Flat Orders.
        collection = self.db["ExtractedOrders"]
        amountOfExtractedFlats = 0
        for hash in input:
            isExtractedFlatOrder = len([ x for x in collection.find({"_id" : hash})]) == 1

            if isExtractedFlatOrder:
                isExtractedFlat = np.all([x["type"] == "Extracted Flat Orders" for x in collection.find({"_id": hash})])
                if isExtractedFlat:
                    amountOfExtractedFlats += 1
                    self.inputDict["flat"] = hash
                    
        return amountOfExtractedFlats == 1



    def make(self):
        """
        run through the steps for optimal extraction
        """

        getPath = lambda x : ([x["path"] for x in self.col.find({"_id" : self.inputDict[x]})])[0]

        # Obtain and return the optimal extracted spectrum 
        return self.extractSpectrum(getPath("flat"), getPath(self.exType.lower()))



    def extractSpectrum(self, flatPath, otherPath):

        fiberAndOrder   = []
        iSpectra  = []
        fSpectra  = []
        oSpectra  = []
        otSpectra = []
        pixels    = []

        fibers, orders = tools.getFibersAndOrders(flatPath)
        flatInfo = tools.getAllExtractedInfo(flatPath)
        otherInfo = tools.getAllExtractedInfo(otherPath)
        
        for o in tqdm(orders):
            for f in fibers:
                xFlat, yFlat, flat = (flatInfo[o])[f]
                xPosition, yPosition, other = (otherInfo[o])[f]

                # We check that for the same order/fiber that the extracted positions for the other and flat image are the same.
                if not np.all(xFlat == xPosition) and np.all(yFlat == yPosition):
                    print("Order: {}, fiber: {} do not have matching coordinates for flat {} and {} {}".format(o, f, flatPath, self.exType, otherPath))

                flats, other, optim = getSpectrum(other, flat, xPosition, yPosition)

                fiberAndOrder.append((o, f))
                fSpectra.append(flats)
                otSpectra.append(other)
                oSpectra.append(optim)
                pixels.append(np.unique(xPosition))

        if self.debug > 2:
            debug.plotOrdersWithSlider(fSpectra, xValues=pixels, yMax=5000)
            debug.plotOrdersWithSlider(otSpectra, xValues=pixels, yMax=3000)
            debug.plotOrdersWithSlider(oSpectra, xValues=pixels, yMax=2)

        if self.exType == "Science":
            return oSpectra, fiberAndOrder, pixels
        elif self.exType == "Etalon":
            print("SUCCES")
            return otSpectra, fiberAndOrder, pixels
        else:
            return None
        


                                          

    def run(self):
        """
        ...
        """
        spectrum, orders, pixels = self.make()
        self.saveImage(spectrum, orders, pixels)
        print("Block Generated!")



    def saveImage(self, spectrum, orders, pixels):
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
            hdr1["fiber"] = fibers[i]

            spect1 = np.array(spectrum[i], dtype=np.float64)
            col1 = fits.Column(name="Spectrum", format='D', array=spect1)

            pixels1 = np.array(pixels[i], dtype=np.int16)
            col2 = fits.Column(name="Pixels", format='J', array=pixels1)

            cols = fits.ColDefs([col1, col2])
            hdu1 = fits.BinTableHDU.from_columns(cols, header=hdr1)

            hdul.append(hdu1)

        hdul.writeto(path, overwrite=True)

        # Add image to the database
        dict = {"_id" : hash, "path" : path, "type" : self.type}
        tools.addToDataBase(dict, overWrite = True)
        

        

        
     
@njit()
def getSpectrum(oFlux, fFlux, xPos, yPos, readout=2000):

    flats = np.zeros_like(np.unique(xPos), dtype=np.float32)
    other = np.zeros_like(np.unique(xPos), dtype=np.float32)
    optim = np.zeros_like(np.unique(xPos), dtype=np.float32)

    for i, x in enumerate(np.unique(xPos)):
        mask = (xPos == x)
        signal = oFlux[mask]
        flat   = fFlux[mask]

        w = 1/(readout*np.ones_like(signal) + signal)

        flats[i] = np.sum(flat)/len(flat)
        other[i] = np.sum(signal)/len(signal)
        if not np.sum(w * flat * flat)  == 0:
            optim[i] = np.sum(w * flat * signal) / np.sum(w * flat * flat)
        else:
            optim[i] = 1

    return flats, other, optim



            


        

        




        

        



if __name__ == "__main__":

    # Extracted flat <-> Extracted science
    hash_list = ["2133e1778dd6f8208763eeeed3a5ae6efd525fabb33a5fdf8809bd77bf89bb2b", "626de973a22fe042b8355bfdb868260e1dc13cbbc411f4baaf6730b813e3a26d"]
    oExtracted = OptimalExtraction(hash_list, debug=3)
    oExtracted.run()
