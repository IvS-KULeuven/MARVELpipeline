from pipeline import PipelineComponent
from pymongo import MongoClient
from numba import njit, jit, vectorize, prange
from tqdm import tqdm
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import hashlib
import os
import tools
import debug
import time


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
        
        fibers, orders = tools.getFibersAndOrders(flatPath)

        flatInfo = tools.getAllExtractedSpectrum(flatPath)
        otherInfo = tools.getAllExtractedSpectrum(otherPath)

        xValues, yValues, fValuesF, fValuesO, self.orders, self.fibers = self.convertToArray(otherInfo, flatInfo, fibers, orders)

        print("Start getOptimalspectrum")
        start = time.time()        
        otSpectra, oSpectra, yPixels, xPixels = getOptimalSpectrum(xValues, yValues, fValuesF, fValuesO)
        end   = time.time()
        print("\tTime: ", end-start, "s")

    
        if self.debug > 2:
            #debug.plotOrdersWithSlider(fSpectra, xValues=xPixels, yMax=5000)
            debug.plotOrdersWithSlider(otSpectra, xValues=xPixels, yMax=3000)
            debug.plotOrdersWithSlider(oSpectra, xValues=xPixels, yMax=2)
        
        if self.exType == "Science":
            return oSpectra, xPixels, yPixels
        elif self.exType == "Etalon":
            return otSpectra, xPixels, yPixels
        else:
            return None
        



    def run(self):
        """
        ...
        """
        spectrum, xPixels, yPixels = self.make()
        self.saveImage(spectrum, xPixels, yPixels)
        print("Block Generated!")



    def saveImage(self, spectrum, xPixels, yPixels):
        """
        Save the image and add it to the database
        """
        hash = hashlib.sha256(bytes("".join(self.input), 'utf-8' )).hexdigest()

        path = self.outputPath + self.getFileName()
        fibers = np.unique(self.fibers)
        orders = np.unique(self.orders)
                  

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
            hdr1["order"] = self.orders[i]
            hdr1["fiber"] = self.fibers[i]

            spect1 = np.array(spectrum[i], dtype=np.float64)
            col1 = fits.Column(name="Spectrum", format='D', array=spect1)

            pixels1 = np.array(xPixels[i], dtype=np.int16)
            col2 = fits.Column(name="xPixels", format='J', array=pixels1)

            pixels2 = np.array(yPixels[i], dtype=np.int16)
            col3 = fits.Column(name="yPixels", format='J', array=pixels2)

            cols = fits.ColDefs([col1, col2, col3])
            hdu1 = fits.BinTableHDU.from_columns(cols, header=hdr1)

            hdul.append(hdu1)

        hdul.writeto(path, overwrite=True)

        # Add image to the database
        dict = {"_id" : hash, "path" : path, "type" : self.type}
        tools.addToDataBase(dict, overWrite = True)
        


    def convertToArray(self, otherInfo, flatInfo, fibers, orders):
        """
        This methods converts the list that is given as input into a numpy array.
        This is done so that numba can be used to optimize the optimal extraction.
        """
        start = time.time()
        print("start convert to array")
        max_row_length = 0
        # Let us first get the maximum length of the rows 
        x_values = []
        y_values = []
        f_values = []
        o_values = []
        ordrs    = []
        fbrs     = []

        for o in orders:
            for f in fibers:
                x, y, flux = (flatInfo[o])[f]
                if len(x) > max_row_length:
                    max_row_length = len(x)

        for o in orders:
            for f in fibers:
                x1, y1, otherFlux = (otherInfo[o])[f]
                x2, y2, flatFlux = (flatInfo[o])[f]

                if not ((x1 == x2).all() and (y1 == y2).all()):
                    print("Flat mask does not correspond with other mask")
                    exit
                x = list(x1)
                y = list(y1)
                otherFlux = list(otherFlux)
                flatFlux  = list(flatFlux)
                

                while (len(x) < max_row_length):
                    x.append(-1)
                    y.append(-1)
                    otherFlux.append(np.nan)
                    flatFlux.append(np.nan)

                x_values += [x]
                y_values += [y]
                f_values += [flatFlux]
                o_values += [otherFlux]
                ordrs.append(o)
                fbrs.append(f)
        end = time.time()
        print("\tTime: ", end-start, "s")
        return np.array(x_values, dtype=np.int16), np.array(y_values, dtype=np.int16), np.array(f_values, np.float64), np.array(o_values, np.float64), np.array(ordrs), np.array(fbrs)

        




        

        
     
@njit()
def getSpectrum(oFlux, fFlux, xPos, yPos, readout=2000):

    # Initialize the arrays that will be added to the output array
    others = np.ones(10560, dtype=np.float32) * -1
    optim  = np.ones(10560, dtype=np.float32) * -1
    yPositionAvg = np.ones(10560, dtype=np.float32) * -1
        
    # Create the optimal order for every xValue
    for i,x in enumerate(np.unique(xPos)):
        if x == -1:
            others[i] = np.nan
            optim[i]  = np.nan
            yPositionAvg[i] = np.nan
            continue

        # Select the pixels that correspond to the xValue
        mask = (xPos == x)
        
        signal = (oFlux)[mask]
        flat   = (fFlux)[mask]
        yCoord = (yPos)[mask]

        w = 1/(readout*np.ones_like(signal) + signal)

        others[i] = np.sum(signal)/len(signal)
        if np.sum(signal) == 0:
            yPositionAvg[i] = np.sum(yCoord) / len(yCoord)
        else:
            yPositionAvg[i] = np.sum(yCoord * signal) / np.sum(signal)
        if not np.sum(w * flat * flat)  == 0:
            optim[i] = np.sum(w * flat * signal) / np.sum(w * flat * flat)
        else:
            optim[i] = 1
    
    return others, optim, yPositionAvg


@njit(parallel=True)
def getOptimalSpectrum(xValues, yValues, flat, other, readout=2000):
    # Initialize the output arrays
    otSpectra = np.empty((66 * 5, 10560))
    oSpectra  = np.empty((66 * 5, 10560))
    yPixels   = np.empty((66 * 5, 10560))
    xPixels   = np.empty((66 * 5, 10560))
    
    # Loop over every order/fiber 
    for line in prange(66*5):
        others, optim, yPositionAvg = getSpectrum(other[line], flat[line], xValues[line], yValues[line])
        # Add the optimal spectrum, avg flux spectrum and average Y value to the output arrays
        otSpectra[line] = others
        oSpectra[line]  = optim
        yPixels[line]   = yPositionAvg
        (xPixels[line])[:len(np.unique(xValues[line]))] = np.unique(xValues[line])

    return otSpectra, oSpectra, yPixels, xPixels
