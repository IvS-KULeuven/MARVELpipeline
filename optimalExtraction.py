from pipeline import PipelineComponent
from pymongo import MongoClient
from numba import njit, jit, vectorize
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
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
        self.outputPath = os.getcwd() + "/Data/ProcessedData/OptimalExtraction"

        if not os.path.isdir(self.outputPath):
            os.mkdir(self.outputPath)

        self.type = "Optimal Extracted"


    def setCollection(self, input):
        """
        returns the collection in the database where the input hashes should be.
        """
        return self.db["ExtractedOrders"]

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
        
        
        spectrum = self.extractSpectrum(getPath("flat"), getPath("science"))
       


    def extractSpectrum(self, flatPath, sciencePath):

        sSpectra = []
        fSpectra = []
        oSpectra = []
        

        for o in tqdm(np.arange(1, 67)):
            for f in np.arange(1, 6):
                if f == 1:
                    continue;
                # Make sure that flatPositions and sciencePositions are the same 
                position = tools.getExtractedPosition(flatPath, o, f)
                haveSamePosition = np.all(position == tools.getExtractedPosition(sciencePath, o, f))
        
                science = tools.getExtractedFlux(sciencePath, o, f)
                flat    = tools.getExtractedFlux(flatPath, o, f)
                
                xPosition, yPosition = zip(*position)
                xPosition = np.array(xPosition)
                yPosition = np.array(yPosition)
            
                flats, scien, optim = getSpectrum(science, flat, xPosition, yPosition)
                sSpectra.append(scien)
                fSpectra.append(flats)
                oSpectra.append(optim)
                                
        debug.plotOrdersWithSlider(fSpectra, yMax=20000)
        debug.plotOrdersWithSlider(sSpectra, yMax=3000)
        debug.plotOrdersWithSlider(oSpectra, yMax=2)
        print("done")



        
     
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
    oExtracted = OptimalExtraction(hash_list)
    oExtracted.make()
