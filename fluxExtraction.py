from pipeline import PipelineComponent
import numpy as np
import os
import tools
from numba import njit, jit, vectorize
import matplotlib.pyplot as plt
from astropy.io import fits
import hashlib

class FluxExtraction(PipelineComponent):
    """
    This is the parent class for fluxScienceExtraction and and fluxEtalonExtraction.
    Both child classes follow a similar alghoritms. The idea behind this parent class
    is to implement methods that can easily be used by both classes.
    """

    def __init__(self, input, extractedType, debug=0):
        """        
        Initialize the component
        """
        self.exType = extractedType
        super().__init__(input)
        self.col = self.setCollection(input)
        self.input = self.createInputDirectory(input)
        self.debug = debug
        self.outputPath = os.getcwd() + "/Data/ProcessedData/ExtractedOrders/"

        if not os.path.isdir(self.outputPath):
            os.mkdir(self.outputPath)

        self.type = "Extracted {} Orders".format(extractedType)



    def setCollection(self, input):
        """
        Where in the database should we look for the input hashes
        """
        return { self.exType.lower(): self.db[self.exType + "Images"], "extracted": self.db["ExtractedOrders"]}



    def createInputDirectory(self, input):
        """
        Create a dictionary to distinguish between hashes
        """
        sortedInput = { self.exType.lower(): [], "extracted": []}
        for hash in input:
            if len([x for x in self.db[self.exType + "Images"].find({"_id" : hash})]) == 1:
                sortedInput[self.exType.lower()].append(hash)
            elif len([x for x in self.db["ExtractedOrders"].find({"_id" : hash})]) == 1:
                sortedInput["extracted"].append(hash)
        return sortedInput



    def checkInput(self, input):
        """
        This function checks that the input has the right format. The input should consist of
        an image from which we want to extract the flux and an image with extracted flat orders.
        This function should be specified for every child class, for this parent class it always
        returns True.
        """

        # 1. We start by checking that the input has the right amount of elements 
        inputAmount = (len(input)==2)
        if not inputAmount:
            print("Input is not in the correct format, we expect a list with 2 elements, instead {} were given".format(len(iput)))
            return False

        # 2. We check that one of the hashes corresponds to the extracted flat orders
        extractedOrders = self.db["ExtractedOrders"]
        amountOfFlatOrders = 0
        for hash in input:
            isExtractedOrders = len([x for x in extractedOrders.find({"_id": hash})]) == 1

            if isExtractedOrders:
                isExtractedFlat = np.all([x["type"] == "Extracted Flat Orders" for x in extractedOrders.find({"_id": hash})])
                if isExtractedFlat:
                    amountOfFlatOrders += 1

        return amountOfFlatOrders == 1


        

    def make(self):
        """
        Runs the algorithm as described in the class definition.
        """
        # 1. We get the image of the Image Frame
        imagePath = [[ x["path"] for x in self.col[self.exType.lower()].find({"_id" : hash})] for hash in self.input[self.exType.lower()]]
        image = tools.getImage((imagePath[0])[0]).astype('float64')

        # 2. We obtain the mask from the Extracted Flat Orders
        mask = self.getMaskFromExtractedFlat()

        # 3. Using this mask we are able to obtain the relavant image
        return self.extractStripes(image, mask)



    def run(self):
        xCoordinates, yCoordinates, fluxValues, orders = self.make()

        self.saveImage(xCoordinates, yCoordinates, fluxValues, orders)
        print("Block Generated")



    def getMaskFromExtractedFlat(self):
        extractedFlatPath = [[x["path"] for x in self.col["extracted"].find({"_id" : hash})] for hash in self.input["extracted"]]
        extractedFlatPath = (extractedFlatPath[0])[0]
        fibers, orders = tools.getFibersAndOrders(extractedFlatPath)

        mask = {}
        for f in fibers:
            for o in orders:
                positions = tools.getExtractedPosition(extractedFlatPath, o, f)
                if f in mask:
                    mask[f].update({o : positions})
                else:
                    mask[f] = {o : positions}
        return mask


    def extractStripes(self, image, mask):
        fluxValues   = []
        orders       = []
        xCoordinates = []
        yCoordinates = []
        
        mask_image = np.zeros(image.shape)
        for f in mask.keys():
            for o in mask[f].keys():
                positions = (mask[f])[o]
                xPositions, yPositions = zip(*positions)
                xPositions = np.array(xPositions)
                yPositions = np.array(yPositions)

                for x, y in positions:
                    mask_image[x, y] = 1
                    
                fluxValues.append([ image[x, y] for x,y in zip(xPositions, yPositions)])
                orders.append((o,f))
                xCoordinates.append(xPositions)
                yCoordinates.append(yPositions)
                    
        if self.debug > 2:
            plt.imshow(mask_image, origin='lower')
            plt.imshow(image, alpha=0.95, origin='lower')
            plt.show()

        return xCoordinates, yCoordinates, fluxValues, orders



    def saveImage(self, xValues, yValues, flux, orders):
        """
        Save the image to a fits file and add the file to the database
        """
        hash = hashlib.sha256(bytes("".join(self.input), 'utf-8')).hexdigest()
        path = self.outputPath + self.getFileName()
        orders, fibers = zip(*orders)
        # Save Extracted Flat Orders as FITS file for primary HDU
        primary_hdr = fits.Header()
        primary_hdr["hash"] = hash
        primary_hdr["path"] = path
        primary_hdr["type"] = self.type
        primary_hdr["orders"] = str(set(orders))
        primary_hdr["fibers"] = str(set(fibers))
        primary_hdr["input"] = str(self.input)


        hdu = fits.PrimaryHDU(header=primary_hdr)
        hdul = fits.HDUList([hdu])

        for i in np.arange(len(flux)):

            hdr1 = fits.Header()
            hdr1["order"]= orders[i]
            hdr1["fiber"]= fibers[i]            

            xValue = np.array(xValues[i], dtype=np.int16)
            col1 = fits.Column(name="X", format='J', array=xValue)

            yValue = np.array(yValues[i], dtype=np.int16)
            col2 = fits.Column(name="Y", format='J', array=yValue)

            fluxValues = np.array(flux[i], dtype=np.float64)
            col3 = fits.Column(name="flux", format='D', array=fluxValues)

            cols = fits.ColDefs([col1, col2, col3])
            hdu1 = fits.BinTableHDU.from_columns(cols, header=hdr1)

            hdul.append(hdu1)

        hdul.writeto(path, overwrite=True)

        # Add image to the database
        dict = {"_id" : hash, "path" : path, "type" : self.type}
        tools.addToDataBase(dict, overWrite = True)        


