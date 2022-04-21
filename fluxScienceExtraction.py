from pipeline import PipelineComponent
import numpy as np
import os
import tools
from numba import njit, jit, vectorize
import matplotlib.pyplot as plt
from astropy.io import fits
import hashlib


class ScienceOrderExtraction(PipelineComponent):
    """
    Docstring
    """
    
    def __init__(self, input, debug=0):
        """
        Docstring
        """
        super().__init__(input)
        self.col = self.setCollection(input)
        self.input = self.createInputDirectory(input)
        self.debug = debug
        self.outputPath = os.getcwd() + "/Data/ProcessedData/ExtractedOrders/"

        if not os.path.isdir(self.outputPath):
            os.mkdir(self.outputPath)

        self.type = "Extracted Science Orders"

    def setCollection(self, input):
        """
        Docstring
        """
        return { "science": self.db["ScienceImages"], "extracted": self.db["ExtractedOrders"]}

    def createInputDirectory(self, input):
        sortedInput = {"science" : [], "extracted" : []}
        for hash in input:
            if len([x for x in self.db["ScienceImages"].find({"_id" : hash})]) == 1:
                sortedInput["science"].append(hash)
            elif len([x for x in self.db["ExtractedOrders"].find({"_id" : hash})]) == 1:
                sortedInput["extracted"].append(hash)
        return sortedInput

    
    def checkInput(self, input):
        """
        Checks that the provided input has the right format. 
        The input should consist out of a calibrated 
        science frame, and one with extracted flat orders. 
        """
        inputAmount = (len(input) == 2)
        if not inputAmount:
            print("Input is not in the correct format, we expect a list with 2 elements, instead {} were given.".format(len(input)))
            return

        scienceFrames   = self.db["ScienceImages"]
        extractedOrders = self.db["ExtractedOrders"]


        for hash in input:
            isScience = len([x for x in scienceFrames.find({"_id" : hash})]) == 1
            isExtractedFOrder = len([x for x in extractedOrders.find({"_id" : hash})]) == 1

            if isScience:
                isCalibrated = np.all([ x["type"] == "Calibrated Science Image" for x in scienceFrames.find({"_id" : hash})])
                if not isCalibrated:
                    print("Science image with hash {} is not Calibrated Image".format(hash))
                    return False
            elif isExtractedFOrder:
                isExtractedFlat = np.all([ x["type"] == "Extracted Flat Orders" for x in extractedOrders.find({"_id" : hash})])
                if not isExtractedFOrder:
                    print("Extracted order image with hash {} is not Extracted Flat Order".format(hash))
                    return False
            else:
                # Is not science of extracted order image
                return False
        return True

    def make(self):
        """
        Docstring
        """
        # 1. We should get the image of the SF
        sciencePath = [[x["path"] for x in self.col["science"].find({"_id" : hash})] for hash in self.input["science"]]
        image = tools.getImage((sciencePath[0])[0]).astype('float64')[0]
        
        
        # 2. We obtain the mask from the Extracted Flat Orders
        mask = self.getMaskFromExtractedFlat()

        # 3. From the mask we extract the stripes from the science image
        return self.extractStripes(image, mask)


    def runComponent(self):
        xCoordinates, yCoordinates, fluxValues, orders = self.make()
        self.saveImage(xCoordinates, yCoordinates, fluxValues, orders)
        print("Block Generated")



    def saveImage(self, xValues, yValues, flux, orders):
        """
        Docstring
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

            fluxValues = np.array(flux[i], dtype=np.int16)
            col3 = fits.Column(name="flux", format='D', array=fluxValues)

            cols = fits.ColDefs([col1, col2, col3])
            hdu1 = fits.BinTableHDU.from_columns(cols, header=hdr1)

            hdul.append(hdu1)

        hdul.writeto(path, overwrite=True)

        # Add image to the database
        dict = {"_id" : hash, "path" : path, "type" : self.type}
        tools.addToDataBase(dict, overWrite = True)


    def getFileName(self):
        """
        Docstring
        """
        return "extracted_science_orders.fits"

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


                    
                
        #mask = tools.getExtractedPosition((extractedFlatPath[0])[0], 1, 1)


                
                    




    


if __name__ == "__main__":
    
    # CalibratedScience <-> Extracted Flat Orders
    hash = ["b0ef6a99bde7cdbc968a46fcd7a57e450a554c548d9cc89d7a9555e7236fe05f", "f90954c06fb77bb3f5594bba9a4da19a98c0ae13bff8eb4cfe5bc30c88b26f66"]
    FExtra = ScienceOrderExtraction(hash, debug=2)
    FExtra.runComponent()    
