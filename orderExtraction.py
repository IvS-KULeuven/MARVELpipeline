import hashlib
import os
import tools

import matplotlib.pyplot as plt
import numpy as np

from numba      import njit, jit, vectorize
from astropy.io import fits
from datetime   import datetime
from pipeline   import PipelineComponent
from database   import DatabaseFromLocalFiles







class OrderExtraction(PipelineComponent):
    """
    Class that performs the extraction of the orders given an image and an order mask.
    """

    def __init__(self, database=None, debug=0, **imageAndMaskHash):
        """
        Initialize the order extraction component

        Input:
            database:      if not None, DatabaseFromLocalFiles object to be used as database.
                           else MongoDB is used as database.

            debug:         0, 1, 2, or 3:  0 meaning no debug output, 3 meaning lots of debug output.

            imageAndMaskHash: hashes of images to be used in the extraction.
                              Given as ImageType/ImageHash (keyword/argument)


        Output:
            None
        """

        super().__init__(database, **imageAndMaskHash)
        imageAndMaskHash = self.inputHashes

        if self.checkSanityOfInputTypes(**imageAndMaskHash):
            self.outputPath = os.getcwd() + "/Data/ProcessedData/ExtractedOrders/"
            self.outputType = "Extracted " + self.extractedType + " Orders"
            self.orderMaskHash = imageAndMaskHash["ExtractedOrders"]
            self.imageHash = imageAndMaskHash[self.extractedType + "Images"]
            self.debug = debug

        else:
            raise Exception("Error: The input hashes do not match the correct type: Aborting")
            exit(1)

        # +=============================================================+
        self.imageCollection = self.db[self.extractedType + "Images"]
        self.orderMaskCollection = self.db["ExtractedOrders"]
        # +=============================================================+

        if not os.path.isdir(self.outputPath):
            os.mkdir(self.outputPath)






    def checkSanityOfInputTypes(self, **imageAndMaskHash):
        """
        This function is ran after we run checkSanityOfInputHashes. This function checks the the
        input types that are given is able to generate a order extracted output file.
        """

        types  = list(imageAndMaskHash.keys())
        values = list(imageAndMaskHash.values())

        # Check that the keys are of the right format. For order extraction there should be one
        # extracted flat image and one calibrated etalon or science frame.

        hasMask = (len(types) == 2) and ("ExtractedOrders" in types)
        hasImages = ("ScienceImages" in types) or ("EtalonImages" in types)
        keysAreCorrect = hasMask and hasImages

        valuesAreCorrect = isinstance(values[0], str) and isinstance(values[1], str)


        # We can set the type of image, e.g. "Science" or "Etalon"
        if "ScienceImages" in types:
            self.extractedType = "Science"
        else:
            self.extractedType = "Etalon"

        if keysAreCorrect and valuesAreCorrect:
            hasExtractedFlat = self.db["ExtractedOrders"].find_one({"_id": imageAndMaskHash["ExtractedOrders"]})["type"] == "Extracted Flat Orders"

            imageType = self.extractedType+"Images"
            hasCalibratedImage = self.db[imageType].find_one({"_id": imageAndMaskHash[imageType]})["type"] == "Calibrated " + self.extractedType + " Image"

            return hasExtractedFlat and hasCalibratedImage
        else:
            return False










    def run(self, outputFileName=None):
        """
        Runs through the steps for the order extraction.

        Input:
            outputFileName: If None, nothing is saved. Otherwise, a string with the name of the outputfile,
                            incl. the extension ".fits".

        Output:
            xCoordinates:  x-coordinates of the pixels that belong to the orders [pix]   TODO: specify: row or column?
            yCoordinates:  y-coordinates of the pixels that belong to the orders [pix]
            fluxValues:    values of the pixels that belong to the orders        [ADU]
            orders:        number of the order to which the pixels belong
        """
        # Get the image

        imagePath = self.imageCollection.find_one({"_id" : self.imageHash})["path"]
        image = tools.getImage(imagePath)

        # Get the order mask, previously derived using a flatfield

        maskOfEachOrder = self.getMaskOfEachOrder()

        # Extract the different orders using the mask

        xCoordinates, yCoordinates, fluxValues, orders = self.extractStripes(image, maskOfEachOrder)

        # If required, save to FITS

        if outputFileName is not None:
            self.saveImageAndAddToDatabase(outputFileName, xCoordinates, yCoordinates, fluxValues, orders)
            print("Orders saved to fits file.")

        # That's it!

        return xCoordinates, yCoordinates, fluxValues, orders





    def getMaskOfEachOrder(self):
        """
        Return a dicitionary with the mask (pixel coordinates) of each order

        Input:
            None

        Output:
            mask: TODO: describe the structure of the mask dictionary
        """
        extractedFlatPath = self.db["ExtractedOrders"].find_one({"_id" : self.orderMaskHash})["path"]
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
        """
        For each order extract all pixels within the stripe of this order, and keep the x-, and y- coordinates
        of this pixel as well as its flux value.

        Input:
            image:
            mask:

        Output:
            xCoordinates:  x-coordinates of the pixels that belong to the orders [pix]   TODO: specify: row or column?
            yCoordinates:  y-coordinates of the pixels that belong to the orders [pix]
            fluxValues:    values of the pixels that belong to the orders        [ADU]
            orders:        number of the order to which the pixels belong
        """
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





    def saveImageAndAddToDatabase(self, outputFileName, xValues, yValues, flux, orders):
        """
        Save the image to a fits file and add the file to the database

        Input:
            outputFileName:
            xValues:
            yValues:
            flux:
            orders:

        Output:
            None

        """

        combinedHash = self.orderMaskHash + self.imageHash
        hash = hashlib.sha256(bytes(combinedHash, 'utf-8')).hexdigest()

        path = self.outputPath + outputFileName
        orders, fibers = zip(*orders)

        # Save Extracted Flat Orders as FITS file for primary HDU

        primary_hdr = fits.Header()
        primary_hdr["hash"] = hash
        primary_hdr["path"] = path
        primary_hdr["type"] = self.outputType
        primary_hdr["orders"] = str(set(orders))
        primary_hdr["fibers"] = str(set(fibers))
        primary_hdr["input"] = str([self.orderMaskHash, self.imageHash])


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
        currentTime = datetime.now()
        dict = {"_id"  : hash,
                "path" : path,
                "type" : self.outputType,
                "date_created" : currentTime.strftime("%d/%m/%Y %H:%M:%S")}
        tools.addToDataBase(dict, self.db, overWrite = True)





if __name__ == "__main__":

    db = DatabaseFromLocalFiles("pipelineDatabase.txt")
    print("")

    # Extract the orders of a science image

    scienceImageHash = "b0ef6a99bde7cdbc968a46fcd7a57e450a554c548d9cc89d7a9555e7236fe05f"
    scienceImagePath = "Data/ProcessedData/CalibratedScience/testFScience.fits"

    orderMaskHash    = "2133e1778dd6f8208763eeeed3a5ae6efd525fabb33a5fdf8809bd77bf89bb2b"
    orderMaskPath    = "Data/ProcessedData/ExtractedOrders/testFMask.fits"

    scienceExtractor1 = OrderExtraction(db, debug=1, ExtractedOrders=orderMaskPath,
                                        ScienceImages=scienceImagePath)
    scienceExtractor2 = OrderExtraction(debug=1, ExtractedOrders=orderMaskHash,
                                        ScienceImages=scienceImageHash)

    scienceExtractor1.run("extractedScienceTestF.fits")
    print("+============================+")
    scienceExtractor2.run("extractedScienceTestD.fits")
    print("+============================+")


    # Extract the orders of an Etalon image (Calibrated Etalon)

    etalonImageHash = "372bd5de4de92bb991f6eb3991a05c423063b67788ab1275cf31d060a347d533"
    etalonImagePath   = "Data/ProcessedData/CalibratedEtalon/testFEtalon.fits"

    orderMaskHash = "2133e1778dd6f8208763eeeed3a5ae6efd525fabb33a5fdf8809bd77bf89bb2b"
    orderMaskPath    = "Data/ProcessedData/ExtractedOrders/testFMask.fits"

    etalonExtractor1 = OrderExtraction(db, debug=1, ExtractedOrders=orderMaskPath,
                                        EtalonImages=etalonImagePath)
    etalonExtractor2 = OrderExtraction(debug=1, ExtractedOrders=orderMaskHash,
                                        EtalonImages=etalonImageHash)

    etalonExtractor1.run("extractedEtalonTestF.fits")
    print("+============================+")
    etalonExtractor2.run("extractedEtalonTestD.fits")

    db.saveToFile("pipelineDatabase.txt")

