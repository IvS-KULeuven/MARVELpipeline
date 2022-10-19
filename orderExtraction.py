import numpy as np
import os
from pipeline import PipelineComponent
import tools
from numba import njit, jit, vectorize
import matplotlib.pyplot as plt
from astropy.io import fits
from datetime import datetime
import hashlib


class OrderExtraction(PipelineComponent):
    """
    Class that performs the extraction of the orders given an image and an order mask.
    """

    def __init__(self, imageHash, extractedType, orderMaskHash, debug=0):
        """
        Initialize the order extraction component

        Input:
            imageHash:      string containing the hashes of the CCD images to process
            extractedType:  string. Type of the CCD image. Can be either "Science" or "Etalon".
            orderMaskHash:  string. The hash of the fits file containing the mask of the orders (derived from a flat)
                            and their derived flux.
            debug:          0, 1, 2, or 3:  0 meaning no debug output, 3 meaning lots of debug output.

        Output:
            None
        """

        self.extractedType = extractedType                             # Image type, e.g. "Science" or "Etalon"
        super().__init__([imageHash, orderMaskHash])                   # Reproducibility: keep track of all hashes used

        isSane = self.checkSanityOfInputHashes(imageHash, extractedType, orderMaskHash)  # will only work after super was initialized
        if not isSane:
            print("Input hashes not sane. Aborting.")                  # TODO: better error message
            exit(1)

        self.imageCollection = self.db[self.extractedType + "Images"]
        self.orderMaskCollection = self.db["ExtractedOrders"]
        self.imageHash = imageHash
        self.orderMaskHash = orderMaskHash
        self.debug = debug
        self.outputPath = os.getcwd() + "/Data/ProcessedData/ExtractedOrders/"

        if not os.path.isdir(self.outputPath):
            os.mkdir(self.outputPath)

        self.outputType = "Extracted {} Orders".format(extractedType)






    def checkSanityOfInputHashes(self, imageHash, imageType, orderMaskHash):
        """
        Check if the imageHash indeed corresponds to an image, and check if the
        orderMaskHash indeed corresponds to a fits file containing the order masks.

        Input:
            imageHash:     string containing the hashes of the images to be order-extracted
            imageType:     type of the image for which the hash is given.
                           Can be "Science" or "Etalon".
            orderMaskHash: hash of the fits file containing the order mask (obtained in a previously
                           step using a flatfield)

        Output:
            isSane: boolean. True if both input hashes are sane, False if at least one of the
                             images is unsane.
        """

        # Check if the CCD image is of the type specified

        if imageType == "Science":
            image = self.db["ScienceImages"].find_one({"_id": imageHash})
            if image["type"] != "Calibrated Science Image":
                return False
        elif imageType == "Etalon":
            image = self.db["EtalonImages"].find_one({"_id": imageHash})
            if image["type"] != "Calibrated Etalon Image":
                return False

        # Check if the orderMaskHash indeed relates to an order mask.

        image = self.db["ExtractedOrders"].find_one({"_id": orderMaskHash})
        if image is None:
            return False
        elif image["type"] != "Extracted Flat Orders":
                return False

        # If we reach this point, nothing abnormal was detected.

        return True






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
        extractedFlatPath = self.orderMaskCollection.find_one({"_id" : orderMaskHash})["path"]
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
        hash = hashlib.sha256(bytes("".join(self.input), 'utf-8')).hexdigest()
        path = self.outputPath + outputFileName
        orders, fibers = zip(*orders)

        # Save Extracted Flat Orders as FITS file for primary HDU

        primary_hdr = fits.Header()
        primary_hdr["hash"] = hash
        primary_hdr["path"] = path
        primary_hdr["type"] = self.outputType
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
        currentTime = datetime.now()
        dict = {"_id"  : hash,
                "path" : path,
                "type" : self.outputType,
                "date_created" : currentTime.strftime("%d/%m/%Y %H:%M:%S")}
        tools.addToDataBase(dict, overWrite = True)





if __name__ == "__main__":

    # Extract the orders of a science image

    scienceImageHash = "b0ef6a99bde7cdbc968a46fcd7a57e450a554c548d9cc89d7a9555e7236fe05f"
    orderMaskHash    = "2133e1778dd6f8208763eeeed3a5ae6efd525fabb33a5fdf8809bd77bf89bb2b"

    extractor = OrderExtraction(scienceImageHash, "Science", orderMaskHash, debug=3)
    outputFileName = "extracted_science_orders.fits"
    dummy = extractor.run(outputFileName)

    # Extract the orders of an Etalon image (Calibrated Etalon)

    etalonImageHashes = "372bd5de4de92bb991f6eb3991a05c423063b67788ab1275cf31d060a347d533"
    orderMaskHash = "2133e1778dd6f8208763eeeed3a5ae6efd525fabb33a5fdf8809bd77bf89bb2b"
    extractor = OrderExtraction(etalonImageHashes, "Etalon", orderMaskHash, debug=3)
    outputFileName = "extracted_etalon_orders.fits"
    dummy = extractor.run(outputFileName)

