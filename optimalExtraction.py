from pipeline import PipelineComponent
from pymongo import MongoClient
from numba import njit, jit, vectorize, prange
from tqdm import tqdm
from astropy.io import fits
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import hashlib
import os
import tools
import debug
import time


class OptimalExtraction(PipelineComponent):
    """
    Class that performs the optimal extraction given an extracted image and an order mask image.
    This reduces the 2D spectrum to a 1D spectrum.
    """

    def __init__(self, extractedImageHash, extractedFlatHash , extractedType, debug=0):
        """
        Initializes the optimal extraction component.

        Input:
            extractedImageHash:  string containing the hash of the image with extracted orders.
            extractedFlatHash:   string containing the hash of the fits file containing the mask
                                 of the orders (derived from a flat) and their flux values.
            extractedType:       string. Type of the extracted order images. Can be either "Science"
                                 or "Etalon".
            debug:               0, 1, 2 or 3: 0 meaning no debug output, 3 meaning lots of debug output.

        Output:
            None
        """

        self.extractedType = extractedType                                # Image type, e.g. "Science" or "Etalon"
        super().__init__([extractedImageHash, extractedFlatHash])

        isSane = self.checkSanityOfInputHashes(extractedImageHash, extractedFlatHash, extractedType)
        if not isSane:
            print("Input hashes not sane. Aborting.")
            exit(1)

        self.imageCollection = self.db["ExtractedOrders"]
        self.extractedImageHash = extractedImageHash
        self.extractedFlatHash  = extractedFlatHash
        self.debug = debug
        self.outputPath = os.getcwd() + "/Data/ProcessedData/OptimalExtraction/"

        if not os.path.isdir(self.outputPath):
            os.mkdir(self.outputPath)

        self.outputType = "Optimal Extracted {}".format(extractedType)






    def checkSanityOfInputHashes(self, extractedImageHash, extractedFlatHash, extractedType):
        """
        Check if the extractedImageHash indeed corresponds to an image of the correct format (given by
        extractedType) and that extractedFlatHash corresponds to a fits file containing the order masks
        and their flux values.

        Input:
            extractedImageHash:   string containing the hashe of the images that we want to optimaly extract.
            extractedFlatHash:    string containing the hash of the fits file containing the mask of the
                                  orders (derived from a flat) and their flux values.
            extractedType:        string. Type of the extracted order images. Can be either "Science" or "Etalon".
            debug:                0, 1, 2 or 3: 0 meaning no debug output, 3 meaning lots of debug output.

        Output:
            isSane: boolean. True if both input hashes correspond to images with the correct type, False if at
                             least one of these hashes does not.
        """

        # Check if the extraded image hash is of the type specified.

        if extractedType == "Science":
            image = self.db["ExtractedOrders"].find_one({"_id": extractedImageHash})
            if image is None:
                return False
            elif image["type"] != "Extracted Science Orders":
                return False
        elif extractedType == "Etalon":
            image = self.db["ExtractedOrders"].find_one({"_id": extractedImageHash})
            if image is None:
                return False
            elif image["type"] != "Extracted Etalon Orders":
                return False

        # check if the extractedFlatHash indeed corresponds to an order mask.

        image = self.db["ExtractedOrders"].find_one({"_id": extractedFlatHash})
        if image is None:
            return False
        elif image["type"] != "Extracted Flat Orders":
            return False

        # If we reach this point, nothing abnormal was detected.

        return True







    def run(self, outputFileName=None):
        """
        Runs through the steps for optimal extraction.

        Input:
            outputFileName: If None, nothing is saved. Otherwise, a string with the name of the outputfile,
                            incl. the extension ".fits".

        Output:
            xCoordinates:   array of x-coordinates of the pixels for every order [pix]
            yCoordinates:   array of y-coordinates of the pixels for every order [pix]
            fluxValues:     array of flux values of the pixels for every order [ADU]
            orders:         array with the fiber/order value
        """

        # Get the image paths

        imagePath = self.imageCollection.find_one({"_id": self.extractedImageHash})["path"]
        flatPath  = self.imageCollection.find_one({"_id": self.extractedFlatHash})["path"]
        #image     = tools.getImage(imagePath)

        # Extract the 1D spectrum

        spectrum, xPixels, yPixels, orders = self.extractSpectrum(flatPath, imagePath)

        # If required, save to FITS
        if outputFileName is not None:
            self.saveImage(spectrum, xPixels, yPixels, orders)
            print("Optimal extracted orders saved to fits file.")

        # That's it!

        print("Block Generated!")
        return spectrum, xPixels, yPixels, orders






    def extractSpectrum(self, flatPath, imagePath):
        """
        Transforms the 2D extracted orders to a 1D optimal extracted spectrum.

        Input:
            flatPath:   path to the image containing the extracted flat mask
            imagePath:  path to the image containing the flux we want to extract.

        Ouptut:
            optimalSpectra:  array containing the optimal extracted spectrum for every order [ADU]
            xPixels:         array containing the x_coordinates for every order              [pix]
            yPixels:         array containing the y_coordinates for every order              [pix]
            fibersAndOrders: zip object containing the fibers and orders
        """

        # Obtain the fiber/orders that are pressent in the extracted orders image

        fibers, orders = tools.getFibersAndOrders(flatPath)

        # Obtain the fluxes and coordinates information from the extracted images

        flatSpectrum  = tools.getAllExtractedSpectrum(flatPath)
        imageSpectrum = tools.getAllExtractedSpectrum(imagePath)

        # Parse the information and store the values into np.arrays.

        xValues, yValues, fValuesF, fValuesI, ordrArray, fbrArray = self.convertToArray(imageSpectrum,
                                                                                flatSpectrum, fibers, orders)

        print("start opimal extraction")
        start = time.time()

        # Apply the optimal extraction alghoritm

        optimalSpectra, imageSpectra, flatSpectra, yPixels, xPixels = getOptimalSpectrum(xValues, yValues,
                                                                                         fValuesF, fValuesI)

        end   = time.time()
        print("\tTime: ", end-start, "s")

        fibersAndOrders = zip(ordrArray, fbrArray)
        imageSpectra, optimalSpectra, flatSpectra, yPixels, xPixels = self.stripUnphysicalValues(imageSpectra,
                                                            optimalSpectra, flatSpectra, yPixels, xPixels)

        if self.debug > 2:
            debug.plotOrdersWithSlider(flatSpectra, xValues=xPixels, yMax=11800)
            debug.plotOrdersWithSlider(imageSpectra, xValues=xPixels, yMax=7000)
            debug.plotOrdersWithSlider(optimalSpectra, xValues=xPixels, yMax=2)

        return optimalSpectra, xPixels, yPixels, fibersAndOrders







    def stripUnphysicalValues(self, image, optimal, flat, yPixels, xPixels):
        """
        Strips the unphysical values that are in the input arrays. For image, optimal and flat these values
        are np.nan and for yPixels and xPixels these are -1.
        """

        stripFloats = lambda x : [ row[~np.isnan(row)] for row in x ]
        stripInt    = lambda x : [ row[row != -1] for row in x ]

        return stripFloats(image), stripFloats(optimal), stripFloats(flat), stripInt(yPixels), stripInt(xPixels)







    def convertToArray(self, imageSpectrum, flatSpectrum, fibers, orders):
        """
        This methods converts the lists that are given as input into a numpy array.
        This is done so that numba can be used to optimize the optimal extraction.
        """

        start = time.time()
        print("start convert to array")

        max_row_length = 0

        # Let us first get the maximum length of the rows

        for o in orders:
            for f in fibers:
                x, y, flux = (flatSpectrum[o])[f]
                if len(x) > max_row_length:
                    max_row_length = len(x)

        # Initialize the np.arrays used to store the data

        outputDimension = (len(fibers)*len(orders), max_row_length)
        x_values = np.full(outputDimension, -1, dtype=np.int16)
        y_values = np.full(outputDimension, -1, dtype=np.int16)

        f_values = np.full(outputDimension, np.nan, np.float64)
        o_values = np.full(outputDimension, np.nan, np.float64)

        ordersArray = np.empty(len(fibers)*len(orders))
        fibersArray = np.empty(len(fibers)*len(orders))

        # For overy fiber/order add information to arrays

        line = 0
        for o in orders:
            for f in fibers:
                x1, y1, imageFlux = (imageSpectrum[o])[f]
                x2, y2, flatFlux = (flatSpectrum[o])[f]

                if not ((x1 == x2).all() and (y1 == y2).all()):
                    print("Flat mask does not correspond with other mask")
                    exit

                (x_values[line])[:len(x1)] = np.array(x1)
                (y_values[line])[:len(y1)] = np.array(y1)

                (f_values[line])[:len(flatFlux)]  = np.array(flatFlux)
                (o_values[line])[:len(imageFlux)] = np.array(imageFlux)

                ordersArray[line] = o
                fibersArray[line] = f

                line+=1

        end = time.time()
        print("\tTime: ", end-start, "s")
        return x_values, y_values, f_values, o_values, ordersArray, fibersArray








    def saveImage(self, spectrum, xPixels, yPixels):
        """
        Save the image and add it to the database

        Input:
            spectrum: optimal extracted flux for every fiber/order
            xPixels:  xPixels for every fiber/order
            yPixels:  yPixels for every fiber/order
        """

        # output hash gets generated and set path for output file

        hash = hashlib.sha256(bytes("".join(self.input), 'utf-8' )).hexdigest()
        path = self.outputPath + self.getFileName()

        # Save the included fiber/orders

        fibers = np.unique(self.fibers)
        orders = np.unique(self.orders)

        # Make header of fits file

        primary_hdr = fits.Header()
        primary_hdr["hash"] = hash
        primary_hdr["path"] = path
        primary_hdr["type"] = self.type
        primary_hdr["orders"] = str(set(orders))
        primary_hdr["fibers"] = str(set(fibers))
        primary_hdr["input"] = str(self.input)

        hdu = fits.PrimaryHDU(header=primary_hdr)
        hdul = fits.HDUList([hdu])

        # Add the flux and x,y coordinates 

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

        currentTime = datetime.now()
        dict = {"_id" : hash, "path" : path, "type" : self.type, "date_created" : currentTime.strftime("%d/%m/%Y %H:%M:%S")}
        tools.addToDataBase(dict, overWrite = True)









@njit()
def getSpectrum(imageFlux, flatFlux, xPositions, yPositions, readout=2000):
    """
    Applies the optimal extraction routine on a certain fiber/order

    Input:
        imageFlux:  flux values image we want to reduce [ADU]
        flatFlux:   flux values of extracted flux image [ADU]
        xPositions: x-positions of the fiber/order      [pix]
        yPositions: y-positions of the fiber/order      [pix]
        readout:    optional parameters. Readout of the CCD.

    Output:
        optimalFluxes:  optimal extracted 1D spectrum                         [ADU]
        avgImageFluxes: average (over constant x) flux values of imageFlux    [ADU]
        avgFlatFluxes:  average (over constant x) flux values of flatFlux     [ADU]
        avgYPosition:   average (over constant x) flux values of y coordinate [pix]

    Note:
        The variables avgImageFluxes and avgFlatFluxes are not needed in the algorithm.
        These values are only returned for debugging purposes.
    """

    # Initialize the arrays that will be added to the output array

    avgImageFluxes = np.full(10560, np.nan, dtype=np.float32)
    avgFlatFluxes  = np.full(10560, np.nan, dtype=np.float32)
    optimalFluxes  = np.full(10560, np.nan, dtype=np.float32)
    avgYPosition   = np.full(10560, -1, dtype=np.float32)

    # Create the optimal order for every unique x value

    for i,x in enumerate(np.unique(xPositions)):

        # when x = -1 the values at this index can be ignored
        if x == -1:
            continue

        # Select the pixels that correspond to the x value

        mask = (xPositions == x)

        image  = (imageFlux)[mask]
        flat   = (flatFlux)[mask]
        ycoord = (yPositions)[mask]

        w = 1/(readout*np.ones_like(image) + image)

        # Calculate the average image and flat fluxes.

        avgImageFluxes[i] = np.sum(image)/len(image)
        avgFlatFluxes[i]  = np.sum(flat)/len(flat)

        # Calculate the average (weigthed) y coordinate value

        if np.sum(image) == 0:
            avgYPosition[i] = np.sum(ycoord) / len(ycoord)
        else:
            avgYPosition[i] = np.sum(ycoord * image) / np.sum(image)

        # Calculate the optimal extracted spectrum

        if not np.sum(w * flat * flat)  == 0:
            optimalFluxes[i] = np.sum(w * flat * image) / np.sum(w * flat * flat)
        else:
            optimalFluxes[i] = 1

    return optimalFluxes, avgImageFluxes, avgFlatFluxes, avgYPosition





@njit(parallel=True)
def getOptimalSpectrum(xValues, yValues, flatFlux, imageFlux, readout=2000):
    """
    Runs the optimal extraction algorithm. This functions uses numba to increase the speed
    and handel parallelization.

    Input:
        xValues:   array with the x-coordinates for every order                  [pix]
        yValues:   array with the y-coordinates for every order                  [pix]
        flatFlux:  array with the flux values of the flat field for every order  [ADU]
        imageFlux: array with the flux values of the image field for every order [ADU]
        readout:   optional parameters. Readout of the CCD.

    Output:
        optimalSpectrum: array with the optimal extracted spectrum for every fiber/order [ADU]
        averageImages:   array with the average (over x-coordinates) flux of the image   [ADU]
        averageFlats:    array with the average (over x-coordinates) flux of the flat    [ADU]
        averageXPixels:  array with the average y-positions for every fiber/order        [pix]
        averageYPixels:  array with the x-positions for every fiber/order                [pix]
    """

    # Initialize the output arrays

    amountOfLines   = len(imageFlux)
    averageImages    = np.empty((amountOfLines, 10560))
    averageYPixels  = np.empty((amountOfLines, 10560))
    averageXPixels  = np.empty((amountOfLines, 10560))
    averageFlats    = np.empty((amountOfLines, 10560))
    optimalSpectrum = np.empty((amountOfLines, 10560))

    # Loop over every order/fiber

    for line in prange(amountOfLines):

        # Calculate the optimal extracted spectrum for one line

        optimalFlux, averageFlux, averageFlat, averageYPosition = getSpectrum(imageFlux[line],
                                                flatFlux[line], xValues[line], yValues[line])
        xPositions = np.unique(xValues[line])

        # Add the optimal spectrum, avg flux spectrum and average y-value to the output arrays

        averageImages[line]   = averageFlux
        averageFlats[line]    = averageFlat
        averageYPixels[line]  = averageYPosition
        optimalSpectrum[line] = optimalFlux

        averageXPixels[line].fill(-1)
        (averageXPixels[line])[:len(xPositions)] = xPositions

    return optimalSpectrum, averageImages, averageFlats, averageYPixels, averageXPixels







if __name__ == "__main__":

    # Optimal extract the orders of an extracted science image
    extractedScienceHash = "b6e9f312efa62bcb90d66b5fee19eae4a6f930e0b3999650fa9c944b5075a4da"
    extractedFlatHash  = "2133e1778dd6f8208763eeeed3a5ae6efd525fabb33a5fdf8809bd77bf89bb2b"


    optimalScience = OptimalExtraction(extractedScienceHash, extractedFlatHash, "Science", debug=3)
    dummy = optimalScience.run()

    # Optimal extract the orders of an extracted etalon image
    extractedEtalonHash = "92fe37cba12963e4ec23a44bd4fb312cb59b2d1c50933dac479fa2181bc333b2"

    optimalEtalon = OptimalExtraction(extractedEtalonHash, extractedFlatHash, "Etalon", debug=3)
    dummy = optimalEtalon.run()
