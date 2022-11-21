import hashlib
import os
import tools
import debug
import time

import matplotlib.pyplot as plt
import numpy as np

from pipeline   import PipelineComponent
from pymongo    import MongoClient
from numba      import njit, jit, vectorize, prange
from database   import DatabaseFromLocalFile
from tqdm       import tqdm
from astropy.io import fits
from datetime   import datetime





class OptimalExtraction(PipelineComponent):
    """
    Class that performs the optimal extraction given an extracted image and an order mask image.
    This reduces the 2D spectrum to a 1D spectrum.
    """

    def __init__(self, database=None, debug=0, **optimalExtractionHash):
        """
        Initializes the optimal extraction component.

        Input:
            database:              If not None, DatabaseFromLocalFile object that is used to
                                   find the input hashes. If is None, mongoDB is used as database.
            debug:                 0, 1, 2 or 3: 0 meaning no debug output, 3 meaning lots of debug output.

            optimalExtractionHash: dictionary with key/value the imagetType/hash(or path) of images used to
                                   determine the optimal extraction. Types should be extracter flat image,
                                   or extracted science (or etalon) image.

        Output:
            None
        """

        super().__init__(database, **optimalExtractionHash)
        optimalExtractionHash = self.inputHashes
        self.extractedFlatHash  = None
        self.extractedImageHash = None

        if self.checkSanityOfInputTypes(**optimalExtractionHash):
            self.outputPath         = os.getcwd() + "/Data/ProcessedData/OptimalExtraction/"
            self.outputType         = f"Optimal Extracted {self.extractedType}"
            self.imageCollection    = self.db["ExtractedOrders"]
            self.debug              = debug
        else:
            raise Exception("Error: The input hashes do not match the correct type: Aborting")
            exit(1)







        if not os.path.isdir(self.outputPath):
            os.mkdir(self.outputPath)








    def checkSanityOfInputTypes(self, **optimalExtractionHash):
        """
        This function is ran after we run checkSanityOfInputHashes. This function checks the the
        input types that are given is able to generate an optimal extracted image.

        Input:
            optimalExtractionHash: dictionary with key/value the imagetType/hash of images used to
                                   determine the optimal extraction. Types should be extracted flat image,
                                   or extracted science (or etalon) image.

        Output:
            isSane: boolean. True if both input hashes correspond to images with the correct type, False if at
                             least one of these hashes does not.
        """

        types = list(optimalExtractionHash.keys())
        values = list(optimalExtractionHash.values())

        # Check that the keys are of the right format. For optimal order extraction there should be one
        # extracted flat image and one extracted etalon or science frame.

        isExtracted = (len(types) == 1) and ("ExtractedOrders" in types)

        # We should check that the values corresponding to Extracted Orders are a list with two hashes.

        if isExtracted and type(values[0]) == list:
            if len(values[0]) == 2:
                typesOfFiles = [(self.db["ExtractedOrders"].find_one({"_id": iHash})["type"], iHash)
                                for iHash in values[0]]

                for iType, iHash in typesOfFiles:
                    # if hash is extracted flat order
                    if iType == "Extracted Flat Orders":
                        self.extractedFlatHash = iHash
                    # if hash is extracted image order
                    elif (iType == "Extracted Science Orders" or iType == "Extracted Etalon Orders"):
                        if iType == "Extracted Science Orders":
                            self.extractedImageHash = iHash
                            self.extractedType = "Science"
                        elif iType == "Extracted Etalon Orders":
                            self.extractedImageHash = iHash
                            self.extractedType = "Science"
                    # if hash is neither extracted image or extracted flat order
                    else:
                        return False

            else:
                return False
        else:
            return False

        # If we reach this point, nothing abnormal was detected. We should check that we have found
        # one extracted image hash and one extracted flat hash

        return (self.extractedImageHash is not None) and (self.extractedFlatHash is not None)









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
            self.saveImageAndAddToDatabase(outputFileName, spectrum, xPixels, yPixels, orders)
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








    def saveImageAndAddToDatabase(self, outputFileName, spectrum, xPixels, yPixels, orders):
        """
        Save the image and add it to the database

        Input:
            outputFileName : string with the name of the file
            spectrum       : optimal extracted flux for every fiber/order
            xPixels        : xPixels for every fiber/order
            yPixels        : yPixels for every fiber/order
            orders         : fiber/orders

        """

        # output hash gets generated and set path for output file

        combinedHash = self.extractedImageHash + self.extractedFlatHash
        hash = hashlib.sha256(bytes(combinedHash, 'utf-8' )).hexdigest()

        path = self.outputPath + outputFileName
        # Save the included fiber/orders
        orders, fibers = zip(*orders)

        # Make header of fits file
        primary_hdr = fits.Header()
        primary_hdr["hash"] = hash
        primary_hdr["path"] = path
        primary_hdr["type"] = self.outputType
        primary_hdr["orders"] = str(set(np.unique(orders)))
        primary_hdr["fibers"] = str(set(np.unique(fibers)))
        primary_hdr["input"] = str([self.extractedImageHash, self.extractedFlatHash])

        hdu = fits.PrimaryHDU(header=primary_hdr)
        hdul = fits.HDUList([hdu])

        # Add the flux and x,y coordinates

        for i in np.arange(len(spectrum)):
            hdr1 = fits.Header()
            hdr1["order"] = orders[i]
            hdr1["fiber"] = fibers[i]

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
        dict = {"_id" : hash,
                "path" : path,
                "type" : self.outputType,
                "date_created" : currentTime.strftime("%d/%m/%Y %H:%M:%S")}
        tools.addToDataBase(dict, self.db, overWrite = True)









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

    db = DatabaseFromLocalFile("pipelineDatabase.txt")
    print("")

    # Optimal extract the orders of an extracted science image

    extractedScienceHash = "b6e9f312efa62bcb90d66b5fee19eae4a6f930e0b3999650fa9c944b5075a4da"
    extractedSciencePath = "Data/ProcessedData/ExtractedOrders/extractedScienceTestF.fits"

    extractedFlatHash = "2133e1778dd6f8208763eeeed3a5ae6efd525fabb33a5fdf8809bd77bf89bb2b"
    extractedFlatPath = "Data/ProcessedData/ExtractedOrders/testFMask.fits"

    optimalScience1 = OptimalExtraction(db, debug=1, ExtractedOrders=[extractedSciencePath, extractedFlatPath])
    optimalScience2 = OptimalExtraction(debug=1, ExtractedOrders=[extractedScienceHash, extractedFlatHash])

    optimalScience1.run("optimal_extracted_science_flux1.fits")
    print("+============================+")
    optimalScience2.run("optimal_extracted_science_flux.fits")
    print("+============================+")


    # Optimal extract the orders of an extracted etalon image
    extractedEtalonHash = "92fe37cba12963e4ec23a44bd4fb312cb59b2d1c50933dac479fa2181bc333b2"
    extractedEtalonPath = "Data/ProcessedData/ExtractedOrders/extractedEtalonTestF.fits"

    optimalEtalon1 = OptimalExtraction(db, debug=1, ExtractedOrders=[extractedEtalonPath, extractedFlatPath])
    optimalEtalon2 = OptimalExtraction(debug=1, ExtractedOrders=[extractedEtalonHash, extractedFlatHash])
    optimalEtalon1.run("optimal_extracted_etalon_flux1.fits")
    print("+============================+")
    optimalEtalon2.run("optimal_extracted_etalon_flux.fits")

    db.save()





