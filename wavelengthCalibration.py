import math
import copy
import numpy as np
import pandas as pd
import h5py
import os
import hashlib
from astropy.io import fits
import statsmodels.api as sm
import lmfit

from datetime          import datetime
from numba             import njit, jit
from matplotlib        import pyplot as plt
from matplotlib.colors import Normalize
from pipeline          import PipelineComponent

import matplotlib.colors as colors


import tools


plt.rc('font',   size=14)          # controls default text sizes
plt.rc('axes',   titlesize=14)     # fontsize of the axes title
plt.rc('axes',   labelsize=14)     # fontsize of the x and y labels
plt.rc('xtick',  labelsize=14)     # fontsize of the tick labels
plt.rc('ytick',  labelsize=14)     # fontsize of the tick labels
plt.rc('legend', fontsize=14)      # legend fontsize
plt.rc('figure', titlesize=14)     # fontsize of the figure title






class WavelengthCalibration(PipelineComponent):
    """
    Class that prerforms the wavelength calibration

    TODO: For new we assume that the 1th fiber image has etalon and 2-5 are science spectra.
          this should be read in from the input file.
    """

    def __init__(self, database=None, debug=0, **optimalScienceHash):
        """
        Initialize the wavelength calibration component

        Input:
            database:      if not None, DatabaseFromLocalFile object to be used as database.
                           else MongoDB is used as database.

            debug:         0, 1, 2, or 3:  0 meaning no debug output, 3 meaning lots of debug output.

            imageAndMaskHash: hashes of images to be used in the wavelength calibration
                              Given as ImageType/ImageHash (keyword/argument)

        Output:
            None
        """

        super().__init__(database, **optimalScienceHash)
        optimalScienceHash = self.inputHashes

        if self.checkSanityOfInputTypes(**optimalScienceHash):
            self.outputPath                  = os.getcwd() + "/Data/ProcessedData/WaveCalibration/"
            self.outputType                  = "Wavelength Calibrated"
            self.optimalExtractedScienceHash = optimalScienceHash["OptimalExtracted"]
            self.debug                       = debug
        else:
            raise Exception("Error: The input hashes do not match the correct type: Aborting")
            exit(1)

        if not os.path.isdir(self.outputPath):
            os.mkdir(self.outputPath)





    def checkSanityOfInputTypes(self, **imageAndMaskHash):
        """
        This function is ran after we run checkSanityOfInputHashes. This function checks that the
        input type that is given is able to generate a wavelength calibration file.
        """


        types  = list(imageAndMaskHash.keys())
        values = list(imageAndMaskHash.values())

        keyIsCorrect = len(types) == 1 and types[0] == "OptimalExtracted"

        if keyIsCorrect:
            isOptimalExtractedScience = self.db["OptimalExtracted"].find_one(
                {"_id": values[0]})["type"] == "Optimal Extracted Science"
            return isOptimalExtractedScience
        else:
            return False






    def run(self, outputFileName=None):
        """
        Runs through the steps for the wavelength calibration

        Input:
            outputFileName: If None, nothing is saved. Otherwise, a string with the name of the outputfile,
                            incl. the extension ".fits".

        Output:
            linecenters_even:
            velocities_even:
            linecenters_odd:
            velocities_odd:
        """

        pathThArSpectrum = self.db["OptimalExtracted"].find_one({"_id": self.optimalExtractedScienceHash})["path"]
        lineParameterGuesses = self.initialThArGaussianParameterGuesses(pathThArSpectrum)


        spectraAllFibersAllOrders = tools.getAllOptimalExtractedSpectrum(pathThArSpectrum)
        fibers, orders = tools.getFibersAndOrders(pathThArSpectrum)
        fiber = 1


        bestFitCoefficients = {}
        fittedGaussianCenters = []
        fittedGaussianStdevs = []
        fittedGaussianAmplitudes = []
        fittedGaussianOffsets = []
        chiSquares = []

        for order in orders:

            if self.debug > 1: print("Wavelength calibration of order {0}".format(order))

            if order in spectraAllFibersAllOrders:
                xobs, yobs, fluxobs = spectraAllFibersAllOrders[order][fiber]         # x is in [pix]
            else:
                print("Order {0} not detected in observed spectrum. Skipping.".format(order))
                continue

            if order in lineParameterGuesses.keys():
                centerGuesses = lineParameterGuesses[order]['center']
                amplitudeGuesses = lineParameterGuesses[order]['amplitude']
                stdevGuesses = 3.0 * np.ones_like(centerGuesses)
                lambdaTrue = lineParameterGuesses[order]['lambdatrue']
            else:
                if self.debug > 1: print("No ThAr lines for order {0}. Skipping order.".format(order))
                continue

            centers, stdevs, amplitudes, offsets, chiSquare, polynomialCoefficients =  \
                crudeWavelengthCalibrationOfOneOrder(xobs, fluxobs, centerGuesses, amplitudeGuesses,
                                                     stdevGuesses, lambdaTrue, self.debug)

            fittedGaussianCenters.append(centers)
            fittedGaussianStdevs.append(stdevs)
            fittedGaussianAmplitudes.append(amplitudes)
            fittedGaussianOffsets.append(offsets)
            chiSquares.append(chiSquare)


            bestFitCoefficients[order] = polynomialCoefficients



        wavelengths = (compute_wavelengths(bestFitCoefficients))

        if self.debug > 2: plot_wavelength_calibrated_spectra(wavelengths)

        if outputFileName is not None:
            self.saveImageAndAddToDatabase(outputFileName, wavelengths, fittedGaussianCenters,
                                           fittedGaussianStdevs, fittedGaussianAmplitudes,
                                           fittedGaussianOffsets, chiSquares, bestFitCoefficients)
            print("Wavelength Calibration saved to fits file.")

        # That's is!

        print("Block Generated!")
        return wavelengths, fittedGaussianCenters, fittedGaussianStdevs, fittedGaussianAmplitudes, \
            fittedGaussianOffsets, chiSquares, bestFitCoefficients







    def initialThArGaussianParameterGuesses(self, pathThArSpectrum):
        """
        Get estimates for the pixel position and the amplitude of the Gaussians used to model ThAr lines.
        This function assumes that the actual distribution of ThAr lines on the CCD does not differ
        significantly from the simulated one of PyEchelle.

        Input:
            pathThArSpectrum: absolute path of the FITS file containing the observed optimally extracted MARVEL spectrum [str]

        Output:
            lineParameterGuess: dictionary. Keys are the order numbers, values are again a dictionary with the following keys:
            center:     numpy array with the guessed centers of the ThAr line          [pix]
            amplitude:  numpy array with the guessed amplitude of the ThAr line        [flux]
            lambdatrue: numpy array with the true (NIST) wavelengths of the ThAr line  [nm]

        Example: lineParameterGuess[90]['amplitude'] contains the estimated amplitudes of the
                         gaussians to be used to model the ThAr lines in order 90, for a wavelength
                         calibration.
        """

        # Get the ThAr line list that we'll be using
        # NIST recommends to use the Ritz wavelength. We skip the few lines that don't have a Ritz wavelength
        # mentioned in the file.

        pathThArLineList = "Data/ThArLineLists/ThAr_linelist_NIST_Nave2017.csv"
        lineList = pd.read_csv(pathThArLineList)
        ritzwavelength = lineList['Ritz_wavelength'].values
        selection = (ritzwavelength != '-')
        ritzwavelength = np.array(ritzwavelength[selection], dtype=np.float64)     # Sometimes there are empty entries in the NIST line list

        # Reads the observed spectrum of a whole CCD, all fibers and all orders

        spectraAllFibersAllOrders = tools.getAllOptimalExtractedSpectrum(pathThArSpectrum)

        # Open the file containing the theoretical distribution of the wavelength on the CCD.
        # This information is produced by PyEchelle.

        pathSimulatedWavelengthDistribution = "Data/MARVEL_2021_11_22_detector1.hdf"
        pyEchelleOutputFile = h5py.File(pathSimulatedWavelengthDistribution, "r")

        lineParameterGuess = {}
        _, orders = tools.getFibersAndOrders(pathThArSpectrum)
        for order in orders:

            if self.debug > 1: print("Processing order #{0}: initial guess of the ThAr line centers".format(order))

            # Extract the PyEchelle info for the relevant order. The ThAr fiber is fiber nr 1.

            fiber = 1
            dataset = pyEchelleOutputFile["wavelengths"]["fiber_{0:1d}".format(fiber)]["order_{0:2d}".format(order)]

            wavelength = np.array(dataset) * 1000.0                       # [nm]
            xsimul = np.arange(len(wavelength))                           # [pix]

            # Extract the observed ThAr spectrum for the relevant order.

            if order in spectraAllFibersAllOrders:
                xobs, yobs, fluxobs = spectraAllFibersAllOrders[order][fiber]    # x is in [pix]
            else:
                print("Order {0} not detected in observed spectrum. Skipping.".format(order))
                continue

            # Select only the lines in the line list relevant to the current order

            selection = (ritzwavelength > wavelength.min()) & (ritzwavelength < wavelength.max())
            if len(ritzwavelength[selection]) == 0:
                if self.debug > 1: print("No lines in NIST line list for order {0} in wavelength range [{1}, {2}] nm".format(order, wavelength.min(), wavelength.max()))
                continue

            # For each ThAr line in the NIST line lines, find the closest wavelength in the PyEchelle simulated
            # wavelength distribution and select the corresponding pixel x-coordinate. Given the latter, find
            # the closest pixel x-coordinate in the observed ThAr spectrum. The result is a link between the
            # wavelength of the NIST ThAr line and the pixel x-coordinate where it's supposed to be in the observed
            # ThAr spectrum

            index = []
            lambdaTrue = []
            for lamda in ritzwavelength[selection]:
                idx = find_index_nearest(wavelength, lamda)
                nearest_x = xsimul[idx]
                candidate_index = find_index_nearest(xobs, nearest_x)
                if candidate_index > 0 and candidate_index < len(xobs):
                    index += [candidate_index]
                    lambdaTrue += [lamda]

            # Remove all elements in 'index' that occur more than once.

            unique = np.array([element for element in index if not index.count(element) > 1])
            lambdaTrue = np.array([lambdaTrue[n] for n in range(len(index)) if not index.count(index[n]) > 1])
            lineCenterGuess = xobs[unique]                                                             # [pix]
            lineAmplitudeGuess = fluxobs[unique]

            # Unfortunately, the NIST line list is not always good, and we need to clean it.
            # Sometimes NIST gives 2 or more lines where there is only one. This will give problems
            # when we try to fit a wavelength solution. I tried using the relative_intensity, but
            # this did not give good results.
            # The cleaning below has the consequence that we eliminate some clearly blended lines,
            # but also removes a few good isolated lines for which NIST erroneously thinks there
            # are two or more very close lines.

            # Only retain those lines who do not have another neighbouring line within 6 pixels,
            # both on the left as well as on the right side.

            isolated_left = np.concatenate([[True], (lineCenterGuess[1:] - lineCenterGuess[:-1]) >= 6])
            isolated_right = np.concatenate([np.abs(lineCenterGuess[:-1] - lineCenterGuess[1:]) >= 6, [True]])
            isolated = isolated_left & isolated_right
            lineCenterGuess = lineCenterGuess[isolated]
            lineAmplitudeGuess = lineAmplitudeGuess[isolated]
            lambdaTrue = lambdaTrue[isolated]

            lineParameterGuess[order] = {'center': lineCenterGuess, 'amplitude': lineAmplitudeGuess, 'lambdatrue': lambdaTrue}

        # That's it!

        return lineParameterGuess




    def plot_assessment(self, unique_orders, linecenters_even, velocities_even, linecenters_odd, velocities_odd):

        cmap = plt.cm.get_cmap('plasma')
        fig, ax = plt.subplots(figsize=(8,12))
        orders = []
        linecenters = []
        velocities = []
        for order in unique_orders:
            orders      += [order] * len(linecenters_even[order])
            linecenters += list(linecenters_even[order])
            velocities  += list(velocities_even[order])
            orders      += [order] * len(linecenters_odd[order])
            linecenters += list(linecenters_odd[order])
            velocities  += list(velocities_odd[order])

        vmin, vmax = min(velocities), max(velocities),
        scatterplot = ax.scatter(linecenters, orders, s=8, c=velocities, vmin=vmin, vmax=vmax, cmap=cmap)
        cb = plt.colorbar(scatterplot, ax=ax, location='bottom', pad=0.08, shrink=0.9, label="velocity (m/s)")
        ax.set(xlabel="Along-dispersion coordinate [pix]", ylabel="order")
        plt.show()






    def saveImageAndAddToDatabase(self, outputFileName, wavelength, fittedGaussianCenters, fittedGaussianStdevs,
                                  fittedGaussianAmplitudes, fittedGaussianOffsets,
                                  chiSquares, bestFitCoefficients):
        """
        Save the image and add it to the database

        Input:
            outputFileName : string with the name of the file
            spectrum       : optimal extracted flux for every fiber/order
            xPixels        : xPixels for every fiber/order
            yPixels        : yPixels for every fiber/order
            orders         : fiber/orders

        TODO: Add flux in there as well
        """

        hash = hashlib.sha256(bytes(self.optimalExtractedScienceHash, 'utf-8')).hexdigest()
        path = self.outputPath + outputFileName

        orders = list(wavelength.keys())
        fibers = list(wavelength[list(orders)[0]].keys())

        # Save wavelength calibration in the FITS file

        primary_hdr = fits.Header()
        primary_hdr["hash"] = hash
        primary_hdr["path"] = path
        primary_hdr["type"] = self.outputType
        primary_hdr["orders"] = str(set(np.unique(orders)))
        primary_hdr["fibers"] = str(set(np.unique(fibers)))
        primary_hdr["input"] = str([self.optimalExtractedScienceHash])

        hdu = fits.PrimaryHDU(header=primary_hdr)
        hdul = fits.HDUList([hdu])

        for i, order in enumerate(orders):
            hdr1 = fits.Header()
            hdr1["order"] = order
            for j, fiber in enumerate(fibers):

                hdr2 = fits.Header()
                hdr2["order"] = order
                hdr2["fiber"] = fiber

                wLength = np.array((wavelength[order][fiber])["lambda"], dtype=np.float64)
                col1    = fits.Column("wavelength", format='D', array=wLength)

                flux = np.array((wavelength[order][fiber])["flux"], dtype=np.float64)
                col2 = fits.Column("flux", format='D', array=flux)

                cols2 = fits.ColDefs([col1, col2])
                hdu2 = fits.BinTableHDU.from_columns(cols2, header=hdr2)

                hdul.append(hdu2)

            bestFit = np.array(bestFitCoefficients[order], dtype=np.float64)
            col3 = fits.Column("best fit coefficients", format='D', array=bestFit)

            gaussianCenter = np.array(fittedGaussianCenters[i], dtype=np.float64)
            col4 = fits.Column("fitted Gaussian centers", format='D', array=gaussianCenter)

            gaussianStdevs = np.array(fittedGaussianStdevs[i], dtype=np.float64)
            col5 = fits.Column("gaussian Stdevs", format='D', array=gaussianStdevs)

            gaussianAmplitudes = np.array(fittedGaussianAmplitudes[i], dtype=np.float64)
            col6 = fits.Column("gaussian amplitude", format='D', array=gaussianAmplitudes)

            gaussianOffsets = np.array(fittedGaussianOffsets[i], dtype=np.float64)
            col7 = fits.Column("gaussian offset", format='D', array=gaussianOffsets)

            chiSquare = np.array(chiSquares[i], dtype=np.float64)
            col8 = fits.Column("chi square", format='D', array=chiSquare)

            col1 = fits.ColDefs([col3, col4, col5, col6, col7, col8])
            hdu1 = fits.BinTableHDU.from_columns(col1, header=hdr1)

            hdul.append(hdu1)

            hdul.writeto(path, overwrite=True)

            # Add image to the database
            currentTime = datetime.now()
            dict = {"_id"  : hash,
                    "path" : path,
                    "type" : self.outputType,
                    "date_created" : currentTime.strftime("%d/%m/%Y %H:%M:%S")}
            tools.addToDataBase(dict, self.db, overWrite = True)












def compute_wavelengths(polynomialCoefficients):
    """
    Compute the wavelengths corresponding to each x-coord pixel value of 1D science spectrum

    INPUT: a dictionary P so that
           P[order] contains a numpy array with the polynomial coefficients of the wavelength calibration

    OUTPUT: a dictionary D so that
            D[order][fiber]["lambda"] contains the wavelengths in [nm]
            D[order][fiber]["flux"] contains the flux values

            Here, fiber runs from 1 to 4 (incl) and order from 30 to 99.
    """

    scienceSpectrum = "Data/ProcessedData/OptimalExtraction/optimal_extracted_science_flux.fits"
    spectraAllFibersAllOrders = tools.getAllOptimalExtractedSpectrum(scienceSpectrum)

    wavelength_calibrated_spectra = {}

    for order in polynomialCoefficients.keys():
        wavelength_calibrated_spectra[order] = {2: {}, 3: {}, 4: {}, 5: {}}
        for fiber in [2,3,4, 5]:
            xobs, yobs, fluxobs = spectraAllFibersAllOrders[order][fiber]                   # xobs is in [pix]
            X = np.vander(xobs, len(polynomialCoefficients[order]), increasing=True)        # Design matrix
            wavelength = X @ polynomialCoefficients[order]                                  # [nm]
            wavelength_calibrated_spectra[order][fiber] = {"lambda": wavelength, "flux": fluxobs}

    # That's it!

    return wavelength_calibrated_spectra




def plot_wavelength_calibrated_spectra(wavelength_calibrated_spectra):

    fig, ax = plt.subplots(figsize=(9,12))
    # for order in wavelength_calibrated_spectra.keys():
    #     for fiber in [2,3,4, 5]:
    #         wavelength = wavelength_calibrated_spectra[order][fiber]['lambda']
    #         flux = wavelength_calibrated_spectra[order][fiber]['flux']
    #         ax.scatter(wavelength, flux, label=f"fiber: {fiber}, order: {order}", s=2)

    order = 80
    for fiber in [2,3,4, 5]:
        wavelength = wavelength_calibrated_spectra[order][fiber]['lambda']
        flux = wavelength_calibrated_spectra[order][fiber]['flux']
        ax.scatter(wavelength, flux, label=f"{fiber}, {order}", s=2)
    plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    plt.show()
































































@njit()
def find_index_nearest(array, value):
    """
    Returns the index of the element in 'array' that is closest to the given 'value'.
    This function **assumes** that array is sorted in ascending order!
    """
    idx = np.searchsorted(array, value, side="right")         # right means array[idx-1] <= value < array[idx]
    if (idx <= 0) or (idx == len(array)): return idx
    left_difference = math.fabs(value - array[idx-1])
    right_difference = math.fabs(value - array[idx])
    if left_difference <= right_difference:
        return idx-1
    else:
        return idx







@njit()
def erf(x):
    """
    Calculates the erf function

    Input:
        x. np.array with the x coordinates

    Output:
        y. np.array with evalutes values
    """

    # save the sign of x
    sign = np.ones_like(x)
    sign[x < 0] = -1

    x = np.absolute(x)

    # constants
    a1 =  0.254829592
    a2 = -0.284496736
    a3 =  1.421413741
    a4 = -1.453152027
    a5 =  1.061405429
    p  =  0.3275911

    # A&S formula 7.1.26
    t = np.ones_like(x) / (1.0 + p*x)
    y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*np.exp(-x*x)

    return sign*y # erf(-x) = -erf(x)







@njit()
def discretizedGaussian(x, offset, amplitude, mu, sigma):
    """
    Compute the value
    \int_a^b (c + A e^{-(x-μ)^2/(2σ^2)}) dx
         = c (b-a) + A σ \sqrt{π/2} [erf((b-μ)/sqrt{2}/σ) - erf((a-μ)/sqrt{2}/σ)]}


    Input:
        x:        1D Numpy array: along-dispersion pixel coordinates [pix]
        offset:   offset of the gaussian to be fitted.
        mu:       mean of the gaussian
        sigma:    standard deviation of the gaussian

    Output:
        flux:     1D Numpy array of length equal to 'x'. y-values of the gaussian.
    """

    a = np.floor(x)                # Lower boundary of the pixel
    b = np.ceil(x)                 # Upper boundary of the pixel
    b[b==a] += 1                   # a == b if x is exactly on a pixel boundary (i.e. an integer number)

    return offset * (b-a)     \
           + amplitude * sigma * np.sqrt(np.pi / 2)  * (erf((b - mu)/(sigma * np.sqrt(2.))) - erf((a - mu)/(sigma * np.sqrt(2.))))





def fitThArLines(xpix, flux, centerGuesses, amplitudeGuesses, stdevGuesses, lambdaTrue):

    """
    Returns the fitted centers [pix] of the given ThAr lines using a non-linear fit of a constant + Gaussian
    as well as the coefficients of the polynomial relating the pixel coordinates [pix] to the wavelengths [nm]

    INPUT:
        xpix:             1D Numpy array: x-pixel coordinates of one order of the the observed ThAr spectrum                  [pix]
        flux:             1D Numpy arary: corresponding flux of the observed ThAr spectrum
        centerGuesses:    1D Numpy array: for each ThAr line to be included in the fit: first guess of its center             [pix]
        amplitudeGuesses: 1D Numpy array: for each ThAr line to be included in the fit: first guess of its height
        stdevGuesses:     1D Numpy array: for each ThAr line to be included in the fit: first guess of its standard deviation [pix]
        lambdaTrue:       1D Numpy array: Ritz wavelengths of the ThAr lines used to do a wavelength calibration              [nm]

    OUTPUT:
        fittedGaussianCenters:    1D Numpy array, same length as centerGuesses. Resulting centers of the Gaussian
                                  fits of each of the given ThAr lines.                                                          [pix]
        fittedGaussianStdevs:     Corresponding standard deviations of the Gaussian fits of each of the given ThAr lines.        [pix]
        fittedGaussianAmplitudes: Corresponding amplitudes of the Gaussian fits of each of the given ThAr lines.
        fittedGaussianOffsets:    Corresponding constant offsets of the Gaussian fits of each of the given ThAr lines.
        chiSquares:               Corresponding chi-square value of the fits for each of the given ThAr lines.
    """

    Ngaussians = len(centerGuesses)
    fittedGaussianCenters = []
    fittedGaussianStdevs = []
    fittedGaussianAmplitudes = []
    fittedGaussianOffsets = []
    chiSquares = []

    for n in range(Ngaussians):
        myModel = lmfit.Model(discretizedGaussian)
        myModel.set_param_hint('offset', value=0.01, min=0.0)
        myModel.set_param_hint('amplitude', value=amplitudeGuesses[n], min=0.0)
        myModel.set_param_hint('mu', value=centerGuesses[n], min=0.0)
        myModel.set_param_hint('sigma', value=stdevGuesses[n], min=0.0)

        regionAroundLine = (xpix >= centerGuesses[n] - 8) & (xpix <= centerGuesses[n] + 8)

        myFit = myModel.fit(flux[regionAroundLine], x=xpix[regionAroundLine])
        fittedGaussianCenters    += [myFit.best_values['mu']]
        fittedGaussianStdevs     += [myFit.best_values['sigma']]
        fittedGaussianAmplitudes += [myFit.best_values['amplitude']]
        fittedGaussianOffsets    += [myFit.best_values['offset']]
        chiSquares += [myFit.chisqr]

    return np.array(fittedGaussianCenters), np.array(fittedGaussianStdevs), np.array(fittedGaussianAmplitudes), \
           np.array(fittedGaussianOffsets), np.array(chiSquares)



def crudeWavelengthCalibrationOfOneOrder(xpix, flux, centerGuesses, amplitudeGuesses, stdevGuesses, lambdaTrue, debug):
    """
    Returns the fitted centers [pix] of the given ThAr lines using a non-linear fit of a constant + Gaussian
    as well as the coefficients of the polynomial relating the pixel coordinates [pix] to the wavelengths [nm]

    INPUT:
        xpix:             1D Numpy array: x-pixel coordinates of one order of the the observed ThAr spectrum                  [pix]
        flux:             1D Numpy arary: corresponding flux of the observed ThAr spectrum
        centerGuesses:    1D Numpy array: for each ThAr line to be included in the fit: first guess of its center             [pix]
        amplitudeGuesses: 1D Numpy array: for each ThAr line to be included in the fit: first guess of its height
        stdevGuesses:     1D Numpy array: for each ThAr line to be included in the fit: first guess of its standard deviation [pix]
        lambdaTrue:       1D Numpy array: Ritz wavelengths of the ThAr lines used to do a wavelength calibration              [nm]
        debug:            0, 1, 2, or 3:  0 meaning no debug output, 3 meaning lots of debug output.

    OUTPUT:
        fittedGaussianCenters:  1D Numpy array, same length as centerGuesses. Resulting centers of the Gaussian
                                fits of each of the given ThAr lines.                                                          [pix]
        chiSquare:              Corresponding chi-square values of the fits for each of the given ThAr lines.
        polynomialCoefficients: Coefficients of a polynomial that relates the x-coordinate [pix] to the wavelength [nm]
                                a[n] contains the coefficient of x^n
    """

    # First fit each ThAr line individually with a constant + Gaussian

    if debug > 1: print("    - fitting the individual ThAr lines")

    fittedGaussianCenters, fittedGaussianStdevs, fittedGaussianAmplitudes, fittedGaussianOffsets, chiSquares = \
                fitThArLines(xpix, flux, centerGuesses, amplitudeGuesses, stdevGuesses, lambdaTrue)

    # Given the fitted ThAr pixel positions and their true wavelengths, fit a polynomial to this relation
    # Use the Akaike Information Criterion (AIC) to determine the optimal degree.

    if debug > 1: print("    - determine a polynomial wavelength solution")

    AIC = []
    polynomialCoefficients = []
    candidatePolynomialDegrees = [2,3,4,5,6]
    for polynomialDegree in candidatePolynomialDegrees:
        X = np.vander(fittedGaussianCenters, polynomialDegree+1, increasing=True)
        olsFit = sm.OLS(lambdaTrue, X).fit()
        AIC += [olsFit.aic]

    bestAICindex = np.argmin(AIC)
    bestPolynomialDegree = candidatePolynomialDegrees[bestAICindex]
    X = np.vander(fittedGaussianCenters, bestPolynomialDegree+1, increasing=True)
    olsFit = sm.OLS(lambdaTrue, X).fit()

    return fittedGaussianCenters, fittedGaussianStdevs, fittedGaussianAmplitudes, fittedGaussianOffsets, chiSquares, olsFit.params



































# def plot(fitttedGaussianCenters, orders, chiSquares):
#     cmap = plt.colormaps["cividis"]
#     fig, ax = plt.subplots(figsize=(8,12))
#     scatterplot = ax.scatter(fittedGaussianCenters, orders, s=8, color=cmap(chiSquares), norm=colors.LogNorm())
#     fig.colorbar(scatterplot, ax=ax, location='bottom', pad=0.08, label="chi-square")
#     ax.set(xlabel="Along-dispersion coordinate [pix]", ylabel="order")
#     plt.show()



# def plotOrderFits(bestFitCoefficients):
#     """

#     """

#     x = np.arange(10500)

#     fig, ax = plt.subplots(figsize=(10,10))
#     for order in bestFitCoefficients.keys():
#         θ = bestFitCoefficients[order]
#         X = np.vander(x, len(θ), increasing=True)
#         y = X @ θ
#         ax.plot(x, y)

#     plt.show()










if __name__ == "__main__":

    # fittedGaussianCenters, fittedGaussianStdevs, fittedGaussianAmplitudes, fittedGaussianOffsets,    \
    # chiSquares, bestFitCoefficients, orders = wavelength_calibration()

    # unique_orders, linecenters_even, velocities_even, linecenters_odd, velocities_odd = assessment()
    # plot_assessment(unique_orders, linecenters_even, velocities_even, linecenters_odd, velocities_odd)

    extracted_science_hash = "28065b0282df4865863b0f4e1aa7d7367cf7cbe2226d642cef5b5a147e3d6e43"
    extracted_science_path = "Data/ProcessedData/OptimalExtraction/optimal_extracted_science_flux.fits"
    item = WavelengthCalibration(debug=3, OptimalExtracted=extracted_science_path)
    item.run("output_wave_calibration.fits")


    hdu_list = fits.open(os.getcwd() + "/Data/ProcessedData/WaveCalibration/output_wave_calibration.fits",
                         memmap=True)

    print(hdu_list.info())

