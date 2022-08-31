
import math
import copy
import numpy as np
import pandas as pd
import h5py
from astropy.io import fits
import statsmodels.api as sm
from lmfit.models import ConstantModel, GaussianModel

from matplotlib import pyplot as plt

from tools import getAllOptimalExtractedSpectrum



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




def initialThArGaussianParameterGuesses(pathThArSpectrum):
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

    spectraAllFibersAllOrders = getAllOptimalExtractedSpectrum(pathThArSpectrum)

    # Open the file containing the theoretical distribution of the wavelength on the CCD.
    # This information is produced by PyEchelle.

    pathSimulatedWavelengthDistribution = "Data/MARVEL_2021_11_22_detector1.hdf"
    pyEchelleOutputFile = h5py.File(pathSimulatedWavelengthDistribution, "r")

    lineParameterGuess = {}

    for order in range(30, 99):

        print("Processing order #{0} to do a crude estimate of the ThAr line centers".format(order))

        # Extract the PyEchelle info for the relevant order. The ThAr fiber is fiber nr 5.
        
        fiber = 5
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
            print("No lines in NIST line list for order {0} in wavelength range [{1}, {2}] nm".format(order, wavelength.min(), wavelength.max()))
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










def wavelengthCalibrationOfOneOrder(xpix, flux, centerGuesses, amplitudeGuesses, stdevGuesses, lambdaTrue):
    """
    Returns the fitted centers of the given ThAr lines using a non-linear fit of a constant + Gaussian                        [pix] 
   
    INPUT:
        xpix:             1D Numpy array: x- pixel coordinates of one order of the the observed ThAr spectrum                 [pix]
        flux:             1D Numpy arary: corresponding flux of the observed ThAr spectrum
        centerGuesses:    1D Numpy array: for each ThAr line to be included in the fit: first guess of its center             [pix]
        amplitudeGuesses: 1D Numpy array: for each ThAr line to be included in the fit: first guess of its height 
        stdevGuesses:     1D Numpy array: for each ThAr line to be included in the fit: first guess of its standard deviation [pix]

    OUTPUT:
        fittedGaussianCenters: 1D Numpy array, same length as centerGuesses. Resulting centers of the Gaussian 
                               fits of each of the given ThAr lines.                                                          [pix]
        polynomialCoefficients: coefficients of a polynomial that relates the x-coordinate [pix] to the wavelength [nm] 

    """

    # First fit each ThAr line individually with a constant + Gaussian 
    
    print("    - fith the individual ThAr lines")

    Ngaussians = len(centerGuesses)
    fittedGaussianCenters = []
    for n in range(Ngaussians):

        model = ConstantModel()
        params = model.make_params()
        params['c'].set(value = 0.01, min=0.0)

        gauss = GaussianModel()
        model += gauss
        params.update(gauss.make_params())

        params['center'].set(value=centerGuesses[n])
        params['amplitude'].set(value=amplitudeGuesses[n])
        params['sigma'].set(value=stdevGuesses[n])

        regionAroundLine = (xpix >= centerGuesses[n] - 8) & (xpix <= centerGuesses[n] + 8) 

        myFit = model.fit(flux[regionAroundLine], params, x=xpix[regionAroundLine])
        fittedGaussianCenters += [myFit.best_values['center']]

    print("    - determine a polynomial wavelength solution")

    fittedGaussianCenters = np.array(fittedGaussianCenters)

    # Given the fitted ThAr pixel positions and their true wavelengths, fit a polynomial to this relation
    # Use the Akaike Information Criterion (AIC) to determine the optimal degree.

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
    
    return fittedGaussianCenters, olsFit.params




 




if __name__ == "__main__":

    pathThArSpectrum = "Data/ProcessedData/OptimalExtraction/optimal_extracted_science_flux.fits"
    lineParameterGuesses = initialThArGaussianParameterGuesses(pathThArSpectrum)

    spectraAllFibersAllOrders = getAllOptimalExtractedSpectrum(pathThArSpectrum)

    fiber = 5                      # ThAr fiber
    bestFitCoefficients = {}

    for order in range(30, 99):

        print("Wavelength calibration of order {0}".format(order))

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
            print("No ThAr lines for order {0}. Skipping order.".format(order))
            continue

        fittedGaussianCenters, polynomialCoefficients = wavelengthCalibrationOfOneOrder(xobs, fluxobs, centerGuesses, amplitudeGuesses, stdevGuesses, lambdaTrue)

        bestFitCoefficients[order] = polynomialCoefficients





def plotOrderFits(bestFitCoefficients):
    """

    """

    x = np.arange(10500)

    fig, ax = plt.subplots(figsize=(10,10))
    for order in bestFitCoefficients.keys():
        θ = bestFitCoefficients[order]
        X = np.vander(x, len(θ), increasing=True)
        y = X @ θ
        ax.plot(x, y)

    plt.show()
    
