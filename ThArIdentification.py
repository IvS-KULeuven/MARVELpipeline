
import math
import copy
import numpy as np
import pandas as pd
import h5py
from astropy.io import fits
from lmfit.models import ConstantModel, GaussianModel

from tools import getAllOptimalExtractedSpectrum
import matplotlib.pyplot as plt



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
    
    Input:
        pathThArSpectrum: absolute path of the FITS file containing the observed optimally extracted MARVEL spectrum [str]

    Output:
        lineParameterGuess: dictionary. Keys are the order numbers, values are again a dictionary with the following keys:
                            center:     numpy array with the guessed centers of the ThAr line  [pix]
                            amplitude:  numpy array with the guessed amplitude of the ThAr line [flux]
                            lambdatrue: numpy array with the true (NIST) wavelengths of the ThAr line [nm]
                            
                            Example: lineParameterGuess[90]['amplitude'] contains the estimated amplitudes of the
                                     gaussians to be used to model the ThAr lines in order 90, for a wavelength
                                     calibration.
    """

    # Get the ThAr line list that we'll be using
    # NIST recommends to use the Ritz wavelength. We skip the few lines that don't have a Ritz wavelength 
    # mentioned in the file.

    pathThArLineList = "/lhome/driess/MARVEL/MARVELpipeline/Data/ThAr_linelist_NIST_Nave2017.csv"
    lineList = pd.read_csv(pathThArLineList)
    ritzwavelength = lineList['Ritz_wavelength'].values
    selection = (ritzwavelength != '-') 
    ritzwavelength = np.array(ritzwavelength[selection], dtype=np.float64)     # Sometimes there are empty entries in the NIST line list

    # Reads the observed spectrum of a whole CCD, all fibers and all orders

    spectraAllFibersAllOrders = getAllOptimalExtractedSpectrum(pathThArSpectrum)
    
    # Open the file containing the theoretical distribution of the wavelength on the CCD.
    # This information is produced by PyEchelle.

    pathSimulatedWavelengthDistribution = "/lhome/driess/MARVEL/MARVELpipeline/Data/MARVEL_2021_11_22_detector1.hdf"
    pyEchelleOutputFile = h5py.File(pathSimulatedWavelengthDistribution, "r")

    lineParameterGuess = {}

    for order_true in range(30, 99):

        print("Processing order #{0}".format(order_true))

        # Extract the PyEchelle info for the relevant order. The ThAr fiber is fiber nr 5.
        
        fiber_true = 5
        dataset = pyEchelleOutputFile["wavelengths"]["fiber_{0:1d}".format(fiber_true)]["order_{0:2d}".format(order_true)]
        wavelength = np.array(dataset) * 1000.0                       # [nm]
        xsimul = np.arange(len(wavelength))                           # [pix]

        # Extract the observed ThAr spectrum for the relevant order.
        # Note that the fiber and order numbers stored in this file deviate from the normal ones.
          
        if order_true in spectraAllFibersAllOrders:
            xobs, yobs, fluxobs = (spectraAllFibersAllOrders[order_true])[fiber_true]    # x is in [pix]
        else:
            print("Order {0} not detected in observed spectrum. Skipping.".format(order_true))
            continue

        # Select only the lines in the line list relevant to the current order
        selection = (ritzwavelength > wavelength.min()) & (ritzwavelength < wavelength.max()) 
        if len(ritzwavelength[selection]) == 0:
            print("No lines in NIST line list for order {0} in wavelength range [{1}, {2}] nm".format(order_true, wavelength.min(), wavelength.max()))
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
        
        lineParameterGuess[order_true] = {'center': lineCenterGuess, 'amplitude': lineAmplitudeGuess, 'lambdatrue': lambdaTrue}

    # That's it!

    return lineParameterGuess









def wavelengthCalibrationOfOneOrder(xpix, flux, centerGuesses, amplitudeGuesses, stdevGuesses, lambdaTrue):
    """

    
    """


    # Set up a model of many gaussians

    model = ConstantModel()
    params = model.make_params()
    params['c'].set(value = 0.01, min=0.0)

    for n in range(len(centerGuesses)):
        prefix = "g{0}_".format(n)
        gauss = GaussianModel(prefix=prefix)
        model += gauss
        params.update(gauss.make_params())

        params[prefix+'center'].set(value=centerGuesses[n])
        params[prefix+'amplitude'].set(value=amplitudeGuesses[n])
        params[prefix+'sigma'].set(value=stdevGuesses[n])

    return model, params

 




if __name__ == "__main__":

    pathThArSpectrum = "Data/ProcessedData/OptimalExtraction/optimal_extracted_science_flux.fits"
    lineParameterGuesses = initialThArGaussianParameterGuesses(pathThArSpectrum)
    
    spectraAllFibersAllOrders = getAllOptimalExtractedSpectrum(pathThArSpectrum)

    order = 90
    centerGuesses = lineParameterGuesses[order]['center']
    amplitudeGuesses = lineParameterGuesses[order]['amplitude']
    stdevGuesses = 3.0 * np.ones_like(centerGuesses)
    lambdaTrue = lineParameterGuesses[order]['lambdatrue']

    xpix, flux = None, None
    model, params = wavelengthCalibrationOfOneOrder(xpix, flux, centerGuesses, amplitudeGuesses, stdevGuesses, lambdaTrue)

    
    spectraAllFibersAllOrders = getAllOptimalExtractedSpectrum(pathThArSpectrum)

    x, y , f = (spectraAllFibersAllOrders[90])[5]

    mask = (f == -1)
    mask = ~mask

    # plt.plot(x[mask], f[mask])
    # plt.show()


