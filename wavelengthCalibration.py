
import math
import copy
import numpy as np
import pandas as pd
import h5py
from astropy.io import fits
import statsmodels.api as sm
import lmfit
from scipy.special import erf

from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.colors as colors

from tools import getAllOptimalExtractedSpectrum



plt.rc('font',   size=14)          # controls default text sizes
plt.rc('axes',   titlesize=14)     # fontsize of the axes title
plt.rc('axes',   labelsize=14)     # fontsize of the x and y labels
plt.rc('xtick',  labelsize=14)     # fontsize of the tick labels
plt.rc('ytick',  labelsize=14)     # fontsize of the tick labels
plt.rc('legend', fontsize=14)      # legend fontsize
plt.rc('figure', titlesize=14)     # fontsize of the figure title




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

        print("Processing order #{0}: initial guess of the ThAr line centers".format(order))

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




def crudeWavelengthCalibrationOfOneOrder(xpix, flux, centerGuesses, amplitudeGuesses, stdevGuesses, lambdaTrue):
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
        fittedGaussianCenters:  1D Numpy array, same length as centerGuesses. Resulting centers of the Gaussian 
                                fits of each of the given ThAr lines.                                                          [pix]
        chiSquare:              Corresponding chi-square values of the fits for each of the given ThAr lines.
        polynomialCoefficients: Coefficients of a polynomial that relates the x-coordinate [pix] to the wavelength [nm] 
                                a[n] contains the coefficient of x^n
    """

    # First fit each ThAr line individually with a constant + Gaussian 
    
    print("    - fith the individual ThAr lines")

    fittedGaussianCenters, fittedGaussianStdevs, fittedGaussianAmplitudes, fittedGaussianOffsets, chiSquares = \
                fitThArLines(xpix, flux, centerGuesses, amplitudeGuesses, stdevGuesses, lambdaTrue)

    # Given the fitted ThAr pixel positions and their true wavelengths, fit a polynomial to this relation
    # Use the Akaike Information Criterion (AIC) to determine the optimal degree.

    print("    - determine a polynomial wavelength solution")

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











    


def wavelength_calibration():
    """ 
    """

    pathThArSpectrum = "Data/ProcessedData/OptimalExtraction/optimal_extracted_science_flux.fits"
    lineParameterGuesses = initialThArGaussianParameterGuesses(pathThArSpectrum)

    spectraAllFibersAllOrders = getAllOptimalExtractedSpectrum(pathThArSpectrum)

    fiber = 5                      # ThAr fiber
    bestFitCoefficients = {}
    orders = []
    fittedGaussianCenters = np.array([])
    fittedGaussianStdevs = np.array([])
    fittedGaussianAmplitudes = np.array([])
    fittedGaussianOffsets = np.array([])
    chiSquares = np.array([])

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

        centers, stdevs, amplitudes, offsets, chiSquare, polynomialCoefficients =  \
                crudeWavelengthCalibrationOfOneOrder(xobs, fluxobs, centerGuesses, amplitudeGuesses, stdevGuesses, lambdaTrue)

        fittedGaussianCenters = np.append(fittedGaussianCenters, centers)
        fittedGaussianStdevs = np.append(fittedGaussianStdevs, stdevs)
        fittedGaussianAmplitudes = np.append(fittedGaussianAmplitudes, amplitudes)
        fittedGaussianOffsets = np.append(fittedGaussianOffsets, offsets)
        chiSquares = np.append(chiSquares, chiSquare)

        bestFitCoefficients[order] = polynomialCoefficients
        orders = orders + len(centers) * [order]

    # That's it!

    return fittedGaussianCenters, fittedGaussianStdevs, fittedGaussianAmplitudes, fittedGaussianOffsets, chiSquares, bestFitCoefficients, orders

    



def assessment():
    """ 
        
    INPUT:

    OUTPUT:
        linecenters_even: 
        velocities_even:
        linecenters_odd:
        velocities_odd: 

    """

    pathThArSpectrum = "Data/ProcessedData/OptimalExtraction/optimal_extracted_science_flux.fits"
    lineParameterGuesses = initialThArGaussianParameterGuesses(pathThArSpectrum)

    spectraAllFibersAllOrders = getAllOptimalExtractedSpectrum(pathThArSpectrum)

    fiber = 5                      # ThAr fiber

    unique_orders = []
    velocities_even = {}
    linecenters_even = {}
    velocities_odd = {}
    linecenters_odd = {}

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

        # Fit the even-indexed ThAr lines and compute a wavelength solution 

        centers_even, stdevs_even, amplitudes_even, offsets_even, chiSquare_even, polynomialCoefficients_even =  \
                crudeWavelengthCalibrationOfOneOrder(xobs, fluxobs, centerGuesses[::2], amplitudeGuesses[::2], stdevGuesses[::2], lambdaTrue[::2])

        linecenters_even[order] = centers_even

        # Fit the odd-indexed ThAr lines and compute a wavelength solution 

        centers_odd, stdevs_odd, amplitudes_odd, offsets_odd, chiSquare_odd, polynomialCoefficients_odd =  \
                crudeWavelengthCalibrationOfOneOrder(xobs, fluxobs, centerGuesses[1::2], amplitudeGuesses[1::2], stdevGuesses[1::2], lambdaTrue[1::2])

        linecenters_odd[order] = centers_odd

        # Verify how well we can predict the even-indexed ThAr lines using the wavelength calibration of the odd-indexed ThAr lines
 
        print("Predicting even-indexed ThAr lines using odd-index ThAr line wavelength calibration")

        X_even = np.vander(centers_even, len(polynomialCoefficients_odd), increasing=True)
        lambda_predicted = X_even @ polynomialCoefficients_odd
        velocity = (lambda_predicted - lambdaTrue[::2]) / lambdaTrue[::2] * 299792458.0
        velocities_even[order] = velocity
        
        # Verify how well we can predict the odd-indexed ThAr lines using the wavelength calibration of the even-indexed ThAr lines

        print("Predicting odd-indexed ThAr lines using even-index ThAr line wavelength calibration")

        X_odd = np.vander(centers_odd, len(polynomialCoefficients_even), increasing=True)
        lambda_predicted = X_odd @ polynomialCoefficients_even
        velocity = (lambda_predicted - lambdaTrue[1::2]) / lambdaTrue[1::2] * 299792458.0
        velocities_odd[order] = velocity

        # Keep the orders for the plot 

        unique_orders += [order]

        # That's it!

    return unique_orders, linecenters_even, velocities_even, linecenters_odd, velocities_odd















def plot(fitttedGaussianCenters, orders, chiSquares):
    cmap = plt.colormaps["cividis"]
    fig, ax = plt.subplots(figsize=(8,12))
    scatterplot = ax.scatter(fittedGaussianCenters, orders, s=8, color=cmap(chiSquares), norm=colors.LogNorm())
    fig.colorbar(scatterplot, ax=ax, location='bottom', pad=0.08, label="chi-square")
    ax.set(xlabel="Along-dispersion coordinate [pix]", ylabel="order")
    plt.show()



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
    




def plot_assessment(unique_orders, linecenters_even, velocities_even, linecenters_odd, velocities_odd):

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
    scatterplot = ax.scatter(linecenters, orders, s=8, c=cmap(velocities), vmin=vmin, vmax=vmax, cmap=cmap)
    cb = plt.colorbar(scatterplot, ax=ax, location='bottom', pad=0.08, shrink=0.9, label="velocity (m/s)")
    ax.set(xlabel="Along-dispersion coordinate [pix]", ylabel="order")
    plt.show()




if __name__ == "__main__":

    # fittedGaussianCenters, fittedGaussianStdevs, fittedGaussianAmplitudes, fittedGaussianOffsets,    \
    # chiSquares, bestFitCoefficients, orders = wavelength_calibration()

    unique_orders, linecenters_even, velocities_even, linecenters_odd, velocities_odd = assessment()
    plot_assessment(unique_orders, linecenters_even, velocities_even, linecenters_odd, velocities_odd)

