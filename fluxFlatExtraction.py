import time
from pipeline import PipelineComponent
from debug import plotOrdersWithSlider, plotGIF
import numpy as np
import tools
import matplotlib.pyplot as plt 
from scipy import ndimage
from tqdm import tqdm
from numba import njit, jit, vectorize
import hashlib
from astropy.io import fits
import os


class FlatOrderExtraction(PipelineComponent):
    """
    DESCRIPTION: This module identifies the correct pixels and their
                 corresponding flux values for every fiber and order.

    INPUT: Input hashes correspond to one Master Flat image.

    ALGORITHM:
    1. Take the middel row of the flat field and identifies the stripes
       by finding the maxima of the flux.
    2. Follow these maxima in the upwards and downwards direction and save
       the relevant (high enough S/N) positions.
    3. Estimate these lines on the flatfield with 5th order polynomial.
    4. Label the lines with their corresponding fiber and order.
    5. For every order, extract the relevant pixels around this polynomial
       fit and save x-values, y-values and flux values.
    """

    def __init__(self, input, debug=0):
        """
        Initialize the componont
        """
        super().__init__(input)
        self.col = self.setCollection(input)
        self.debug = debug
        self.outputPath = os.getcwd() + "/Data/ProcessedData/ExtractedOrders/"

        if not os.path.isdir(self.outputPath):
            os.mkdir(self.outputPath)

        self.type = "Extracted Flat Orders"



    def setCollection(self, input):
        """
        Where in the database we should look for the input hashes
        """
        return self.db["FlatImages"]



    def checkInput(self, input):
        isCorrect = []
        collection = self.db["FlatImages"]
        for hash in input:
            instances = collection.find({"_id" : hash})
            isCorrect.append(np.all([( x["type"] == "Master Flat Image") for x in instances]))
        return np.all(isCorrect)



    def make(self):
        """
        Runs the algorithm described after the class definition.
        """
        # 1. We should get the image
        for hash in self.input:
            instances = self.col.find({"_id": hash})
            paths = [ x["path"] for x in instances ]
        image = tools.getImage(paths[0]).astype('float64')

        # 2. Locate the orders on the return a polynomial fit of them
        if (self.debug > 2):
            print("Find the stripes:")
        x_values, polynomials = self.getStripes(image, debug=True)

        # 3. Identify the stripes
        if (self.debug > 2):
            print("\nIdentify the Stripes")
        id_p = self.identifyStripes(image, polynomials, x_values)

        # 4. Extract the stripes
        if (self.debug > 2):
            print("\nExtract the Stripes")
        return self.extractFlatStripes(image, id_p)



    def getStripes(self, image, deg_polynomial=5, median_filter=1, gauss_filter_sigma=3.,  min_peak=0.05, debug=False):
        """
        DESCRIPTION: Find the maxima values of the middle line of the flat image. Follows the image to the top and bottem 
                     and saves the position of the brightest (relevant) pixels. Returns a polynomial fit of degree 
                     deg_polynomials of the lines.

        INPUT: - image: Flat field image 
               - deg_polynomial: degree of the polynomial to fit the lines with 
               - median_filter=1, gauss_filter_sigma=3 = value of median filter and gaussian filter sigma that we apply to 
                                                         image to reduce noise. 
               - min_peak: maxima of the middle line lower then min_peak * max_value will not be taken into account in order 
                           to avoid false maxima. 
                           
        """
        start = time.time()
        nx, ny = image.shape

        # Apply median filter and gaussian filter to the image  to reduce noise and improve algorithm stability.
        image = ndimage.median_filter(image, median_filter)
        image = ndimage.gaussian_filter(image, gauss_filter_sigma)

        # Central row of the CCD
        centRow  = image[int(nx/2),:]
        peaks    = np.r_[True, (centRow[1:-1] >= centRow[:-2]) & (centRow[1:-1] > centRow[2:]), True]

        # Identify the maxima along the central row, as these are the stripes
        # We only keep maxima that are larger then min_peak * maximum_central_row
        peak_idx = np.arange(ny)[np.logical_and(peaks, centRow > min_peak * np.max(centRow))]

        # Exclude peaks too close to the border
        farFromBorder = lambda x : np.logical_and(x > 5, x < ny - 5)
        peak_idx = peak_idx[farFromBorder(peak_idx)]

        orders = [0]*len(peak_idx)
        values = [0]*len(peak_idx)

        first_fiber_idx = peak_idx[::5]
        final_fiber_idx = np.concatenate((np.zeros(1), peak_idx[4::5]))[:-1]

        for m, row_max in tqdm(enumerate(peak_idx)):
            if (row_max in first_fiber_idx):
                dark_idx = int((row_max + final_fiber_idx[first_fiber_idx == row_max])/2)

            values[m], orders[m] = followOrders(row_max, dark_idx, image)

        x_values = [ np.arange(nx)[ordr != 0] for ordr in orders]

        if self.debug > 2:
            print("\t{} number of peaks found".format(len(peak_idx)))
            plt.plot(centRow, "r")
            plt.plot(peak_idx, [centRow[i] for i in peak_idx], "bo")
            plt.show()           

            plotOrdersWithSlider(values)

            # Plot the selected pixels
            plt.imshow(image, origin='lower')
            for i in np.arange(330):
                ords = orders[i]
                mask = (ords != 0)
                plt.plot(ords[mask], np.arange(nx)[mask], alpha=0.5)
            plt.show()
        

        polynomials = [np.poly1d(np.polyfit(np.arange(nx)[order != 0], order[order != 0], deg_polynomial)) for order in orders]
        end = time.time()

        if self.debug > 2:
            xx = np.arange(nx)
            for x, p in zip(x_values, polynomials):
                plt.plot(p(x), x)
            plt.imshow(image, alpha=1, origin='lower')
            plt.xlim([0, nx])
            plt.ylim([0, ny])
            plt.show()
        return x_values, polynomials



    def run(self):
        """
        Ruins the alghoritm in main and saves the values in fits file and adds the image to the database. 
        """
        xCoordinates, yCoordinates, fluxValues, orders = self.make()
        self.saveImage(xCoordinates, yCoordinates, fluxValues, orders)
        print("Block Generated!")



    def saveImage(self, xValues, yValues, flux, orders):
        """
        Save the image and add to the database. 
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



    def getFileName(self):
        return "extracted_flat_orders.fits"



    def identifyStripes(self, image, polynomials, xValues, positions=None, selected_fibers=None):
        """
        Identify the stripes with their correct fiber and order label.
        """
        start = time.time()
        nx, ny = image.shape

        useAllFibers = selected_fibers is None

        # We store the expected positions of the orders on the CCD 
        # THIS MUST BE CHANGED LATER ON, TESE VALUES SHOULD BE INCLUDED IN THE FITS FILE,
        # FOR NOW THIS IS OBTAINED FROM THE tools.py MODULE
        positions = tools.getPositionsOfOrders()

        yPositions = np.array(positions[2])
        orders    = np.array(positions[1])
        fibers    = np.array(positions[0])

        if selected_fibers is None:
            selected_fibers = np.unique(fibers)

        for f in np.unique(fibers):
            if f not in selected_fibers:
                idx = (fibers != f)
                yPositions = yPositions[idx]
                fibers     = fibers[idx]
                orders     = orders[idx]

        # Store the oberserved y_values of the middle column.
        yObserved = np.array([np.poly1d(p)(int(nx/2)) for p in polynomials])

        shift_calculated = getShift(yPositions, yObserved, useAllFibers=useAllFibers)

        if self.debug > 2:
            plt.figure()
            plt.title("Stripe positions from configuration file")
            plt.imshow(image, origin='lower', vmin=np.min(image), vmax=0.5 * np.max(image))
            plt.plot(yObserved, np.repeat(nx / 2, len(yObserved)), 'b+', label='original stripe positions')
            plt.plot(yPositions + shift_calculated, np.repeat(nx / 2, len(yPositions)), 'g+',
                     label=f'corrected stripe positions ({shift_calculated:+.2f})')
            plt.legend()
            plt.show()

        polynomials = np.array([ p.coefficients for p in polynomials ])

        if self.debug > 2:
            end = time.time()
            print("\tTime for the function identifyStripes is {}".format(end-start))
        return identify(yPositions, yObserved, polynomials, xValues, fibers, orders, shift_calculated)



    def extractFlatStripes(self, flat, p_id):
        """
        Extract the relevant values of the image for every order and fiber.
        """
        start = time.time()
        nx, ny = flat.shape

        index_fiber = np.zeros_like(flat, dtype=np.int8)
        index_order = np.zeros_like(flat, dtype=np.int8)

        cross_widths = {}

        previous_idx = 0
        for f in p_id.keys():
            for o, (xx, p) in p_id[f].items():
                xSize = np.size(xx)
                y = np.poly1d(p)(xx)
                #cross_width = getCrossOrderWidth(flat, int(np.poly1d(p)([int(nx)/2])), previous_idx)
                cross_width = 6

                slit_indices_y = np.arange(-cross_width, cross_width+1).repeat(xSize).reshape((2 * cross_width+1, xSize))
                slit_indices_x = np.tile(xx, 2 * cross_width+1).reshape((2 * cross_width+1, xSize))                

                indices = np.rint(slit_indices_y + y).astype(int)

                # Exclude the indices that are no longer on the CCD 
                valid_indices = np.logical_and(indices < ny, indices > 0)
                index_fiber[slit_indices_x[valid_indices], indices[valid_indices]] = f
                index_order[slit_indices_x[valid_indices], indices[valid_indices]] = o
                previous_idx = int(np.poly1d(p)([int(nx)/2]))

                if f in cross_widths:
                    cross_widths[f].update({o: cross_width})
                else:
                    cross_widths[f] = {o: cross_width}

        # image with only values within stripes, 0 elsewhere
        # TODO: Check if we need cleaned image
        cleaned_image = np.where(index_order > 0, flat, 0)

        if self.debug > 2:
            fig, ax = plt.subplots(1, 3)
            ax[0].imshow(index_fiber, origin='lower')
            ax[1].imshow(index_order, origin='lower')
            ax[2].imshow(cleaned_image, origin='lower')
            plt.show()

        if self.debug > 1:
            plt.imshow(flat, origin='lower')
            plt.imshow(index_fiber, alpha=0.3, origin='lower')
            plt.show()


        if self.debug > 2:
            end = time.time()
            print("\tTime for first part of the function extractFlatStripes is {}".format(end-start))
            start = time.time()

        xCoordinates, yCoordinates, fluxValues, orders = extractStripes(flat, index_fiber, index_order)

        if self.debug > 2:
            end = time.time()
            print("\tTime for second part of the function extractFlatStripes is {}".format(end-start))
        return xCoordinates, yCoordinates, fluxValues, orders
        



    
#@njit()
def extractStripes(image, fiber_indexes, order_indexes):
    """
    Docstring
    """

    xCoordinates = [] 
    yCoordinates = []
    fluxValues   = []
    orders       = []
    
    nRow, nCol = image.shape
    columns    = np.repeat(np.arange(nRow), nCol).reshape((nRow, nCol))
    for o in np.arange(33, np.max(order_indexes)+1):
        order_mask = order_indexes == o

        for f in np.arange(1, np.max(fiber_indexes[order_mask])+1):

            fiber_mask = fiber_indexes == f
            combined_mask = fiber_mask * order_mask
            xValues, yValues, flux = extractSingleStripe2(image, combined_mask, columns)
            xCoordinates.append(xValues)
            yCoordinates.append(yValues)
            fluxValues.append(flux)
            orders.append((o,f))

    return xCoordinates, yCoordinates, fluxValues, orders



@njit()
def extractSingleStripe2(image, mask, columns):
    xCol = np.extract(mask, columns)
    yCol = np.extract(mask, np.transpose(columns))
    flux = np.extract(mask, image)
    return xCol, yCol, flux








@njit()
def getSignalToNoiseSinglePixel(signal, background):
    """
    Docstring
    """
    gain = tools.getGain(" ")
    tExposure = tools.getExposureTime(" ")
    darkRate = tools.getDarkRate(" ")

    return (signal * gain) / np.sqrt((signal * gain + background * gain + tExposure * darkRate))



        
            


@njit()
def followOrders(max_row_0, dark_column, image):
    """
    Starting at the central peak, walk right/left and select the brightest pixel
    """

    nx, ny = image.shape
    value = np.zeros(nx)
    order = np.zeros(nx)

    row_max  = max_row_0
    column = int(nx/2)
    
    getNeighbourIndices = lambda x : np.array([x-1, x, x+1]) if (x>1 and x<ny-2) else (np.array([x-1, x]) if x>1 else np.array([x, x+1]) )


    # Add center value to order/value
    value[column] = image[column, row_max]
    order[column] = row_max
    dark_value    = image[column, dark_column]
    
    # Walk to right left
    while column+1 < nx:
        column += 1
        rows   = getNeighbourIndices(row_max)
        values = np.array([image[column, row] for row in rows])

        row_max = rows[values == np.max(values)][0]

        value[column] = image[column, row_max]
        order[column] = row_max
        dark_value    = image[column, dark_column]


        if (row_max == 1) or (row_max == nx):
            break

        if getSignalToNoiseSinglePixel(value[column], dark_value) < 20:
            break

    # Reset column and row_max and walk to the left
    column = int(nx/2)
    row_max = max_row_0

    while column > 0:
        column += -1

        rows = getNeighbourIndices(row_max)
        values = np.array([image[column, row] for row in rows])

        row_max = rows[values == np.max(values)][0]

        value[column] = image[column, row_max]
        order[column] = row_max
        dark_value    = image[column, dark_column]

        if (row_max == 1) or (row_max == nx):
            break

        if getSignalToNoiseSinglePixel(value[column], dark_value) < 20:
            break

    return value, order




@njit()
def getShift(positions, observed, useAllFibers=True):
    """
    Docstring
    """
    shifts = np.linspace(-200, 200, 20000)
    distanceForShift = lambda shift : np.array([ np.min( np.abs(positions + shift - y)) for y in observed])
    distanceForAllShifts = np.array( [ distanceForShift(shift).sum() for shift in shifts])

    # Keeps track of the orders that have been identified
    used = np.zeros_like(positions)

    # This is important to correctly label the fiber IDs when not all fibers are used. We weigh
    # closer to the initial guess. 
    if not useAllFibers:
        distanceForAllShifts += np.abs(shifts) * 2
        
    shift = shifts[np.argmin(distanceForAllShifts)]
    return shift

    

def identify(positions, observed, polynomials, xValues, fibers, orders, shift):
    """
    Docstring
    """
    p_id = {}
    # Keeps track of the orders that have been identified
    used = np.zeros_like(positions)

    for i, (x, p) in enumerate(zip(xValues, polynomials)):
        closest_stripe_idx = np.argmin( np.abs(positions + shift - observed[i]) )
        if np.abs(positions[closest_stripe_idx] + shift - observed[i]) < 7:
            if used[closest_stripe_idx] == 0:
                used[closest_stripe_idx] = 1
                fiber = fibers[closest_stripe_idx]
                order = orders[closest_stripe_idx]
                
                if fiber in p_id:
                    p_id[fiber].update({order: (x, p)})
                else:
                    p_id[fiber] = {order: (x, p)}
            else:
                print("WARNING: Stripe at {} could not be identified unambiguously".format(observed[i]))
        else:
            print("Stripe at {} could not be identified.".format(observed[i]))
    return p_id





@njit()
def extractSingleStripe(y, image, slit_height=15):
    """
    Docstring
    """
    ny, nx = image.shape

    split_indices_y = np.arange(-slit_height, slit_height).repeat(nx).reshape((2 * slit_height, nx))
    split_indices_x = np.arange(nx).repeat(2 * slit_height).reshape((nx, 2 * slit_height)).transpose()

    indices = np.empty_like(split_indices_y)
    np.round_(split_indices_y + y, 0, indices)

    valid_indices = np.logical_and(indices < ny, indices > 0)


    ind_img = np.zeros_like(image.transpose())
    for i in np.arange(valid_indices.size):
        if valid_indices.flat[i]:
            ind_img[indices.flat[i], split_indices_x.flat[i]] = 1

    return ind_img.transpose()
    



#@njit()
def getCrossOrderWidth(image, idx, idx_previous):
    """
    Docstring
    """
    # Get the optimal slith height (<20, >0) that maximaizes S/N

    nx, ny = image.shape
    centralRow = image[int(nx/2),:]

    background = centralRow[int((idx_previous + idx)/2)]
    previousSignalToNoise = 0
    noise = []
    for width in np.arange(15):
        
        signal     = sum(centralRow[idx-width:idx+width+1])
        background = background * (width * 2 + 1)


        currentSignalToNoise = getSignalToNoise(signal, background, (2 * width + 1))
        noise.append(currentSignalToNoise)
        if (currentSignalToNoise >= previousSignalToNoise):
            previousSignalToNoise = currentSignalToNoise
        else:
            return width + 2
    return 15



@njit()
def getSignalToNoise(signal, background, Npixels):
    """
    Docstring
    """
    gain = tools.getGain(" ")
    tExposure = tools.getExposureTime(" ")
    darkRate = tools.getDarkRate(" ")

    return (signal * gain) / np.sqrt((signal * gain + background * gain + tExposure * darkRate * Npixels))







if __name__ == "__main__":

    masterflat_hash = ["641f2b56be1a8a86848c29abbd81858ddd15e42359ae84473050962efb9dea06"]
    fluxExtraction = FlatOrderExtraction(masterflat_hash, debug=1)


    fluxExtraction.run()





