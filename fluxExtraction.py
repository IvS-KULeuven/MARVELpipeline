import time
from pipeline import PipelineComponent
from debug import plotOrdersWithSlider, plotGIF
import numpy as np
import tools
import matplotlib.pyplot as plt 
from scipy import ndimage
from tqdm import tqdm
from numba import njit, jit, vectorize
import numba


class OrderExtraction(PipelineComponent):


    def __init__(self, input, debug=0):
        super().__init__(input)
        self.col = self.setCollection(input)
        self.debug = debug
    
    def setCollection(self, input):
        return self.db["FlatImages"]

    def checkInput(self, input):
        isCorrect = []
        collection = self.db["FlatImages"]
        for hash in input:
            instances = collection.find({"_id" : hash})
            isCorrect.append(np.all([( x["type"] == "Master Flat Image") for x in instances]))
        return np.all(isCorrect)


    def make(self):

        # 1. We should get the image
        for hash in self.input:
            instances = self.col.find({"_id": hash})
            paths = [ x["path"] for x in instances ]
        image = tools.getImage(paths[0]).astype('float64')

        # 2. Locate the orders on the return a polynomial fit of them
        print("Find the stripes:")
        polynomials = self.getStripes(image, debug=True)

        # 3. Identify the stripes
        print("\nIdentify the Stripes")
        id_p = self.identifyStripes(image, polynomials )

        # 4. Extract the stripes
        print("\nExtract the Stripes")
        flat_stripes, index_fiber, index_order = self.extractFlatStripes(image, id_p)
        



    def getStripes(self, image, deg_polynomial=5, median_filter=1, gauss_filter_sigma=3.,  min_peak=0.125, debug=False):
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

        previous_idx = 0
        for m, row_max in tqdm(enumerate(peak_idx)):
            values[m], orders[m] = followOrders(row_max, previous_idx, image)
            previous_idx = row_max

        if self.debug > 2:
            print("{} number of peaks found".format(len(peak_idx)))
            plt.plot(centRow, "r")
            plt.plot(peak_idx, [centRow[i] for i in peak_idx], "bo")
            plt.show()           

            plotOrdersWithSlider(values)

            # Plot the selected pixels
            plt.imshow(image, origin='lower')
            for i in np.arange(320):
                ords = orders[i]
                mask = (ords != 0)
                plt.plot(ords[mask], np.arange(nx)[mask], alpha=0.5)
            plt.show()
        

        polynomials = [np.poly1d(np.polyfit(np.arange(nx)[order != 0], order[order != 0], deg_polynomial)) for order in orders]
        end = time.time()
        
        if self.debug > 2:
            print("Time for the function GetStripes is {}".format(end-start))
            xx = np.arange(nx)
            for p in polynomials:
                plt.plot(p(xx), xx)
            plt.imshow(image, alpha=1, origin='lower')
            plt.xlim([0, nx])
            plt.ylim([0, ny])
            plt.show()

        return polynomials









    def identifyStripes(self, image, polynomials, positions=None, selected_fibers=None):
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

        # Store the oberserved y_values of the middle column
        yObserved = np.array([np.poly1d(p)(int(nx/2)) for p in polynomials])

        shift_calculated = getShift(yPositions, yObserved, useAllFibers=useAllFibers)



        if self.debug > 2:
            
            plt.figure()
            plt.title("Stripe positions from configuration file")
            plt.imshow(image, origin='lower', vmin=np.min(image), vmax=0.5 * np.max(image))
            plt.plot(yObserved, np.repeat(nx / 2, len(yPositions)), 'b+', label='original stripe positions')
            plt.plot(yPositions + shift_calculated, np.repeat(nx / 2, len(yPositions)), 'g+',
                     label=f'corrected stripe positions ({shift_calculated:+.2f})')
            plt.legend()
            plt.show()

        polynomials = np.array([ p.coefficients for p in polynomials ])
        end = time.time()
        print("Time for the function identifyStripes is {}".format(end-start))
        return identify(yPositions, yObserved, polynomials, fibers, orders, shift_calculated)


    def extractFlatStripes(self, flat, p_id, slit_height=10):
        start = time.time()
        print("extract flat field stipes...")
        nx, ny = flat.shape
        xx = np.arange(nx)

        index_fiber = np.zeros_like(flat, dtype=np.int8)
        index_order = np.zeros_like(flat, dtype=np.int8)

        slit_indices_y = np.arange(-slit_height, slit_height).repeat(nx).reshape((2 * slit_height, nx))
        slit_indices_x = np.tile(np.arange(nx), 2 * slit_height).reshape((2 * slit_height, nx))

        for f in p_id.keys():
            for o, p in p_id[f].items():
                y = np.poly1d(p)(xx)
                indices = np.rint(slit_indices_y + y).astype(int)
                # Exclude the indices that are no longer on the CCD 
                valid_indices = np.logical_and(indices < ny, indices > 0)
                index_fiber[slit_indices_x[valid_indices], indices[valid_indices]] = f
                index_order[slit_indices_x[valid_indices], indices[valid_indices]] = o

        # image with only values within stripes, 0 elsewhere
        cleaned_image = np.where(index_order > 0, flat, 0)

        if self.debug > 2:
            fig, ax = plt.subplots(1, 3)
            ax[0].imshow(index_fiber, origin='lower')
            ax[1].imshow(index_order, origin='lower')
            ax[2].imshow(cleaned_image, origin='lower')
            plt.show()

        end = time.time()
        print("Time for first part of the function extractFlatStripes is {}".format(end-start))
        start = time.time()
        flat_stripes = self.extractStripes(flat, p_id)
        end = time.time()
        print("Time for second part of the function extractFlatStripes is {}".format(end-start))

        return flat_stripes, index_fiber, index_order
        

    def extractStripes(self, image, p_id):
        
        start = time.time()
        stripes = {}
        nx, ny = image.shape
        xx = np.arange(nx)

        images = []
        i = 0

        
        for f in p_id.keys():

            print(i)
            i += 1
            ind_img = np.zeros_like(image.transpose())

            for o, p in p_id[f].items():
                y  = np.poly1d(p)(xx)
                ind_img += extractSingleStripe(y, image)
            images.append(ind_img)
        if self.debug > 2:
            plotGIF(images, image)

        return images




        

        


    
@njit()
def getSignalToNoise(signal, background):
    gain = tools.getGain(" ")
    tExposure = tools.getExposureTime(" ")
    darkRate = tools.getDarkRate(" ")

    return (signal * gain) / np.sqrt((signal * gain + background * gain + tExposure * darkRate))

        
            


@njit()
def followOrders(max_row_0, previous_row, image):
    # Starting at the central peak, walk right/left and select the brightest pixel
    
    nx, ny = image.shape
    value = np.zeros(nx)
    order = np.zeros(nx)
    
    row_max  = max_row_0
    column = int(nx/2)

    getNeighbourIndices = lambda x : np.array([x-1, x, x+1]) if (x>1 and x<ny-2) else (np.array([x-1, x]) if x>1 else np.array([x, x+1]) )


    # Add center value to order/value
    value[column] = image[column, row_max]
    order[column] = row_max
    dark_value    = image[column, int((row_max + previous_row)/2)]
    # Walk to right left
    while column+1 < nx:
        column += 1

        rows   = getNeighbourIndices(row_max)
        values = np.array([image[column, row] for row in rows])

        row_max = rows[values == np.max(values)][0]
        value[column] = image[column, row_max]
        order[column] = row_max
        dark_value    = image[column, int((row_max + previous_row)/2)]
        


        if (row_max == 1) or (row_max == nx):
            break
        
        if getSignalToNoise(value[column], dark_value) < 100:
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
        dark_value    = image[column, int((row_max + previous_row)/2)]        

        if (row_max == 1) or (row_max == nx):
            break

        if getSignalToNoise(value[column], dark_value) < 100:
            break

    return value, order


@njit()
def getShift(positions, observed, useAllFibers=True):

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

    

def identify(positions, observed, polynomials, fibers, orders, shift):
    print("Beginning of identify")

    p_id = {}
            
    # Keeps track of the orders that have been identified
    used = np.zeros_like(positions)

    for i, p in enumerate(polynomials):
        closest_stripe_idx = np.argmin( np.abs(positions + shift - observed[i]) )
        if np.abs(positions[closest_stripe_idx] + shift - observed[i]) < 7:
            if used[closest_stripe_idx] == 0:
                used[closest_stripe_idx] = 1
                fiber = fibers[closest_stripe_idx]
                order = orders[closest_stripe_idx]
                
                if fiber in p_id:
                    p_id[fiber].update({order: p})
                else:
                    p_id[fiber] = {order: p}
            else:
                print("WARNING: Stripe at {} could not be identified unambiguously".format(observed[i]))
        else:
            print("Stripe at {} could not be identified.".format(observed[i]))
    print("End of identify")
    return p_id





@njit()
def extractSingleStripe(y, image, slit_height=10):

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
    








if __name__ == "__main__":
    hash = ["e5577b3329caa7c3b98143a6dc0593b89ab5e11ee95b4d065f2cee210f81644e"]
    FExtra = OrderExtraction(hash, debug=3)
    FExtra.make()
