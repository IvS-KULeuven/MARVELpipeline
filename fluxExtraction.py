from pipeline import PipelineComponent
from debug import plotOrdersWithSlider
import numpy as np
import tools
import matplotlib.pyplot as plt 
from scipy import ndimage
from tqdm import tqdm



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
        image = tools.getImage(paths[0])

        # get Stripes
        print("Get the stripes:")
        polynomials = self.getStripes(image, debug=True)

        # identify stripes
        id_p = self.identifyStripes(image, polynomials, selected_fibers=[1])
        



    def getStripes(self, image, deg_polynomial=5, median_filter=1, gauss_filter_sigma=3.,  min_peak=0.25, debug=False):
        nx, ny = image.shape

        # Apply median filter and gaussian filter to the image  to reduce noise and improve algorithm stability.
        image = ndimage.median_filter(image, median_filter)
        image = ndimage.gaussian_filter(image, gauss_filter_sigma)

        # Central column
        centRow  = image[int(nx/2),:]
        peaks    = np.r_[True, (centRow[1:-1] >= centRow[:-2]) & (centRow[1:-1] > centRow[2:]), True]
        
        peak_idx = np.arange(ny)[np.logical_and(peaks, centRow > min_peak * np.max(centRow))]

        # Exclude peaks too close to the border
        farFromBorder = lambda x : np.logical_and(x > 5, x < ny - 5)
        peak_idx = peak_idx[farFromBorder(peak_idx)]

        orders = [0]*len(peak_idx)


        for m, row_max in tqdm(enumerate(peak_idx)):
            orders[m] = self.followOrders(row_max, image)

        if self.debug > 2:
            print("{} number of peaks found".format(len(peak_idx)))
            plt.plot(centRow, "r")
            plt.plot(peak_idx, [centRow[i] for i in peak_idx], "bo")
            plt.show()            

            plotOrdersWithSlider(orders)
        

        polynomials = [np.poly1d(np.polyfit(np.arange(nx), order, deg_polynomial)) for order in orders]
        return polynomials






    def followOrders(self, max_row_0, image):
        nx, ny = image.shape
        order = np.zeros(nx)

        row_max  = max_row_0
        column = int(nx/2)

        # getNeighbourIndices = lambda x : np.array([x-1, x, x+1]) if (x>1 and x<ny-2) else (np.array([x-1, x]) if x>1 else np.array([x, x+1]) )

        def getNeighbourIndices(x):

            indicesToConsided = np.array([x-2, x-1, x, x+1, x+2])

            if (x>2 and x<ny-3):
                return indicesToConsided
            elif (x>2):
                return indicesToConsided[:2*ny - x - 5]
            else:
                return indicesToConsided[x:]


        # Add center value to order
        order[column] = image[column, row_max]

        # Walk to right left
        while column+1 < nx:

            column += 1

            rows   = getNeighbourIndices(row_max)
            values = image[column, rows]
            
            row_max = rows[values == np.max(values)][0]
            order[column] = image[column, row_max]


        # Reset column and row_max and walk to the left
        column = int(nx/2)
        row_max = max_row_0
        
        while column > 0:
            column += -1

            rows = getNeighbourIndices(row_max)
            values = image[column, rows]

            row_max = rows[values == np.max(values)][0]
            order[column] = image[column, row_max]

        return order



    def identifyStripes(self, image, polynomials, positions=None, selected_fibers=None):

        print("Identify stripes ...")
        nx, ny = image.shape
        p_id = {}

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
        yObserved = []
        for i, p in enumerate(polynomials):
            yObserved.append(np.poly1d(p)(int(nx / 2)))

                
        # What are these shifts for?
        shifts = np.linspace(-200, 200, 2000)
        distanceForShift = lambda shift : [ np.min( np.abs(yPositions + shift - y)) for y in yObserved ]
        totalDistance = [ np.array(distanceForShift(shift)).sum() for shift in shifts]

        shift_calculated = shifts[np.argmin(totalDistance)]
        print("shift: ", shift_calculated)

        

        # print(yObserved)

        
            










if __name__ == "__main__":
    hash = ["e5577b3329caa7c3b98143a6dc0593b89ab5e11ee95b4d065f2cee210f81644e"]
    FExtra = OrderExtraction(hash, debug=2)
    FExtra.make()
