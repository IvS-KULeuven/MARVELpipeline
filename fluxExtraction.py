from pipeline import PipelineComponent
import numpy as np
import tools
import matplotlib.pyplot as plt 
from scipy import ndimage


class FluxExtraction(PipelineComponent):
    # Rename: Order extraction 


    def __init__(self, input):
        super().__init__(input)
        self.col = self.setCollection(input)
    
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
        print("Get the stripes:")

        # 1. We should get the image
        for hash in self.input:
            instances = self.col.find({"_id": hash})
            paths = [ x["path"] for x in instances ]
        image = tools.getImage(paths[0])

        # getStripes
        self.getStripes(image)



    def getStripes(self, image, deg_polynomial=5, median_filter=1, gauss_filter_sigma=3.,  min_peak=0.25):
        nx, ny = image.shape


        # Apply median filter and gaussian filter to the image  to reduce noise and improve algorithm stability.
        image = ndimage.median_filter(image, median_filter)
        image = ndimage.gaussian_filter(image, gauss_filter_sigma)

        # Central column
        centCol     = image[:, int(nx/2)]
        peaks    = np.r_[True, (centCol[1:-1] >= centCol[:-2]) & (centCol[1:-1] > centCol[2:]), True]

        peak_idx = np.arange(ny)[np.logical_and(peaks, centCol > min_peak * np.max(centCol))]
        # Exclude peaks too close to the border
        farFromBorder = lambda x : np.logical_and(x > 5, x < ny - 5)
        peak_idx = peak_idx[farFromBorder(peak_idx)]

        #============================================================================================
        print("{} number of peaks found".format(len(peak_idx)))

        plt.plot(centCol)
        plt.plot([ x if y else 0 for x, y in zip(centCol, peaks)])
        plt.show()
        #============================================================================================

        orders = [0]*len(peak_idx)
        for m, row_max in enumerate(peak_idx):
            orders[m] = self.followOrders(row_max, image)
        print(orders[0])
        #============================================================================================
        plt.plot(orders[0])
        plt.show()
        #============================================================================================

        polynomials = [np.poly1d(np.polyfit(np.arange(nx), order, deg_polynomial)) for order in orders]
        print(polynomials)





    def followOrders(self, max_row_0, image):
        nx, ny = image.shape
        order = np.zeros(nx)

        row_max  = max_row_0
        column = int(nx/2)

        getNeighbourIndices = lambda x : np.array([x-1, x, x+1]) if (x>2 and x<ny-2) else (np.array([x-1, x]) if x>2 else np.array([x, x+1]) )

        # Add center value to order
        order[column] = image[row_max, column]

        # Walk to right left
        while column+1 < nx:
            column += 1
            rows   = getNeighbourIndices(row_max)
            values = image[rows, column]
            
            row_max = rows[values == np.max(values)]
            order[column] = image[row_max, column]
            if (image[row_max, column] == 0):
                print("0")

        # Reset column and row_max and walk to the left
        column = int(nx/2)
        row_max = max_row_0
        
        while column > 0:
            column += -1
            rows = getNeighbourIndices(row_max)
            values = image[rows, column]

            row_max = rows[values == np.max(values)]
            order[column] = image[row_max, column]

        return order

        
            










if __name__ == "__main__":
    hash = ["e5577b3329caa7c3b98143a6dc0593b89ab5e11ee95b4d065f2cee210f81644e"]
    FExtra = FluxExtraction(hash)
    FExtra.make()
