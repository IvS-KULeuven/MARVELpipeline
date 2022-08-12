from optimalExtraction import OptimalExtraction
import numpy as np

class OptimalScienceExtraction(OptimalExtraction):
    """
    This is where we use the optimal extraction methond on the Extracted Science Images. 
    """
    def __init__(self, input, debug=0):
        """
        Initializes the component
        """
        super().__init__(input, "Science", debug=debug)


    def getFileName(self):
        return "optimal_extracted_science_flux.fits"

    def checkInput(self, input):
        if not super().checkInput(input):
            return False

        else:
            collection = self.db["ExtractedOrders"]
            amountOfScience = 0
            for hash in input:
                isExtractedScienceOrder = len([x for x in collection.find({"_id" : hash})]) == 1
                if isExtractedScienceOrder:
                    isExtractedScience = np.all([x["type"] == "Extracted Science Orders" for x in collection.find({"_id": hash})])
                    if isExtractedScience:
                        amountOfScience += 1
                        self.inputDict["science"] = hash
        return amountOfScience == 1





if __name__ == "__main__":
    hash_list = ["2133e1778dd6f8208763eeeed3a5ae6efd525fabb33a5fdf8809bd77bf89bb2b", "626de973a22fe042b8355bfdb868260e1dc13cbbc411f4baaf6730b813e3a26d"]    
    test = OptimalScienceExtraction(hash_list, debug=3)
    test.run()
