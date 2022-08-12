from optimalExtraction import OptimalExtraction
import numpy as np


class OptimalEtalonExtraction(OptimalExtraction):
    """
    This is where we use the optimal extraction methond on the Extracted Etalon Images.
    """
    def __init__(self, input, debug=0):
        """
        Initializes the component
        """
        super().__init__(input, "Etalon", debug=debug)



    def getFileName(self):
        return "optimal_extracted_etalon_flux.fits"



    def checkInput(self, input):

        if not super().checkInput(input):
            return False
        else:
            collection = self.db["ExtractedOrders"]
            amountOfEtalon = 0
            for hash in input:
                isExtractedEtalonOrder = len([ x for x in collection.find({"_id" : hash})]) == 1
                if isExtractedEtalonOrder:
                    isExtractedEtalon = np.all([x["type"] == "Extracted Etalon Orders" for x in collection.find({"_id" : hash})])
                    if isExtractedEtalon:
                        amountOfEtalon += 1
                        self.inputDict["etalon"] = hash
        return amountOfEtalon == 1






if __name__ == "__main__":
    hash_list = ["2133e1778dd6f8208763eeeed3a5ae6efd525fabb33a5fdf8809bd77bf89bb2b", "f9900df2548dba0b438594192d6df0083ee954facfa3869cad7344262de23533"]
    test = OptimalEtalonExtraction(hash_list, debug=3)
    test.run()
