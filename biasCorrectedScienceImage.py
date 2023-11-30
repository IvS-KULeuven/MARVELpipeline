from pipeline   import PipelineComponent
from database import DatabaseFromLocalFile

import os
import tools
import hashlib
import numpy as np




class BiasCorrectedScienceFrames(PipelineComponent):
    """
    TODO: documentation. What is this class for?
    """

    def __init__(self, database=None, **inputScienceHashes):
        super().__init__(database, **inputScienceHashes)
        inputScienceHashes = self.inputHashes

        if self.checkSanityOfInputTypes(**inputScienceHashes):
            self.outputPath = os.getcwd() + "/Data/ProcessedData/BiasCorrectedScience/"
            self.type = "Bias Corrected Science Image"
            self.masterBiasHash = inputScienceHashes["BiasImages"]
            self.rawScienceHash = inputScienceHashes["ScienceImages"]
        else:
            raise Exception("Error: The input hashes do not match the correct type: Aborting")
            exit(1)

        if not os.path.isdir(self.outputPath):
            os.mkdir(self.outputPath)







    def checkSanityOfInputTypes(self, **inputScienceHashes):
        """
        This function is ran after we run checkSanityOfInputHashes. This function checks the the
        input types that are given is able to generate a bias corrected science output file.
        """
        types  = list(inputScienceHashes.keys())
        values = list(inputScienceHashes.values())

        # Check that the keys are of the right format. For a bias corrected science image, these should be science
        # images, dark images and bias images.

        correctKeysFormat = ["ScienceImages", "BiasImages"]
        keysAreCorrect = (len(types) == 2) and np.all([key in types for key in correctKeysFormat])

        if keysAreCorrect:
            valuesAreCorrectType = np.all([ isinstance(value, str) for value in values])
        else:
            return False

        if valuesAreCorrectType:

            isMasterBiasImage = self.db["BiasImages"].find_one({"_id": inputScienceHashes["BiasImages"]})["type"] == "Master Bias Image"
            isRawScienceImage = self.db["ScienceImages"].find_one({"_id": inputScienceHashes["ScienceImages"]})["type"] == "Raw Science Image"
        else:
            return False

        return isMasterBiasImage and isRawScienceImage








    def run(self, outputFileName=None):
        """
        We run through the algorithm to create the bias corrected science images.

        Input:
            outputFileName: If None, nothing is saved. Otherwise, a string with the name of the outputfile,
                            incl. the extension ".fits".

        Output:
            calibratedScience:     bias corrected science image [ADU]
        """

        # Get all the paths of the files corresponding to these hashes
        biasPath    = self.db["BiasImages"].find_one({"_id": self.masterBiasHash })["path"]
        sciencePath = self.db["ScienceImages"].find_one({"_id": self.rawScienceHash})["path"]

        # Get all the fits files corresponding to these hashes
        bias    = tools.getImage(biasPath)
        science = tools.getImage(sciencePath)

        biasCorrectedScience = science - bias

        # Add offset so that all the values in the MasterFlat are positive
        if np.min(biasCorrectedScience) < 0:
            biasCorrectedScience = biasCorrectedScience - np.min(biasCorrectedScience)


        if outputFileName is not None:
            self.saveImageAndAddToDatabase(biasCorrectedScience, outputFileName)
            print("Bias Corrected science image saved to fits file")

        # That's it!

        print("Block generated!")
        return biasCorrectedScience






    def getHashOfOutputfile(self):
        """
        The function returns the hash id for the output file.
        This hash is made from the hash files that are used as input.

        Ouput:
           hash. string containing the output hash
        """
        combinedHashes =  self.masterBiasHash + self.rawScienceHash
        hash = hashlib.sha256(bytes(combinedHashes, 'utf-8')).hexdigest()
        return hash



if __name__ == "__main__":
    databaseName = "pipelineDatabase.txt"
    print("Creating a local database file with the name: ", databaseName)

    db = DatabaseFromLocalFile(databaseName)

    # bias corrected Science Image
    rawScienceHash = "d2ad34be4443931a5128f3e035ce12c869fb801c017f856cb3c68c1e180d1267"

    master_bias_pathF = "Data/ProcessedData/MasterBias/testsFBias.fits"


    print(" ")
    calibration1 = BiasCorrectedScienceFrames(db, ScienceImages=rawScienceHash, BiasImages=master_bias_pathF)

    calibration1.run("testFScience.fits")
    db.save()
