from pipeline   import PipelineComponent
from database   import DatabaseFromLocalFile
import os
import tools
import numpy as np
import hashlib






class MasterFlat(PipelineComponent):

    def __init__(self, database=None, **flatAndBiasHashes):

        super().__init__(database, **flatAndBiasHashes)
        flatAndBiasHashes = self.inputHashes

        if self.checkSanityOfInputTypes(**flatAndBiasHashes):
            self.outputPath = os.getcwd() + "/Data/ProcessedData/MasterFlat/"
            self.type = "Master Flat Image"
            self.masterBiasHash = flatAndBiasHashes["BiasImages"]
            self.rawFlatHashes  = flatAndBiasHashes["FlatImages"]
        else:
            raise Exception("Error: The input hashes do not match the correct type: Aborting")
            exit(1)

        if not os.path.isdir(self.outputPath):
            os.mkdir(self.outputPath)






    def checkSanityOfInputTypes(self, **flatAndBiasHashes):
        """
        This function is ran after we run checkSanityOfInputHashes. This function checks the the
        input types that are given is able to generate a master flat output file.
        """

        types  = list(flatAndBiasHashes.keys())
        #values = list(flatAndBiasHashes.values())

        # Check that the keys are of the right format. For a master dark these should only be DarkImages

        keysAreCorrect = (len(types)) == 2 and ("BiasImages" in types) and ("FlatImages" in types)

        if keysAreCorrect:
            valuesAreCorrect = (len(flatAndBiasHashes["FlatImages"]) > 1) and isinstance(flatAndBiasHashes["BiasImages"], str)
        else:
            return False

        if valuesAreCorrect:
            areRawFlatImages  = np.all([ self.db["FlatImages"].find_one({"_id": hash})["type"] == "Raw Flat Image"
                                        for hash in flatAndBiasHashes["FlatImages"] ])
            isMasterBiasImage = self.db["BiasImages"].find_one({"_id": flatAndBiasHashes["BiasImages"]})["type"] == "Master Bias Image"
        else:
            return False

        # If we get here, the input is sane if the hashes correspond to raw flat images or one master bias image

        return areRawFlatImages and isMasterBiasImage









    def run(self, outputFileName=None):
        """
        We run through the alghoritm to create the master flat images.

        Input:
            outputFileName: If None, nothing is saved. Otherwise, a string with the name of the outputfile,
                            incl. the extension ".fits".

        Output:
            masterFlat:     master flat image [ADU]
        """

        # Get all the paths of the files corresponding to these hashes
        flatPaths = [ (self.db["FlatImages"].find_one({"_id": hash}))["path"] for hash in self.rawFlatHashes ]
        biasPath  = self.db["BiasImages"].find_one({"_id": self.masterBiasHash})["path"]

        # Get all the fits files corresponding to these hashes
        flats = tools.getImages(flatPaths)
        bias  = tools.getImage(biasPath)

        # Use the image in the fits files, and use mean_combining to obtain the the master image
        masterFlat = np.median(flats, axis=0) - bias

        # Add offset so that all the values in the MasterFlat are positive
        if np.min(masterFlat) < 0:
                  masterFlat = masterFlat -np.min(masterFlat)

        if outputFileName is not None:
            self.saveImageAndAddToDatabase(masterFlat, outputFileName)
            print("Master flat image saved to fits file")

        # That's it!

        print("Block generated!")
        return masterFlat







    def getHashOfOutputfile(self):
        """
        The function returns the hash id for the output file.
        This hash is made from the hash files that are used as input.

        Ouput:
           hash. string containing the output hash
        """
        hash = hashlib.sha256(bytes("".join(self.rawFlatHashes) + self.masterBiasHash, 'utf-8')).hexdigest()
        return hash


if __name__ == "__main__":
    databaseName = "pipelineDatabase.txt"
    print("Creating a local database file with the name: ", databaseName)

    db = DatabaseFromLocalFile(databaseName)
    print("")

        # Mater Flat Image
    raw_flat_hashes = ["e003c1a1f66829c11bef1cfa35d4391a87b89cb0e35dcbdbde0c9acfff19f64b",
                       "844a3574b341fec134c3c6bd7464d4dfa1a19da94c477f42e65484be90676d08",
                       "9b2902d98bf14d627bf60182b6974eb86d6b313652a65188724c4a9c521ab0c8",
                       "e25e03dacdfcaa1da19a5f95261def79c1cc806588bca507446264b8fb0fe24a"]

    master_bias_pathF = "Data/ProcessedData/MasterBias/testsFBias.fits"
    masterF1 = MasterFlat(db, FlatImages=raw_flat_hashes, BiasImages=master_bias_pathF)
    masterF1.run("testsFFlat.fits")
    db.save()
