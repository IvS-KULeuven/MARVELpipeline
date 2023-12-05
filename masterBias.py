from pipeline   import PipelineComponent
from database   import DatabaseFromLocalFile

import yaml
import os
import tools
import numpy as np
import hashlib







class MasterBias(PipelineComponent):
    """
    Class that creates the master bias image. Such an image is obtained by taking the median of multiple
    raw bias images.
    """

    def __init__(self, database=None, **input):
        super().__init__(database, **input)
        input = self.inputHashes

        if self.checkSanityOfinputTypes(**input):

            self.outputPath = os.getcwd()
            self.type = "Master Bias Image"
            self.rawBiasHashes = input["BiasImages"]
        else:
            raise Exception("Error: The input hashes do not match the correct type: Aborting")
            exit(1)


        if not os.path.isdir(self.outputPath):
            os.mkdir(self.outputPath)





    def checkSanityOfinputTypes(self, **input):
        """
        This function is ran after we run checkSanityOfInputHashes. This function checks the the
        input types that are given is able to generate a master bias output file.
        """
        types  = list(input.keys())
        values = list(input.values())

        # Check that the keys are of the right format. For a master bias these should only be BiasImages

        keysAreCorrect = (len(types) == 1) and types[0] == "BiasImages"

        valuesAreCorrect = (len(values) == 1) and (len(values[0]) > 1)

        if keysAreCorrect and valuesAreCorrect:
            areRawImages = np.all([ self.db[types[0]].find_one({"_id": hash})["type"] == "Raw Bias Image" for hash in values[0]])
        else:
            return False

        # If we get here, the input is sane if the hashes correspond to raw Images

        return areRawImages






    def run(self, outputFileName=None):
        """
        We run through the alghoritm to create the master bias images.

        Input:
            outputFileName: If None, nothing is saved. Otherwise, a string with the name of the outputfile,
                            incl. the extension ".fits".

        Output:
            masterBias:     master bias image [ADU]
        """

        # Get all the paths of the files corresponding to these hashes
        paths = [ (self.db["BiasImages"].find_one({"_id": hash}))["path"] for hash in self.rawBiasHashes]

        # Get all the fits files corresponding to these hashes
        biases = tools.getImages(paths)

        # Use the image in the fits files, and use mean_combining to obtain the the master image
        masterBias = np.median(biases, axis=0)


        if outputFileName is not None:
            self.saveImageAndAddToDatabase(masterBias, outputFileName)
            print("Master bias image saved to fits file")

        # That's it!

        print("Block generated!")

        return masterBias



    def getHashOfOutputfile(self):
        """
        This function returns the hash id for the output file.
        This hash is made from the hash files that are used as input.

        Ouput:
            hash. string containing the output hash
        """
        hash = hashlib.sha256(bytes("".join(self.rawBiasHashes), 'utf-8')).hexdigest()
        return hash





if __name__ == "__main__":

    params = yaml.safe_load(open("params.yaml"))["rawBiasImage"]

    databaseName = "pipelineDatabase.txt"
    print("Creating a local database file with the name: ", databaseName)

    db = DatabaseFromLocalFile(databaseName)
    print(" ")

    # Master Bias Image
    raw_bias_path = params["path"]

    masterB = MasterBias(db, BiasImages=raw_bias_path)


    masterB.run(params["outpath"])

    db.save()



