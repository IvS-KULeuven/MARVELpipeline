# This is the parent class of the all pipeline components + components that generate master images
# TODO: Move the child components to their seperate file



import os
#from pymongo import MongoClient
from astropy.io import fits
import tools
import numpy as np
import hashlib
from datetime import datetime

from database import DatabaseFromLocalFile

#client = MongoClient()






class PipelineComponent():
    """
    MARVEL pipeline building components

    The PipelineComponent represents one isolated component of the pipeline.
    Every end product of a component has a hash code associated to it, which
    is also added to the MARVEL database. A component is created using a list
    a hashes that represents other component-outputs as input.
    PipelineComponents are isolated components in that they do not depend on
    other components. This implies that two components that have the same
    input should also have the same output (within the same version of the pipeline).

    Once we have created a component the output is generated and added to the
    MARVEL database by running the runComponent() method on a PipelineComponent instance.

    Example:
        masterComponent = MasterBias(BiasImages=hash_list)
        masterComponent.run()
    """

    def __init__(self, database=None, **inputHashes):
        """
        Initializes the component of the pipeline. This functions is called after a component
        is generated. This checks that the input hashes are sensible.

        Input:
            database:    database object that is used to store and access information form the fits files.
            inputHashes: imageType/hash (or path) of the FITS file that is given as input.
        """

        # if isinstance(database, DatabaseFromLocalFile):
        #     self.db = database
        # else:
        #     if not "databaseMARVEL" in client.list_database_names():
        #         raise Exception("MARVEL database does not exist.")
        #     else:
        #         self.db = client["databaseMARVEL"]
        self.db = database
        inputHashes = self.convertInputToHash(**inputHashes)


        if not (self.checkSanityOfInput(**inputHashes)):
            keysInDatabase = np.array([key in self.db.list_collection_names() for key in inputHashes.keys()])

            if not np.all(keysInDatabase):
                unrecognizedKeys = np.array(list(inputHashes.keys()))[~keysInDatabase]
                unrecognizedKeys = list(unrecognizedKeys)
                keys = str(unrecognizedKeys)[1:-1]
                raise Exception(f"Keyvalues: {keys} not found in database")
            else:
                raise Exception("One of the input hashes/paths is not found in database")
        self.inputHashes = inputHashes

        if not os.path.isdir(os.getcwd() + "/Data/ProcessedData"):
            os.mkdir(os.getcwd() + "/Data/ProcessedData")







    def checkSanityOfInput(self, **inputHash):
        """
        Check if the image hashes indeed corresponds to an image of the correct format (given by
        imageType).

        Input:
             **inputHash: imageType/hash(es) of the images given as input.

        Output:
            isSane: boolean. True if both input hashes correspond to images with the correct type, False if at
                             least one of these hashes does not.
        """

        for imageType in inputHash:
            if imageType not in self.db.list_collection_names():
                return False
            else:
                # If there are multiple hashes corresponding to one imagetype

                if isinstance(inputHash[imageType], list):
                    for hash in inputHash[imageType]:
                        image = self.db[imageType].find_one({"_id": hash})
                        if image is None:
                            return False

                # If there is one hash corresponding to one imageType

                else:
                    image = self.db[imageType].find_one({"_id":  inputHash[imageType]})
                    if image is None:
                        return False

        # If we get here, the hashes were found at the correct location

        return True






    def convertInputToHash(self, **inputHashes):
        """
        Checks the input. If the inputHashes contain paths to files (instead of image hashes),
        convert the path to the corresponding hash.

        Input:
            inputHashes: inputType/hash (or path)

        Output:
            inputHashes: inputType/hash
        """

        for key, value in inputHashes.items():
            if isinstance(value, list):
                inputHashes[key] = [tools.convertPathToHash(item, self.db)
                                    if os.path.isfile(item) else item for item in value]

            else:
                if os.path.isfile(value):
                    inputHashes[key] = tools.convertPathToHash(value, self.db)

        return inputHashes







    def saveMultipleImagesAndAddToDatabase(self, images, fileNames, imageHashes, **keywords):
        """
        Save multiple images and add them all to the database.

        input:
            image: list of np.arrays of the images we want to save to fits files
            fileNames: paths of the fits files
            keyswords: dictionary with option parameters we want to add to the fits header
        """
        for file, (path, hash) in zip(images, zip(fileNames, imageHashes)):
            self.saveImageAndAddToDatabase(file, path, imageHash=hash, **keywords)
        



    def saveImageAndAddToDatabase(self, image, fileName, imageHash=None, **keywords):
        """
        Save the image and add it to the database.

        input:
            image: np.array of the image that we want to save in fits files
            filname: path of the fits file
            keywords: dictionary with optional parameters we want to add to the fits header
        """
        hash = self.getHashOfOutputfile(imageHash)
        path = self.outputPath + "/" + fileName
        num_row, num_col = image.shape

        # Save Master Image as FITS file
        hdr = fits.Header()
        hdr["hash"] = hash
        hdr["path"] = path
        hdr["type"] = self.type
        hdr["rows"] = num_row
        hdr["cols"] = num_col

        for key, value in keywords.items():
            hdr[key] = value

        
        hdu = fits.PrimaryHDU(image, header=hdr)
        hdu.writeto(path, overwrite=True)
        
        # Add Image to database
        currentTime = datetime.now()
        dict = {"_id" : hash, "path" : path, "type" : self.type, "date_created" : currentTime.strftime("%d/%m/%Y %H:%M:%S")}
        tools.addToDataBase(dict, self.db, overWrite = True)










































class BiasCorrectedEtalonImage(PipelineComponent):
    """
    The component calculates the Bias Corrected Etalon Image. This  is caluculated
    by substracting a master bias image from the etalon calibration image. Finaly,
    if their are negative values in the output file we subtstract from this image
    its minimum value so that the new minimum is 0.
    """

    def __init__(self, database=None, **inputEtalonHashes):
        """
        Initialize the bias corrected etalon image.
        """

        super().__init__(database, **inputEtalonHashes)
        inputEtalonHashes = self.inputHashes

        if self.checkSanityOfInputTypes(**inputEtalonHashes):
            self.outputPath = os.getcwd()
            self.type       = "Bias Corrected Etalon Image"
            self.rawEtalonHash  = inputEtalonHashes["EtalonImages"]
            self.masterBiasHash = inputEtalonHashes["BiasImages"]
        else:
            raise Exception("Error: The input hashes do not match the correct type: Aborting")
            exit(1)

        if not os.path.isdir(self.outputPath):
            os.mkdir(self.outputPath)






    def checkSanityOfInputTypes(self, **inputEtalonHashes):
        """
        This function is ran after we run checkSanityOfInputHashes. This function checks the the
        input types that are given is able to generate a bias corrected etalon output file.
        """
        types  = list(inputEtalonHashes.keys())
        values = list(inputEtalonHashes.values())

        # Check that the keys are of the right format. For a bias corrected etalon image, these should be etalon
        # images, dark images and bias images.

        correctKeysFormat = ["EtalonImages", "BiasImages"]
        keysAreCorrect = (len(types) == 2) and np.all([key in types for key in correctKeysFormat])

        if keysAreCorrect:
            valuesAreCorrectType = np.all([ isinstance(value, str) for value in values])
        else:
            return False

        if valuesAreCorrectType:

            isMasterBiasImage = self.db["BiasImages"].find_one({"_id": inputEtalonHashes["BiasImages"]})["type"] == "Master Bias Image"
            isRawEtalonImage = self.db["EtalonImages"].find_one({"_id": inputEtalonHashes["EtalonImages"]})["type"] == "Raw Etalon Image"
        else:
            return False

        return isMasterBiasImage and isRawEtalonImage










    def run(self, outputFileName=None):
        """
        Runs through the steps to get the bias corrected etalon image.

        Input:
            outputFileName: If None, nothing is saved. Otherwise, a string with the name of the outputfile,
                            incl. the extension ".fits".

        Ouput:
            array with the flux values of the bias corrected etalon image.
        """

        # Get all the paths of the files corresponding to the input hashes

        etalonPath = self.db["EtalonImages"].find_one({"_id": self.rawEtalonHash})["path"]
        biasPath   = self.db["BiasImages"].find_one({"_id": self.masterBiasHash})["path"]

        # Get all the images that correspond to these paths

        etalon = tools.getImage(etalonPath)
        bias   = tools.getImage(biasPath)

        # Use these images to get the bias corrected etalon image.

        calibratedEtalon = etalon - bias
        if np.min(calibratedEtalon) < 0:
            calibratedEtalon = calibratedEtalon - np.min(calibratedEtalon)

        # If required, save to FITS
        if outputFileName is not None:
            self.saveImageAndAddToDatabase(calibratedEtalon, outputFileName)
            print("Bias corrected etalon image saved to fits file")

        # That's it!

        print("Block generated!")
        return calibratedEtalon





    def getHashOfOutputfile(self):
        """
        The function returns the hash id for the output file.
        This hash is made from the hash files that are used as input.

        Ouput:
           hash. string containing the output hash
        """
        combinedHashes = self.masterBiasHash + self.rawEtalonHash
        hash = hashlib.sha256(bytes(combinedHashes, 'utf-8')).hexdigest()
        return hash

























"""
Whenever we want to execute a component from the pipeline.
1. We give input as a list of hashes
2. We use the DB in order to check that the format of these files is correct
3. From the DB we obtain the path to the files
4. Every "type" Component will be its own subclass of pipelineComponent
==============================================================================
=> Check that master dark works well
=> Master Dark and scaled to exposure time of the SF
=> Check that calibrated master frames works
"""


if __name__ == "__main__":



    databaseName = "pipelineDatabase.txt"
    print("Creating a local database file with the name: ", databaseName)

    db = DatabaseFromLocalFile(databaseName)





    # # Calibrated Etalon Image
    # raw_etalon_hash  = "cd883d3e471c1be5176a74565f1d1b7a29fbe8a52c6e8ca1556a44cdac25c3b1"

    # master_bias_pathF = "Data/ProcessedData/MasterBias/testsFBias.fits"
    # master_bias_pathD = "Data/ProcessedData/MasterBias/testsDBias.fits"

    # print(" ")
    # calibratedEtalon1 = CalibratedEtalonImage(db, EtalonImages=raw_etalon_hash, BiasImages=master_bias_pathF)
    # calibratedEtalon2 = CalibratedEtalonImage(EtalonImages=raw_etalon_hash, BiasImages=master_bias_pathD)
    # calibratedEtalon1.run("testFEtalon.fits")
    # calibratedEtalon2.run("testDEtalon.fits")


    db.save()
