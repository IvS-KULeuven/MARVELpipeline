# This is the parent class of the all pipeline components + components that generate master images
# TODO: Move the child components to their seperate file



import os
from pymongo import MongoClient
from astropy.io import fits
import tools
import numpy as np
import hashlib
from datetime import datetime
import matplotlib.pyplot as plt
from database import DatabaseFromLocalFile

client = MongoClient()






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
            inputHashes: imageType/hash (or path) of the fits file that is given as input.
        """

        if isinstance(database, DatabaseFromLocalFile):
            self.db = database
        else:
            if not "databaseMARVEL" in client.list_database_names():
                raise Exception("MARVEL database does not exist.")
            else:
                self.db = client["databaseMARVEL"]
        inputHashes = self.convertInputToHash(**inputHashes)




        if not (self.checkSanityOfInput(**inputHashes)):
            keysInDatabase = np.array([key in self.db.list_collection_names() for key in inputHashes.keys()])



            if not np.all(keysInDatabase):
                unrecognizedKeys = np.array(list(inputHashes.keys()))[~keysInDatabase]
                unrecognizedKeys = list(unrecognizedKeys)
                keys = str(unrecognizedKeys)[1:-1]
                raise Exception(f"Keyvalues:  {keys} not found in database")
            else:
                raise Exception("One of the input hashes/paths is not found in databse")
        self.inputHashes = inputHashes







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
                if type(inputHash[imageType]) == list:
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

            if type(value) == list:
                inputHashes[key] = [tools.convertPathToHash(item, self.db)
                                    if os.path.isfile(item) else item for item in value]
            else:
                if os.path.isfile(value):
                    inputHashes[key] = tools.convertPathToHash(value, self.db)

        return inputHashes







    def saveImageAndAddToDatabase(self, image, fileName):
        """
        Save the image and add it to the database
        """
        hash = self.getHashOfOutputfile()
        path = self.outputPath + fileName

        # Save Master Image as FITS file
        hdr = fits.Header()
        hdr["hash"] = hash
        hdr["path"] = path
        hdr["type"] = self.type

        hdu = fits.PrimaryHDU(image, header=hdr)
        hdu.writeto(path, overwrite=True)

        # Add Image to database
        currentTime = datetime.now()
        dict = {"_id" : hash, "path" : path, "type" : self.type, "date_created" : currentTime.strftime("%d/%m/%Y %H:%M:%S")}
        tools.addToDataBase(dict, self.db, overWrite = True)






class MasterBias(PipelineComponent):
    """
    Class that creates the master bias image. Such an image is obtained by taking the median of multiple
    raw bias images.
    """

    def __init__(self, database=None, **input):
        super().__init__(database, **input)
        input = self.inputHashes

        if self.checkSanityOfinputTypes(**input):

            self.outputPath = os.getcwd() + "/Data/ProcessedData/MasterBias/"
            self.type = "Master Bias Image"
            self.rawBiasHashes = input["BiasImages"]
        else:
            raise Exception("Error: The input hashes do not match the correct type: Aborting")
            exit(1)





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
            areRawImages = np.all([ self.db[types[0]].find_one({"_id": hash})["type"] == "Raw Bias Image"
                                    for hash in values[0]])
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








# class MasterDark(PipelineComponent):
#     """
#     Class that creates the master dark image. Such an image is obtained by taking the median of multiple
#     raw dark images.
#     """


#     def __init__(self, database=None, **darkHashes):
#         super().__init__(database, **darkHashes)
#         darkHashes = self.inputHashes

#         if self.checkSanityOfinputTypes(**darkHashes):

#             self.outputPath = os.getcwd() + "/Data/ProcessedData/MasterDark/"
#             self.type = "Master Dark Image"
#             self.rawDarkHashes = darkHashes["DarkImages"]
#         else:
#             raise Exception("Error: The input hashes do not match the correct type: Aborting")
#             exit(1)








#     def checkSanityOfinputTypes(self, **input):
#         """
#         This function is ran after we run checkSanityOfInputHashes. This function checks the the
#         input types that are given is able to generate a master bias output file.
#         """

#         types  = list(input.keys())
#         values = list(input.values())

#         # Check that the keys are of the right format. For a master dark these should only be DarkImages

#         keysAreCorrect = (len(types)) == 1 and types[0] == "DarkImages"

#         valuesAreCorrect = (len(values) == 1) and (len(values[0]) > 1)

#         if keysAreCorrect and valuesAreCorrect:
#             areRawImages = np.all([ self.db[types[0]].find_one({"_id": hash})["type"] == "Raw Dark Image"
#                                     for hash in values[0]])
#         else:
#             return False

#         # If we get here, the input is sane if the hashes correspond to raw Images

#         return areRawImages







#     def run(self, outputFileName=None):
#         """
#         We run through the alghorim to create the master dark image.

#         Input:
#             outputFileName: If None, nothings is saved. Otherwise, a string with the name fo the outputfile,
#                             incl. the extension ".fits".

#         Output:
#             masterDark:     master dark image [ADU]
#         """

#         # Get all the paths of the files corresponding to these hashes
#         paths = [ (self.db["DarkImages"].find_one({"_id": hash}))["path"] for hash in self.rawDarkHashes]

#         # Get all the fits files corresponding to these hashes
#         darks = tools.getImages(paths)

#         # Use the image in the fits files, and use mean_combining to obtain the the master image
#         masterDark = np.median(darks, axis=0)

#         if outputFileName is not None:
#             self.saveImageAndAddToDatabase(masterDark, outputFileName)
#             print("Master dark image saved to fits file")

#         # That's it!

#         print("Block generated!")
#         return masterDark






#     def getHashOfOutputfile(self):
#         """
#         This function returns the hash id for the output file.
#         This hash is made from the hash files that are used as input.

#         Ouput:
#             hash. string containing the output hash
#         """
#         hash = hashlib.sha256(bytes("".join(self.rawDarkHashes), 'utf-8')).hexdigest()
#         return hash







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






    def checkSanityOfInputTypes(self, **flatAndBiasHashes):
        """
        This function is ran after we run checkSanityOfInputHashes. This function checks the the
        input types that are given is able to generate a master flat output file.
        """

        types  = list(flatAndBiasHashes.keys())
        values = list(flatAndBiasHashes.values())

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
            self.outputPath = os.getcwd() + "/Data/ProcessedData/BiasCorrectedEtalon/"
            self.type       = "Bias Corrected Etalon Image"
            self.rawEtalonHash  = inputEtalonHashes["EtalonImages"]
            self.masterBiasHash = inputEtalonHashes["BiasImages"]
        else:
            raise Exception("Error: The input hashes do not match the correct type: Aborting")
            exit(1)








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

    db = DatabaseFromLocalFile("pipelineDatabase.txt")
    print("")

    # Master Bias Image
    raw_bias_hashes = ['f2753cdd3370f5596a6c574a5e837fb2837e04bc4b2bbb1dc4bd26f270849d45',
                       '6d6187e691a99e49aa54fef06bfa289867ee36bb0c70ad1845f3a2ec1354b0f6',
                       '2ae6d82c0628c26af27553d7ba33c1121e32c755ef93c4fad13e66fb475c2127',
                       '4e1b127813b2b8bd6a6f90285fbbc87cbdba697503606ea96eddbe5ec4affdbc',
                       '3e947c38dee83746b746031d8dae57fa2d6a6f31c7fb0be31ad7f14b1f37b99b',
                       '1902162ddc57a095c10a4571922cc3d76ead46eedeed3eefaac88c882097172a',
                       '53fa9f81ffba0b3916bcb90603e49becb4e77eef0e73d6f8073132d8b585c703']

    masterB1 = MasterBias(db, BiasImages=raw_bias_hashes)
    masterB2 = MasterBias(BiasImages=raw_bias_hashes)
    masterB1.run("testsFBias.fits")
    masterB2.run("testsDBias.fits")
    print("")


    # # Master Dark Image
    # raw_dark_hashes =  ['6128a1e361aca3d5b17e366511efc75b0aeffbada070cbc4a3aebdd1cb1d66db',
    #                     '649d0b67d7be70ef286d86c9629bd017716cd3667fde46beeab35dcd27a98f0c',
    #                     'f53e4b7837347cdcddf0bf41a1cd5ac40f7594c4561ac8b71b08fb4da541f1f5']
    # masterD1 = MasterDark(db, DarkImages=raw_dark_hashes)
    # masterD2 = MasterDark(DarkImages=raw_dark_hashes)
    # masterD1.run("testFDark.fits")
    # masterD2.run("testDDark.fits")
    # print("")



    # Flat
    raw_flat_hashes = ["9c43630b8c8865f9040ebf8938ece78b72849b4435a897e16291eed222801305",
                       "e68a7f29ce87fb61d2aa58add658d18e08c78f679f3bbcc43071672c351fa6d6",
                       "cd1e1ffd95b22875a79163ab977e5f96bb0eab9d0f22574374176e1c5ed605ee",
                       "74bb4c8de06386600d1f99f6bdc390aa84edd99c11fbc648c12d0a039f4dee47",
                       "6a6a2a048ea9c1c2fffd5fb3ca0f26df77866973300cf2de223d63dd9df32f93"]
    master_bias_pathF = "Data/ProcessedData/MasterBias/testsFBias.fits"
    master_bias_pathD = "Data/ProcessedData/MasterBias/testsDBias.fits"
    masterF1 = MasterFlat(db, FlatImages=raw_flat_hashes, BiasImages=master_bias_pathF)
    masterF2 = MasterFlat(FlatImages=raw_flat_hashes, BiasImages=master_bias_pathD)

    masterF1.run("testsFFlat.fits")
    masterF2.run("testsDFlat.fits")


    # bias corrected Science Image
    rawScienceHash = "d99dff18a15ab83b51af4ecdca9dc99c2069bdc17eb6a5d1cdb67e7e86c92e4a"

    master_bias_pathF = "Data/ProcessedData/MasterBias/testsFBias.fits"
    master_bias_pathD = "Data/ProcessedData/MasterBias/testsDBias.fits"


    print(" ")
    calibration1 = BiasCorrectedScienceFrames(db, ScienceImages=rawScienceHash,
                                              BiasImages=master_bias_pathF)
    calibration2 = BiasCorrectedScienceFrames(ScienceImages=rawScienceHash,
                                           BiasImages=master_bias_pathD)

    calibration1.run("testFScience.fits")
    calibration2.run("testDScience.fits")


    # Calibrated Etalon Image
    raw_etalon_hash  = "e0ac021d19ce5520d0ba92df1d5aadf6541f35f76a84121576828287937ca508"

    master_bias_pathF = "Data/ProcessedData/MasterBias/testsFBias.fits"
    master_bias_pathD = "Data/ProcessedData/MasterBias/testsDBias.fits"

    print(" ")
    calibratedEtalon1 = BiasCorrectedEtalonImage(db, EtalonImages=raw_etalon_hash,
                                                 BiasImages=master_bias_pathF)
    calibratedEtalon2 = BiasCorrectedEtalonImage(EtalonImages=raw_etalon_hash,
                                              BiasImages=master_bias_pathD)
    calibratedEtalon1.run("testFEtalon.fits")
    calibratedEtalon2.run("testDEtalon.fits")


    db.save()
