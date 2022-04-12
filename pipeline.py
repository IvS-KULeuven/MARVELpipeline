# This is the parent class of the all pipeline components + components that generate master images 
# TODO: Move the child components to their seperate file  



import os
from pymongo import MongoClient
from astropy.io import fits
import tools
import numpy as np
import hashlib

client = MongoClient()



"""
MARVEL pipeline building blocks

The PipelineComponent represents one isolated block of the pipeline. Every end product of a block has a hash code associated 
to it, which is also added to the MARVEL database. A block is created using a list a hashes that represents other block-outputs
 as input. PipelineComponents are isolated blocks in that they do not depend on other components. This implies that two blocks 
that have the same input should also have the same output (within the same version of the pipeline).

Once we have created a block the output is generated and added to the MARVEL database by running the runComponent() method on a PipelineComponent instance. 

Example:
    masterBlock = MasterBias(hash_list)
    masterBlock.runComponent()
"""


class PipelineComponent():

    def __init__(self, input):

        if not "databaseMARVEL" in client.list_database_names():
            raise Exception("MARVEL database does not exist.")
        else:
            self.db = client["databaseMARVEL"]
            if (self.checkInput(input)):
                self.input = input
            else:
                raise Exception("Input is not correct format")

    def checkInput(self, input):
        return True

    def runComponent(self):
        img = self.make()
        self.saveImage(img)
        print("Block Generated!")

    def make(self):
        ...

    def saveImage(self, image):
        hash = hashlib.sha256(bytes("".join(self.input), 'utf-8')).hexdigest()
        path = self.outputPath + self.getFileName()

        # Save Master Image as FITS file
        hdr = fits.Header()
        hdr["hash"] = hash
        hdr["path"] = path
        hdr["type"] = self.type

        hdu = fits.PrimaryHDU(image, header=hdr)
        hdu.writeto(path, overwrite=True)


        # Add Mast Bias Image to database
        dict = {"_id" : hash, "path" : path, "type" : self.type}
        tools.addToDataBase(dict, overWrite = True)





class MasterBias(PipelineComponent):

    def __init__(self, input):
        super().__init__(input)
        self.col = self.setCollection(input)
        self.outputPath = os.getcwd() + "/Data/ProcessedData/MasterBias/"
        self.type = "Master Bias Image"


    def setCollection(self, input):
        return self.db["BiasImages"]

    def checkInput(self, input):
        isCorrect = []
        collection = self.db["BiasImages"]
        for hash in input:
            instances = collection.find({"_id" : hash})
            isCorrect.append(np.all([( x["type"] == "Raw Bias Image") for x in instances]))
        return np.all(isCorrect)


    def make(self):
        # Get all the paths of the files corresponding to these hashes
        paths = [[x["path"] for x in self.col.find({"_id" : hash})] for hash in self.input]
        paths = [ x[0] for x in paths]

        # Get all the files fits files corresponding to these hashes
        biases = tools.getImages(paths)
 
        # Use the image in the fits files, and use mean_combining to obtain the the master image
        MasterBias = np.median(biases, axis=0)
        return MasterBias

    def getFileName(self):
        return "master_bias.fits"



        











class MasterDark(PipelineComponent):


    def __init__(self, input):
        super().__init__(input)
        self.col = self.setCollection(input)
        self.outputPath = os.getcwd() + "/Data/ProcessedData/MaterDark/"
        self.type = "Master Dark Image"


    def setCollection(self, input):
        return self.db["DarkImages"]

    def checkInput(self, input):
        isCorrect = []
        collection = self.db["DarkImages"]
        for hash in input:
            instances = collection.find({"_id" : hash})
            isCorrect.append(np.all([( x["type"] == "Raw Dark Image") for x in instances]))

        return np.all(isCorrect)
    
    def make(self):
        # Get all the paths of the files corresponding to these hashes
        paths = [[x["path"] for x in self.col.find({"_id" : hash})] for hash in self.input]
        paths = [ x[0] for x in paths]

        # Get all the files fits files corresponding to these hashes
        darks = tools.getImages(paths)
 
        # Use the image in the fits files, and use mean_combining to obtain the the master image
        MasterDark = np.median(darks, axis=0)
        return MasterDark

    def getFileName(self):
        return "master_dark.fits"





class MasterFlat(PipelineComponent):

    def __init__(self, input):
        super().__init__(input)
        self.col = self.setCollection(input)
        self.outputPath = os.getcwd() + "/Data/ProcessedData/MasterFlat/"
        self.type = "Master Flat Image"


    def setCollection(self, input):
        return self.db["FlatImages"]


    def checkInput(self, input):
        isCorrect = []
        collection = self.db["FlatImages"]
        for hash in input:
            instances = collection.find({"_id" : hash})
            isCorrect.append(np.all([( x["type"] == "Raw Flat Image") for x in instances]))

        return np.all(isCorrect)

    def make(self):
        # Get all the paths of the files corresponding to these hashes
        paths = [[x["path"] for x in self.col.find({"_id" : hash})] for hash in self.input]
        paths = [ x[0] for x in paths]

        # Get all the files fits files corresponding to these hashes
        flats = tools.getImages(paths)
 
        # Use the image in the fits files, and use mean_combining to obtain the the master image
        MasterFlat = np.median(flats, axis=0)
        return MasterFlat


    def getFileName(self):
        return "master_flat.fits"








class CalibratedScienceFrames(PipelineComponent):

    def __init__(self, input):
        super().__init__(input)
        self.col = self.setCollection(input)
        self.input = self.createInputDirectory(input)
        self.outputPath = os.getcwd() + "/Data/ProcessedData/CalibratedScience/"
        self.type = "Calibrated Science Image"

    def setCollection(self, input):
        return { "science": super().db["ScienceImages"], "dark": super().db["DarkImages"], "bias": super().db["BiasImages"]}

    def checkInput(self, input):
        # We should have as input, one raw science image, one master bias frame and one master dark frame
        inputAmount = (len(input) == 3)
        if not inputAmount:
            print("Input is not in the correct format, we expect a list with 3 elements, instead {} were given.".format(len(input)))
            return 

        scienceFrames = self.db["ScienceImages"]
        darkFrames    = self.db["DarkImages"]
        biasFrames    = self.db["BiasImages"]

        for hash in input:
            isScience = len([ x for x in scienceFrames.find({"_id" : hash})]) == 1
            isDark = len([ x for x in darkFrames.find({"_id" : hash})]) == 1
            isBias = len([ x for x in biasFrames.find({"_id" : hash})]) == 1

            if isScience:
                isRawImage = np.all([ x["type"] == "Raw Science Image" for x in scienceFrames.find({"_id" : hash})])
                if not isRawImage:
                    print("Science image with hash {} is not Raw Image".format(hash))
                    return False
            elif isDark:
                isRawImage = np.all([ x["type"] == "Master Dark Image" for x in darkFrames.find({"_id" : hash})])
                if not isRawImage:
                    print("Dark image with hash {} is not Master Image".format(hash))
                    return False
            elif isBias:
                isRawImage = np.all([ x["type"] == "Master Bias Image" for x in biasFrames.find({"_id" : hash})])
                if not isRawImage:
                    print("Bias image with hash {} is not Master Image".format(hash))
                    return False
            else:
                # Is not science, dark or bias image
                return False
        return True


    def createInputDirectory(self, input):
        sortedInput = {"science" : [], "bias" : [], "dark" : []}
        for hash in input:
            if len([x for x in self.db["ScienceImages"].find({"_id" : hash})]) == 1:
                sortedInput["science"].append(hash)
            elif len([x for x in self.db["BiasImages"].find({"_id" : hash})]) == 1:
                sortedInput["bias"].append(hash)
            elif len([x for x in self.db["DarkImages"].find({"_id" : hash})]) == 1:
                sortedInput["dark"].append(hash)
        return sortedInput



    def make(self):
        # Get all the paths of the files corresponding to these hashes
        sciencePath = [ [x["path"] for x in self.col["science"].find({"_id" : hash})] for hash in self.input["science"] ]
        darkPath    = [ [x["path"] for x in self.col["dark"].find({"_id" : hash})] for hash in self.input["dark"] ]
        biasPath    = [ [x["path"] for x in self.col["bias"].find({"_id" : hash})] for hash in self.input["bias"] ]

        # Get all the files fits files corresponding to these hashes
        science = tools.getImages(sciencePath)[0]
        dark = tools.getImages(darkPath)[0]
        bias = tools.getImages(biasPath)[0]

        # Use the image in the fits files, and generate the calibrated science frames.
        return science - dark - bias

    def getFileName(self):
        return "calibrated_science_frames.yaml"












"""
Whenever we want to execute a block from the pipeline.
1. We give input as a list of hashes
2. We use the DB in order to check that the format of these files is correct
3. From the DB we obtain the path to the files
4. Every "type" Block will be its own subclass of pipelineBlock
==============================================================================
=> Check that master dark works well
=> Master Dark and scaled to exposure time of the SF
=> Check that calibrated master frames works
"""


if __name__ == "__main__":
    # hash_list = ['900515b630b9787e4ab5032242228ac21dda537355de1a790270a7871a407652', 'a4880b4ce9510679210561bff9432ac6c42b243d995c7d4adf23902d77757662', '8afd899dc50f8532c9a44d8637ced940eaa702e30e3ec7c8e3adde5119b7c07b', '9b5f55279418cf54c90432d45319b6e7d5c74d9b03ccd444483f4ac0893c4b83', '21f289581264659eaeb76ec2163bc521575a444b1a00628a870302c847242053', '80baf7fd974a26b3d3f66627b4e150aa814e9bdd92b190dc2715d03257e60633', 'a947d1cdd89cfc7241bf0fa528357bc7a6f34e48fd9d1820d059ec270abcff49', 'fbfd5d76908a512a07cce993688c583b4b68729d5baba14b6107998619f4f16b', '652c6fcce89348d34195f06fddce1d35e63d3fdd5c633637d80b13d6f4c6a49d', '12b9a41a8c577a369b1513035426edf5372210c96051a4b18a7c868add87a30c']
    # masterB = MasterBias(hash_list)
    # masterB.runComponent()
    
    # calibration = CalibratedScienceFrames(["5e323082bef426f425d01c05436f68f09a3d19158590a410cf975e3670b38663", "642393b786ede43436654b9fd477d7f72cbd30b389edb6cafc5ac5e45fc94554", "84762ec2e7393bfffa3ea505692d33b599a5dfa2cda1efa3f791726c20c2df5d"])
    # calibration.runComponent()
                                          


    # Flat
    hash_list = ["1440a24e48f0e043da26dfbb26595f263aee42ecfeb13b540c38b6418a8351d3", "40e031729d8131b9605b7feecd8f96bc95febcb39faa40f5b427d814de726eff", 
                 "9fb0578e82514a3b73b7fef554363ea16b1bf8aed1347511b4475fe4f9d79b6e", "7b3e72779ae2fc29b22761cddf4cf3d2b4d2803882e87fdfc7fa99c506abcf88", 
                 "50359821613eca4a88daed91be15356011a87d98a506bbe0d072970a4f764534"]
    masterF = MasterFlat(hash_list)
    masterF.runComponent()
    





    
