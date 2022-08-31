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
    masterBlock.run()
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



    def run(self):
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


        # Add Master Bias Image to database
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
        self.outputPath = os.getcwd() + "/Data/ProcessedData/MasterDark/"
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
        self.input = self.createInputDirectory(input)
        self.outputPath = os.getcwd() + "/Data/ProcessedData/MasterFlat/"
        self.type = "Master Flat Image"



    def setCollection(self, input):
        return {"flat" : self.db["FlatImages"], "bias" : self.db["BiasImages"]}



    def checkInput(self, input):
        # We should have as input, one master Bias image, and more then one flat images.
        biasFrames = self.db["BiasImages"]
        flatFrames = self.db["FlatImages"]
        
        for hash in input:
            isBias = len([x for x in biasFrames.find({"_id" : hash})]) == 1
            isFlat = len([x for x in flatFrames.find({"_id" : hash})]) == 1

            if isBias:
                isMasterImage = np.all([x["type"] == "Master Bias Image" for x in biasFrames.find({"_id" : hash})])
                if not isMasterImage:
                    print("Bias image with hash {} is not Master Image".format(hash))
                    return False
            elif isFlat:
                isRawImage = np.all([x["type"] == "Raw Flat Image" for x in flatFrames.find({"_id" : hash})])
                if not isRawImage:
                    print("Flat image with hash {} is not Image".format(hash))
                    return False
            else:
                print("Image with hash {} is not bias or flat image".format(hash))
        return True



    def createInputDirectory(self, input):
        sortedInput = {"flat" : [], "bias" : []}
        for hash in input:
            if len([x for x in self.db["BiasImages"].find({"_id" : hash})]) == 1:
                sortedInput["bias"].append(hash)
            elif len([x for x in self.db["FlatImages"].find({"_id" : hash})]) == 1:
                sortedInput["flat"].append(hash)
        return sortedInput



    def make(self):

        # Get all the paths of the files corresponding to these hashes
        flatPaths = [([x["path"] for x in self.col["flat"].find({"_id" : hash})])[0] for hash in self.input["flat"]]
        biasPaths = [([x["path"] for x in self.col["bias"].find({"_id" : hash})])[0] for hash in self.input["bias"]]

        # Get all the files fits files corresponding to these hashes
        flats = tools.getImages(flatPaths)
        bias  = tools.getImages(biasPaths)

        # Use the image in the fits files, and use mean_combining to obtain the the master image
        MasterFlat = np.median(flats, axis=0) - np.median(bias, axis=0)

        # Add offset so that all the values in the MasterFlat are positive
        if np.min(MasterFlat) < 0:
                  MasterFlat = MasterFlat -np.min(MasterFlat)
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
        return { "science": self.db["ScienceImages"], "dark": self.db["DarkImages"], "bias": self.db["BiasImages"]}



    def checkInput(self, input):
        # We should have as input, one raw science image, one master bias frame and one master dark frame
        inputAmount = (len(input) == 3)
        if not inputAmount:
            print("Input is not in the correct format, we expect a list with 3 elements, instead {} were given.".format(len(input)))
            return False

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
        science = tools.getImages(sciencePath[0])
        dark = tools.getImages(darkPath[0])
        bias = tools.getImages(biasPath[0])

        # Use the image in the fits files, and generate the calibrated science frames.
        CalibratedScience = np.median(science - bias, axis=0)

        if np.min(CalibratedScience) < 0:
            CalibratedScience = CalibratedScience - np.min(CalibratedScience)
        return CalibratedScience



    def getFileName(self):
        return "calibrated_science_frames.fits"












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

    # Master Bias Image 
    raw_bias_hashes = ['f2753cdd3370f5596a6c574a5e837fb2837e04bc4b2bbb1dc4bd26f270849d45', '6d6187e691a99e49aa54fef06bfa289867ee36bb0c70ad1845f3a2ec1354b0f6',
                       '2ae6d82c0628c26af27553d7ba33c1121e32c755ef93c4fad13e66fb475c2127', '4e1b127813b2b8bd6a6f90285fbbc87cbdba697503606ea96eddbe5ec4affdbc', 
                       '3e947c38dee83746b746031d8dae57fa2d6a6f31c7fb0be31ad7f14b1f37b99b', '1902162ddc57a095c10a4571922cc3d76ead46eedeed3eefaac88c882097172a', 
                       '53fa9f81ffba0b3916bcb90603e49becb4e77eef0e73d6f8073132d8b585c703']
    masterB = MasterBias(raw_bias_hashes)
    masterB.run()



    # Master Dark Image
    raw_dark_hashes =  ['e547b0390ddcc6e0ec3b32bb85f2abf7c8f9f869edb45c068ec90e693883300c', '8161836c875b139b922fa3b0ca3dedd38a22474421846d6018ae1cdc0913cd86', 
                        '6846f7a8550ecf62f09ba4258097b5e2876ce5d70f831636b7144560958cdbc8']
    masterD = MasterDark(raw_dark_hashes)
    masterD.run()



    # Flat
    # 5flats <-> 1 master bias 
    raw_flat_hashes = ["9c43630b8c8865f9040ebf8938ece78b72849b4435a897e16291eed222801305", "e68a7f29ce87fb61d2aa58add658d18e08c78f679f3bbcc43071672c351fa6d6",
                      "cd1e1ffd95b22875a79163ab977e5f96bb0eab9d0f22574374176e1c5ed605ee", "74bb4c8de06386600d1f99f6bdc390aa84edd99c11fbc648c12d0a039f4dee47",
                       "6a6a2a048ea9c1c2fffd5fb3ca0f26df77866973300cf2de223d63dd9df32f93", "9b0c4e6cad3771c6f1a74186f2e7a3fa689a85a43f15b73067c23b6e8c64aa0d"]

    masterF = MasterFlat(raw_flat_hashes)
    masterF.run()



    # Calibrated Science Image
    raw_scienceimage_hashes = ["9b0c4e6cad3771c6f1a74186f2e7a3fa689a85a43f15b73067c23b6e8c64aa0d",
                               "edf0482526e97a6eb087d267b12af6fa7756c9709fb09c4f136a65ffc48ebaf1",
                               "3704c9e675cb7e438c4f2eb4e097b923f46e31d38f4466ef6d620cd643356735"]
    calibration = CalibratedScienceFrames(raw_scienceimage_hashes)
    calibration.run()
                                          


    




    
