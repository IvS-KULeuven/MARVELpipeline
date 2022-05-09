# This function creates/resets a database from the raw input images.
 

from pymongo import MongoClient
import hashlib
import os
from glob import glob



def createDataBase(makeIfExist=False):

    client = MongoClient()

    # Check if database already exist

    if ("databaseMARVEL" in client.list_database_names() and not makeIfExist):
        print("There already exist a Database for MARVEL. To explitly remake the database run:")
        print("createDataBase(makeIfExist=True)")
        return

    db = client["databaseMARVEL"]

    # If database already exist, we make sure it is empty before we start 

    if ("databaseMARVEL" in client.list_database_names() and makeIfExist):
        print("Clearing existing database")
        db["BiasImages"].drop()
        db["DarkImages"].drop()
        db["ScienceImages"].drop()
        db["FlatImages"].drop()

    pathToRaw  = os.getcwd() + "/Data/RawData/"
    pathToPros = os.getcwd() + "/Data/ProcessedData/Data/"

    # Add the Bias images to the DataBase

    biasImageDirectory = pathToRaw + "CalibrationImages/Bias"
    biasImagePaths = glob(biasImageDirectory + "/*.fits")
    print("Found {0} bias frames to be inserted in MongoDB".format(len(biasImagePaths)))
    if len(biasImagePaths) != 0:
        biasImageBaseNames = [os.path.splitext(os.path.basename(imagePath))[0] for imagePath in biasImagePaths]
        biasCollection = db["BiasImages"]
        biasImages = [addBiasImage(biasImageBaseNames[n], biasImagePaths[n], biasCollection) for n in range(len(biasImagePaths))]
        biasCollection.insert_many(biasImages)


    # Add the Dark images to the DataBase
    
    darkImageDirectory = pathToRaw + "CalibrationImages/Dark"
    darkImagePaths = glob(darkImageDirectory + "/*.fits")
    print("Found {0} dark frames to be inserted in MongoDB".format(len(darkImagePaths)))
    if len(darkImagePaths) != 0:
        darkImageBaseNames = [os.path.splitext(os.path.basename(imagePath))[0] for imagePath in darkImagePaths]
        darkCollection = db["DarkImages"]
        darkImages = [addDarkImage(darkImageBaseNames[n], darkImagePaths[n], darkCollection) for n in range(len(darkImagePaths))]
        darkCollection.insert_many(darkImages)

    # Add the Flat images to the Database

    flatImageDirectory = pathToRaw + "CalibrationImages/Flat"
    flatImagePaths = glob(flatImageDirectory + "/*.fits")
    print("Found {0} flat frames to be inserted in MongoDB".format(len(flatImagePaths)))
    if len(flatImagePaths) != 0:
        flatImageBaseNames = [os.path.splitext(os.path.basename(imagePath))[0] for imagePath in flatImagePaths]
        flatCollection = db["FlatImages"]
        flatImages = [addFlatImage(flatImageBaseNames[n], flatImagePaths[n], flatCollection) for n in range(len(flatImagePaths))]
        flatCollection.insert_many(flatImages)


    # Add the Science images to the DataBase

    scienceImageDirectory = pathToRaw + "ScienceFrames"
    scienceImagePaths = glob(scienceImageDirectory + "/*.fits")
    print("Found {0} science frames to be inserted in MongoDB".format(len(scienceImagePaths)))
    if len(scienceImagePaths) != 0:
        scienceImageBaseNames = [os.path.splitext(os.path.basename(imagePath))[0] for imagePath in scienceImagePaths]
        scienceCollection = db["ScienceImages"]
        scienceImages = [addScienceImage(scienceImageBaseNames[n], scienceImagePaths[n], scienceCollection) for n in range(len(scienceImagePaths))]
        scienceCollection.insert_many(scienceImages)



def addBiasImage(imageName, path, collection):
    """
    Does something for which a doc string is required
    """
    hashInput = imageName + "Raw Bias Image"
    hash = hashlib.sha256(bytes(hashInput, 'utf-8')).hexdigest()
    return {"_id" : hash, "path" : path, "type" : "Raw Bias Image"}



def addDarkImage(imageName, path, collection):
    hashInput = imageName + "Raw Dark Image"
    hash = hashlib.sha256(bytes(hashInput, 'utf-8')).hexdigest()
    return {"_id": hash, "path": path, "type": "Raw Dark Image"}



def addFlatImage(imageName, path, collection):
    hashInput = imageName + "Raw Flat Image"
    hash = hashlib.sha256(bytes(hashInput, 'utf-8')).hexdigest()
    return {"_id": hash, "path": path, "type": "Raw Flat Image"}



def addScienceImage(imageName, path, collection):
    hashInput = imageName + "Raw Science Image"
    hash = hashlib.sha256(bytes(hashInput, 'utf-8')).hexdigest()
    return {"_id": hash, "path": path, "type": "Raw Science Image"}

    
    






if __name__ == "__main__":

    createDataBase(makeIfExist=True)
    
    
    




"""
NOTES: 
1. We should check what happens when we add things that are already in the collections. Are the duplicated or not? 
2. We should generate a _uniqe_ hash code, that is reproducible (stored in _id)
3. What information should Dark and Bias images have: - date 
"""
