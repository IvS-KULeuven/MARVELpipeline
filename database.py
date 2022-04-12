# This function creates/resets a database from the raw input images.
# 



from pymongo import MongoClient
import hashlib
import os


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
        db["BiasImages"].drop()
        db["DarkImages"].drop()
        db["ScienceImages"].drop()
        db["FlatImages"].drop()

    pathToRaw  = os.getcwd() + "/Data/RawData/"
    pathToPros = os.getcwd() + "/Data/ProcessedData/Data/"


    # Add the Bias images to the DataBase

    biasImagePath  = pathToRaw + "CalibrationImages/Bias"
    biasImages     = os.listdir(biasImagePath)
    if len(biasImages) != 0:
        biasCollection = db["BiasImages"]
        biases = [addBiasImage(image, biasImagePath + "/" + image, biasCollection) for image in biasImages]
        biasCollection.insert_many(biases)
    

    # Add the Dark images to the DataBase
    
    darkImagePath  = pathToRaw + "/CalibrationImages/Dark"
    darkImages     = os.listdir(darkImagePath)
    if len(darkImages) != 0:
        darkCollection = db["DarkImages"]
        darks = [addDarkImage(image, darkImagePath + "/" + image , darkCollection) for image in darkImages]
        darkCollection.insert_many(darks)

    # Add the Flat images to the Database

    flatImagePath  = pathToRaw + "/CalibrationImages/Flat"
    flatImages     = os.listdir(flatImagePath)
    flatCollection = db["FlatImages"]
    flats = [addFlatImage(image, flatImagePath + "/" + image, flatCollection) for image in flatImages]
    flatCollection.insert_many(flats)


    # Add the Science images to the DataBase

    scienceImagePath = pathToRaw + "/ScienceFrames"
    scienceImages = os.listdir(scienceImagePath)
    if len(scienceImages) != 0:
        scienceCollection = db["ScienceImages"]
        sciences = [addScienceImage(image, scienceImagePath + "/" + image, scienceCollection) for image in scienceImages]
        scienceCollection.insert_many(sciences)



def addBiasImage(image, path, collection):
    # Should think about what information form the file should be included in the hasInput
    hashInput = image + "Raw Bias Image"
    hash = hashlib.sha256(bytes(hashInput, 'utf-8')).hexdigest()
    
    return {"_id" : hash, "path" : path, "type" : "Raw Bias Image"}


def addDarkImage(image, path, collection):
    # Should think about what information form the file should be included in the hasInput
    hashInput = image + "Raw Dark Image"
    hash = hashlib.sha256(bytes(hashInput, 'utf-8')).hexdigest()
    return {"_id": hash, "path": path, "type": "Raw Dark Image"}


def addFlatImage(image, path, collection):
    # Should think about what information form the file should be included in the hasInput
    hashInput = image + "Raw Flat Image"
    hash = hashlib.sha256(bytes(hashInput, 'utf-8')).hexdigest()
    return {"_id": hash, "path": path, "type": "Raw Flat Image"}


def addScienceImage(image, path, collection):

    hashInput = image + "Raw Science Image"
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
