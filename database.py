
from pymongo import MongoClient
import hashlib
import os
from glob import glob
from datetime import datetime


def createHash(imageName, imageType):
    """
    Create a SHA256 hash of an image based on the filename and the image type.
    
    Input:
        imageName: base filename of the image, without extension. E.g. "science_0004"
        imageType: string containing the image type. Supported image types are
                   "Raw Bias Image", "Raw Dark Image", "Raw Flat Image", "Raw Etalon Image", 
                   and "Raw Science Image".

    Output:
        hash: string containing the SHA256 hash
    """
    supportedImageTypes = ["Raw Bias Image", "Raw Dark Image", "Raw Flat Image", "Raw Etalon Image", "Raw Science Image"]
    if imageType not in supportedImageTypes:
        print(f"Error: '{imageType}' is not one of the supported image types")
        exit(1)

    hashInput = imageName + imageType
    hash = hashlib.sha256(bytes(hashInput, 'utf-8')).hexdigest()
    return hash





def initializeDataBaseWithRawImages(clearIfExist=False):
    """ 
    Create and initialize the MongoDB database with the information on the raw CCD fits images.
    Only information is stored, not the actual fits files.
    It is assumed that a MongoDB daemon is running, e.g. execute "mongod --dbpath ~/MongoDB" on the prompt.

    Input:
        clearIfExist: clear all existing information if the database was already previously initialized.

    Output:
        None
    """
    client = MongoClient()

    # Check if database already exist

    if ("databaseMARVEL" in client.list_database_names() and not clearIfExist):
        print("There already exist a Database for MARVEL. To explitly remake the database run:")
        print("createDataBase(clearIfExist=True)")
        return

    db = client["databaseMARVEL"]

    # If database already exist, make sure it is emptied before we start 

    if ("databaseMARVEL" in client.list_database_names() and clearIfExist):
        print("Clearing existing database")
        db["BiasImages"].drop()
        db["DarkImages"].drop()
        db["ScienceImages"].drop()
        db["FlatImages"].drop()
        db["EtalonImages"].drop()

    pathToRaw  = os.getcwd() + "/Data/RawData/"
    print(f"Searching for raw images in {pathToRaw}")

    for imageType in ["Bias", "Dark", "Flat", "Etalon", "Science"]:
        if imageType == "Science":
            imagePaths = glob(pathToRaw + "/ScienceFrames/*.fits")
        else:
            imagePaths = glob(pathToRaw + f"/CalibrationImages/{imageType}/*.fits")

        Nimages = len(imagePaths)
        print(f"Found {Nimages} {imageType} frames to be inserted in MongoDB")
        if len(imagePaths) != 0:
            imageBaseNames = [os.path.splitext(os.path.basename(imagePath))[0] for imagePath in imagePaths]
            print(imageBaseNames)
            collection = db[imageType+"Images"]
            hashes = [createHash(imageBaseNames[n], f"Raw {imageType} Image") for n in range(Nimages)]
            currentTime = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
            imageInfo = [{"_id":          hashes[n], 
                          "path":         imagePaths[n],
                          "type":         f"Raw {imageType} Image",
                          "date_created": currentTime
                        } for n in range(Nimages)]
            collection.insert_many(imageInfo)





if __name__ == "__main__":

    initializeDataBaseWithRawImages(clearIfExist=True)
    
