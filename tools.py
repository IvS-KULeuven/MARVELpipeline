# We use this module for opperations that will be commonly used in different pipelineModules




from astropy.io import fits
import numpy as np
from pymongo import MongoClient

client = MongoClient()
db = client["databaseMARVEL"]




def saveFigure(figure):
    ...


def getImages(paths):
    return np.array([getImage(path) for path in paths])



def getImage(path):
    hdul = fits.open(path)
    return hdul[0].data




def addToDataBase(dict, overWrite=False):
    # We should have proper error catching for:
    # 1. We want to add somthing that already exist
    # 2. image is not in the dictornary
    # 3. Check that dict is in the right format 

    images     = {"Master Dark Image": "DarkImages", "Master Bias Image": "BiasImages", "Master Flat Image": "FlatImages", "Calibrated Science Image" : "ScienceImages"}
    typeImage  = dict["type"]
    collection = db[images[typeImage]]

    isInCollection = np.all([x == dict for x in collection.find({"_id" : dict["_id"]})])

    if not isInCollection:
        collection.inster_one(dict)

    elif overWrite:
        collection.delete_many({"_id": dict["_id"]})
        collection.insert_one(dict)
    else:
        print("Document is already in database")
    
    
        




if __name__ == "__main__":
    path = "/lhome/driess/MARVEL/MARVELpipe/Data/RawData/CalibrationImages/Bias/bias_0001.fits"
    print(getImage(path))
