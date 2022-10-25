
import hashlib
import os
import pandas as pd
import numpy  as np
from pymongo  import MongoClient
from glob     import glob
from datetime import datetime
from io       import StringIO






class DatabaseFromLocalFiles():
    """
    This class is used for locally running the pipeline without needing to install
    mongodb. This is only advisable if you are dealing with a small amount of fits
    files. The database object in this case is a pandas dataframe.
    """

    def __init__(self, initializationFilePath=None):
        """
        Initializes the database object.

        Input:
            initializationFilePath: If None, the database gets created from the files that are in
                                that raw data directory. Otherwise, a string that points to a
                                file from wich the database gets created.
        """

        if initializationFilePath is None:

            # Create database from the files in the raw data directory
            dataDirectoryPath = os.getcwd() + "/Data/RawData/"
            self.localDataBase = self.makeDataBaseLocally(dataDirectoryPath)


        else:

            # Check that file exists
            if not os.path.isfile(initializationFilePath):
                print(f"Input hash: {initializationFilePath} does not correspond to vallid file")
            else:
                self.localDataBase = self.loadFromFile(initializationFilePath)






    def makeDataBaseLocally(self, path):
        """
        The database gets created from the files that are in the directory that path
        points to.

        Input:
            path: string that points to the paths where the raw files are located.

        Ouptut:
            database: pandas dataframe with the information of the local raw fils.
        """
        database = {}


        # Search for the images at the path
        print(f"Searching for raw images in {path}")
        for imageType in ["Bias", "Dark", "Flat", "Etalon", "Science"]:
            if imageType == "Science":
                imagePaths = glob(path + "ScienceFrames/*.fits")

            else:
                imagePaths = glob(path + f"CalibrationImages/{imageType}/*.fits")

            Nimages = len(imagePaths)
            print(f"Found {Nimages} {imageType} frames to be inserted")
            if Nimages != 0:
                imageBaseNames = [os.path.splitext(os.path.basename(imagePath))[0] for imagePath in imagePaths]

                hashes = [createHash(imageBaseNames[n], f"Raw {imageType} Image") for n in range(Nimages)]
                currentTime = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
                dat = { "_id": hashes, "path": imagePaths, "type": [f"Raw {imageType} Image"]*Nimages, "date_created": [currentTime]*Nimages}
                df = pd.DataFrame.from_dict(data=dat)
                database[imageType+"Images"] = Collection(df)


        return database







    def __getitem__(self, item):
        if not item in self.localDataBase.keys():
            dat = { "_id": [], "path": [], "type": [], "date_created": []}
            self.localDataBase[item] = Collection(pd.DataFrame.from_dict(data=dat))
        return self.localDataBase[item]



    def saveToFile(self, outputFilePath):
        """
        Save the information in the database object to text file.

        Input
            outputFilePath. string with the path of the file
        """

        keys = list(self.localDataBase.keys())
        values = list(self.localDataBase.values())
        values = [value.df.to_csv() for value in values]

        fileAsDictionary = dict(zip(keys, values))

        file = open(outputFilePath, "w")
        file.write(str(fileAsDictionary))
        file.close()





    def loadFromFile(self, filePath):
        """
        Initialize the database object with the data from a file.

        Input
            filepath. string with the path of the file
        """

        # Read in the dictionary string
        with open(filePath, "r") as f:
            lines = f.readlines()

        # Sanity check
        if not len(lines) == 1:
            raise Exception("File is not readable")
            exit(1)

        # Initialize database dictionary object
        database = {}

        # Seperate the keys and values from the dictionary string
        listOfLines = lines[0].split("'")[1::2]
        keys = listOfLines[::2]
        values = listOfLines[1::2]

        # Add the keys/value in a dictionary
        for key, value in zip(keys, values):

            # Generate the dataframe to initialize the collection object
            headers = (value.split("\\n")[:-1:])[0]
            headers = np.array(headers.split(",")[1::])

            data    = (value.split("\\n")[:-1:])[1::]
            data    = np.array([val.split(",")[1::] for val in data])

            dat = dict(zip(headers, data.transpose()))
            df  = pd.DataFrame.from_dict(data=dat)

            # Add to dictionary object
            database[key] = Collection(df)

        # Return the database property with the values from the file

        print("File loaded")
        return database


















class Collection:
    """
    Simulate objects simlar to collection in mongoDB.
    """

    def __init__(self, dataFrame):
        self.df = dataFrame




    def find_one(self, filter):
        """
        Finds the value in the database collection for the key

        Input:
            filter: with {key : values} to find in database.
        """
        keys, values = zip(*filter.items())
        key = keys[0]
        value = values[0]

        if not key in self.df.columns:
            print(f"{key} is not in the database")
            return None

        # Find the data in the database

        matchedData = self.df.loc[self.df[key] == value]

        if matchedData.empty:
            return None

        else:
            return self.df.loc[self.df[key] == value].to_dict('records')[0]



    def insert_one(self, document):
        """
        Adds the dictionary to the database collection object

        Input: document. dictionary items with keywords the columns values of the dataframe.

        Note:  we assume the value that we want to add is not already present in the database
        """

        # check that keywords equal the column values of the dataframe
        if list(document.keys()).sort() == list(self.df.columns).sort():

            # add this to the database
            df1 = pd.DataFrame(data=document, index=[len(self.df)])
            self.df = pd.concat([self.df, df1])

        else:
            print("WARNING: dictionary format is not in the right format.")
            return




    def delete_many(self, filter):
        """
        Deletes the row that in which the dictionary key/value can be matched.

        Input: filter. dictionary item that contains criteria of the line to delete
        """

        keys, values = zip(*filter.items())
        key = keys[0]
        value = values[0]

        self.df = self.df[self.df[key] != value]



    def __str__(self):
        return "{}".format(self.df)












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
        if Nimages != 0:
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

    #    initializeDataBaseWithRawImages(clearIfExist=True)
    db = DatabaseFromLocalFiles()
    db.saveToFile("Test.txt")
    db.loadFromFile("Test.txt")
    db.saveToFile("Test2.txt")
