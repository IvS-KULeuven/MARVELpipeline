# This is the parent class of the all pipeline components + components that generate master images



import os
from pathlib import Path
#from pymongo import MongoClient
from astropy.io import fits
import tools
import numpy as np
import hashlib
from datetime import datetime














class PipelineComponent():
    """
    MARVEL pipeline building components

    The PipelineComponent represents one isolated component of the pipeline.
    A component is created using a list of paths that represents other
    component-outputs.

    Once we have created a component the output is output is stored in a file
    and saved.

    Example:
        masterComponent = MasterBias(BiasImages=path_list)
        masterComponent.run()
    """


    def __init__(self, **inputPaths):
        """
        Initializes the component of the pipeline. This functions is called after a component
        is generated. This checks that the input hashes are sensible.

        Input:
            inputHashes: imageType/path of the FITS file that is given as input.
        """
        
        if not (self.checkSanityOfInput(**inputPaths)):
            keysExists = np.array([key in self.imageTypes for key in inputPaths.keys()])

            # First we check if the keys that are given are part of ImageTypes

            if not np.all(keysExists):
                unrecognizedKeys = np.array(list(input.keys()))[~keysExists]
                unrecognizedKeys = list(unrecognizedKeys)
                keys = str(unrecognizedKeys)[1:-1]
                raise Exception(f"Keyvalues: {keys} not found in database")

            else:

                # If the keys are correct, something must be wrong with the values.

                wrongPaths = []

                for key in inputPaths.keys():
                    if isinstance(inputPaths[key], str):
                        if not os.path.isfile(inputPaths[key]):
                            wrongPaths.append(inputPaths[key])
                    elif isinstance(inputPaths[key], list):
                        for path in inputPaths[key]:
                            if not os.path.isfile(path):
                                wrongPaths.append(path)
                    else:
                        raise Exception(f"Value of {key} is in the wrong format.")
                if len(wrongPaths) == 1:
                    raise Exception(f"The path: {wrongPaths[0]} is not found")
                else:
                    paths = " ".join(wrongPaths)
                    raise Exception(f"The paths: {paths} are not found")

        self.inputPaths = inputPaths








    def checkSanityOfInput(self, **inputPaths):
        """
        Check if the image hashes indeed corresponds to an image of the correct format (given by
        imageType).

        Input:
             **inputPaths: imageType/paths of the images given as input.

        Output:
            isSane: boolean. True if both input hashes correspond to images with the correct type, False if at
                             least one of these hashes does not.
        """

        for imageType in inputPaths:
            if imageType not in self.imageTypes:
                return False
            else:

                # If there are multiple hashes corresponding to one imagetype

                if isinstance(inputPaths[imageType], list):
                    for path in inputPaths[imageType]:
                        isPath = os.path.isfile(path)
                        if not isPath:
                            return False

                # If there is one hash corresponding to one imageType

                else:
                    isPath = os.path.isfile(inputPaths[imageType])
                    if not isPath:
                        return False

        # If we get here, the hashes were found at the correct location

        return True








    def saveMultipleImages(self, images, fileNames, imageHashes, **keywords):
        """
        Save multiple images

        input:
            image: list of np.arrays of the images we want to save to fits files
            fileNames: paths of the fits files
            keyswords: dictionary with option parameters we want to add to the fits header
        """
        for file, (path, hash) in zip(images, zip(fileNames, imageHashes)):
            self.saveImage(file, path, imageHash=hash, **keywords)
        



    def saveImage(self, image, fileName, imageHash=None, **keywords):
        """
        Save the image

        input:
            image: np.array of the image that we want to save in fits files
            filname: path of the fits file
            keywords: dictionary with optional parameters we want to add to the fits header
        """
        hash = self.getHashOfOutputfile(imageHash)

        path = fileName
        num_row, num_col = image.shape

        # Save Master Image as FITS file
        hdr = fits.Header()
        hdr["hash"] = hash
        hdr["path"] = path
        hdr["rows"] = num_row
        hdr["cols"] = num_col

        for key, value in keywords.items():
            hdr[key] = value
            
        output = Path(path).parent.absolute()
        if not os.path.exists(output):
            os.makedirs(output)


        hdu = fits.PrimaryHDU(image, header=hdr)
        hdu.writeto(path, overwrite=True)





