from pipeline   import PipelineComponent
import yaml
import os
import tools
import numpy as np
import hashlib
import time







class MasterFlat(PipelineComponent):
    """
    Class that creates the master flat image. This image is
    obtained by taking the median of multiple flat images.
    """

    def __init__(self, debug=0, **flatAndBiasPaths):

        self.imageTypes = ["FlatImages", "BiasImages"]
        super().__init__(**flatAndBiasPaths)
        flatAndBiasPaths = self.inputPaths
        self.debug = debug

        self.masterBiasPath = flatAndBiasPaths["BiasImages"]
        self.rawFlatPath  = flatAndBiasPaths["FlatImages"]


















    def run(self, outputFileName=None):
        """
        We run through the alghoritm to create the master flat images.

        Input:
            outputFileName: If None, nothing is saved. Otherwise, a string with the name of the outputfile,
                            incl. the extension ".fits".

        Output:
            masterFlat:     master flat image [ADU]
        """

        # Get all the fits files corresponding to these hashes

        flats = tools.getImages(self.rawFlatPath)
        bias  = tools.getImage(self.masterBiasPath)
        stdBias = tools.getStdBias(self.masterBiasPath)
        meanBias = np.mean(bias)

        # Use the image in the fits files, and use mean_combining to obtain the the master image

        masterFlat = np.median(flats, axis=0) - bias

        # Add offset so that all the values in the MasterFlat are positive

        if np.min(masterFlat) < 0:
                  masterFlat = masterFlat -np.min(masterFlat)

        if outputFileName is not None:
            self.saveImageAndAddToDatabase(masterFlat, outputFileName, std_bias=stdBias, m_bias=meanBias)
            if (self.debug > 1):
                print("Master flat image saved to fits file")

        # That's it!

        return masterFlat







    def getHashOfOutputfile(self, imageHash=None):
        """
        The function returns the hash id for the output file.
        This hash is made from the hash files that are used as input.

        Ouput:
           hash. string containing the output hash
        """
        hash = hashlib.sha256(bytes("".join(self.rawFlatPath) + self.masterBiasPath, 'utf-8')).hexdigest()
        return hash


if __name__ == "__main__":

    t1 = time.time()

    params   = yaml.safe_load(open("params.yaml"))

    root     = (params["configuration"])["rootFolder"]
    f_params = params["rawFlatImage"]
    b_params = params["rawBiasImage"]

    # Mater Flat Image
    raw_flat_path = [ root+path for path in f_params["path"] ]
    master_bias_path = root + b_params["outputpath"]

    masterF = MasterFlat(FlatImages=raw_flat_path, BiasImages=master_bias_path)
    masterF.run(root+f_params["outputpath"])

    t2 = time.time()

    print(f"[{t2-t1}]")
