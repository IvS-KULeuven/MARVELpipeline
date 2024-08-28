import tools
import numpy as np
import hashlib
import time
import yaml

from pipeline import PipelineComponent



    

class MasterBias(PipelineComponent):
    """
    Class that creates the master bias image. Such an image is obtained by taking the median of multiple
    raw bias images.
    """

    def __init__(self, debug=0, **input):
        self.debug = debug
        self.imageTypes = ["BiasImages"]
        super().__init__(**input)
        input = self.inputPaths
        self.rawBiasPaths = input["BiasImages"]











    def run(self, outputFileName=None):
        """
        We run through the alghoritm to create the master bias images.

        Input:
            outputFileName: If None, nothing is saved. Otherwise, a string with the name of the outputfile,
                            incl. the extension ".fits".

        Output:
            masterBias:     master bias image [ADU]
        """
        
        # Get all the fits files corresponding to the input 
        biases = tools.getImages(self.rawBiasPaths)

        # Use the image in the fits files, and use mean_combining to obtain the the master image
        masterBias = np.median(biases, axis=0)
        stdBias = np.mean([np.std(b) for b in biases])

        if outputFileName is not None:
            self.saveImage(masterBias, outputFileName, std_bias=stdBias)
            if (self.debug >= 1):
                print("Master bias image saved to fits file")

        # That's it!

        return masterBias



    def getHashOfOutputfile(self, imageHash):
        """
        This function returns the hash id for the output file.
        This hash is made from the hash files that are used as input.

        Ouput:
            hash. string containing the output hash
        """
        hash = hashlib.sha256(bytes("".join(self.rawBiasPaths), 'utf-8')).hexdigest()
        return hash





if __name__ == "__main__":

    t1 = time.time()
    
    params = yaml.safe_load(open("params.yaml"))

    root  = (params["Configuration"])["rootFolder"]
    paths = (params["MasterBiasImage"])["inputpath"]
    
    paths = [ root+path for path in paths]
    output = root + (params["MasterBiasImage"])["outputpath"]

    # Master Bias Image

    masterB = MasterBias(BiasImages=paths)
    masterB.run(output)

    t2 = time.time()
    print(f"[{t2-t1}]")

