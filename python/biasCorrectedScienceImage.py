from pipeline   import PipelineComponent
import yaml
import tools
import hashlib
import time
import numpy as np




class BiasCorrectedScienceFrames(PipelineComponent):
    """
    Class that creates the bias corrected science image. These images are
    created by subtracting the median bias image from a science image.
    """

    def __init__(self, debug=0, **inputSciencePaths):
        self.imageTypes = ["BiasImages", "ScienceImages"]
        super().__init__(**inputSciencePaths)
        self.masterBiasPath  = inputSciencePaths["BiasImages"]
        self.rawSciencePaths = inputSciencePaths["ScienceImages"]



    def run(self, outputFileName=None):
        """
        We run through the algorithm to create the bias corrected science images.

        Input:
            outputFileName: If None, nothing is saved. Otherwise, a string with the name of the outputfile,
                            incl. the extension ".fits".

        Output:
            calibratedScience:     bias corrected science image [ADU]
        """

        # Get all the fits files corresponding to these hashes

        bias                = tools.getImage(self.masterBiasPath)
        self.masterBiasHash = tools.getHash(self.masterBiasPath)


        scienceImages = [tools.getImage(path) for path in self.rawSciencePaths]
        meanBias = np.mean(bias)

        biasSubtractedScienceImages = [ s - bias for s in scienceImages]

        # Add offset so that all the values in the MasterFlat are positive

        biasSubtractedScienceImages = [ s - np.min(s) for s in biasSubtractedScienceImages if np.min(s) < 0]

        # Convert float64 values in science images to int32 

        biasSubtractedScienceImages = [s.astype(np.int32) for s in biasSubtractedScienceImages]

        if outputFileName is not None:
            self.saveMultipleImages(biasSubtractedScienceImages, outputFileName, imageHashes=self.rawSciencePaths, m_bias=meanBias)


        # That's it!
        return biasSubtractedScienceImages






    def getHashOfOutputfile(self, rawScienceHash):
        """
        The function returns the hash id for the output file.
        This hash is made from the hash files that are used as input.

        Ouput:
           hash. string containing the output hash
        """

        combinedHashes =  self.masterBiasHash + rawScienceHash
        hash = hashlib.sha256(bytes(combinedHashes, 'utf-8')).hexdigest()
        return hash



if __name__ == "__main__":

    t1 = time.time()

    params   = yaml.safe_load(open("params.yaml"))

    root     = (params["configuration"])["rootFolder"]
    s_params = params["rawScienceImage"]
    b_params = params["rawBiasImage"]
    
    # Bias corrected Science Image

    rawSciencePath = [root+path for path in s_params["path"]]
    master_bias_path = root + b_params["outputpath"]

    calibration = BiasCorrectedScienceFrames(ScienceImages=rawSciencePath, BiasImages=master_bias_path)
    calibration.run([root+path for path in s_params["outputpath"]])

    t2 = time.time()

    print(f"[{t2-t1}]")
