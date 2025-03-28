import errno
import yaml
import os
import sys
from pathlib import Path





def verifyInputfiles(filename, rootFolderRawData):
    """
    Read in the user defined input file. This function checks that the
    input file is sensible and returns a dict with the values from the file.

    Parameters:
        filename: Name of the input file we should parse.

    Output:
        dict with the values of the input file.

    Example:
        readInputfile("inputfile.yaml")
    """
    # check the file exicts

    file_exist = os.path.exists(filename)

    if not file_exist:
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), filename)

    # If file exists open the file

    with open(filename, 'r') as stream:
        try:
            yaml_input = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)


    # Check that the file has the right keys
    
    keys = set(yaml_input.keys())
    expectedKeys = {'rawScienceImages', 'rawFlatImages', 'rawBiasImages', 'rawDarkImages', 'rawThArImages'}
    if not keys == expectedKeys:
        raise Exception(f"\n\nExpected keys: {expectedKeys}\n Keys found in file are: {keys}")
        
    # Parse the values and check that the files exist

    for keyValue in expectedKeys:
        
        # The value might be a list or a string.
        
        # If they are a string we replace the string by a list containing the single string
        # if the string is a path that exisits. Else we raise an exception.

        # If they are a list, we check that every item in the list is a path that exists. If not
        # we replace the list with a stripped list containg only those paths that exisits. If no
        # existing path remains, we raise an exception.

        if isinstance(yaml_input[keyValue], list):

            stripped_list = [rootFolderRawData + "/" + path for path in yaml_input[keyValue]    \
                             if os.path.exists(rootFolderRawData + "/" + path) ]    

            if len(stripped_list) == 0:
                raise Exception(f"File(s) for {keyValue} do not exists.")

            elif len(stripped_list) != len(yaml_input[keyValue]):
                print(f"WARNING: Not every file specified in {keyValue} exists.\nWe continue with the files\n\n{stripped_list}")
                yaml_input[keyValue] = stripped_list

        elif isinstance(yaml_input[keyValue], str):

            if not os.path.exists(yaml_input[keyValue]):
                raise Exception(f"File {yaml_input[keyValue]} for {keyValue} does not exists.")

            else:
                yaml_input[keyValue] = [ yaml_input[keyValue] ]

    return yaml_input
                









def createCalibrationOutputPath(inputPaths):
    """
    Function that generates one output path given a list of calibration file. This can
    be either Bias, Flatfield, Dark, or ThAr image files

    Parameters:
        inputPaths: List of paths that correspond to raw bias/flat images

    Output:
        string: output path for a master calibration image. 
                Each type of calibration file (Bias, Flat, Dark, ThAr) has a unique output path. 
    """

    # Combine the sorted list in one string to generate a unique hash
    
    inputPaths.sort()

    # We use the date of the earliest file

    fileStem = Path(inputPaths[0]).stem
    parts = fileStem.split("_")
    dateTime = parts[0]
    imageType = parts[1]
    fiberType = parts[2]
    exposureTime = parts[3]
    date = dateTime.split("T")[0]
    
    # Specify the type of file

    if imageType == "BBBBB":
        outputType = "MasterBias"
    elif imageType == "FFFFF":
        outputType = "MasterFlat"
    elif imageType == "DDDDD":
        outputType = "MasterDark"
    elif imageType.startswith("T"):
        outputType = "MasterThAr"

    # Create the output path
    
    outputPath = "/ProcessedData/" + outputType + "/" + date + "_" + outputType 
    if outputType == "MasterThAr":
        outputPath += "_" + fiberType 
    elif outputType == "MasterDark":
        outputPath += "_" + exposureTime 

    outputPath += ".fits"

    return outputPath







def createScienceOutputPath(inputPaths):
    """
    Function that generates one output path given a list of science images.

    Parameters:
        inputPaths: List of paths that correspond to raw science images.

    Output:
        list of paths to be used to save bias subtracted science images
    """

    root = "/ProcessedData/BiasAndDarkSubtractedScience/"
    output_paths = []

    for path in inputPaths:
        fileStem = Path(path).stem
        output_paths.append(root + fileStem + "_biasdark_subtracted.fits")

    return output_paths








def createMaskOutputPathPrefix(masterFlatPath):
    """
    Create the first part of the path for the mask and the smooth master flat images

    Parameters:
        masterFlatPath: path of the master flat file.

    """

    # Root of the file

    root = "/ProcessedData/OrderMask/"

    # Get date of the files

    fileStem = Path(masterFlatPath).stem
    dateTime = fileStem.split("_")[0]
    output_path = root + dateTime

    return output_path

    








def createScience2DordersOutputPath(scienceImagePaths):
    """
    Create a list with the paths to the science mask.

    Parameters:
        scienceImagePaths: List with (absolute or relative) paths to bias subtracted science images.
    """
    output_paths = []
    root = "/ProcessedData/ExtractedOrders/Science/"
    for path in scienceImagePaths:
        fileStem = Path(path).stem.replace("_biasdark_subtracted", "_2d_orders")
        output_paths.append(root + fileStem + ".fits")

    return output_paths
    








def createThAr2DordersOutputPath(masterThArPath):
    """
    Create a list with the paths to the ThAr 2D orders

    Parameters:
        masterThArPath: string: path to the masterThArPath

    Output: 
        string: path to the FITS file with the 2D orders of the master ThAr image
    """
    root = "/ProcessedData/ExtractedOrders/ThAr/"
    outputPath = root + Path(masterThArPath).stem + "_2d_orders.fits"

    return outputPath
    







def createScienceOptimalOrderExtractionOutputPath(bias_subtracted_paths):
    """
    Create the path for optimal order extraction output.

    Parameters:
        bias_subtracted_paths: List of paths to the fits files with the two dimensional orders

    Output:
        list of strings: paths to fits files with the extracted 1D orders
    """

    output_paths = [] 
    root = "/ProcessedData/OptimalExtraction/"
    for path in bias_subtracted_paths:
        fileStem = Path(path).stem.replace("_biasdark_subtracted", "_1d_orders")
        output_paths.append(root + fileStem + ".fits")

    return output_paths
    
    
    



def createOneDimOrderExtractedMasterFlatPath(master_flat_path):
    """
    Create the path for the 1D order-extracted master flat field.

    Parameters:
        master_flat_path: (relative) path to the 2D master field image

    Output:
        path the FITS file containing the 1D order-extracted master flat field
    """

    root = "/ProcessedData/MasterFlat/"
    output_path = Path(master_flat_path).stem.replace("MasterFlat", "1d_order_MasterFlat")
    output_path = root + output_path + ".fits"
    return output_path






def createThArOptimalOrderExtractionOutputPath(thar_master_path):
    """
    Create the path for optimal order extraction output.

    Parameters:
        twoDimOrdersPaths: List of paths to the fits files with the two dimensional orders

    Output:
        string: path to fits files with the extracted 1D orders of the ThAr master file
    """

    root = "/ProcessedData/OptimalExtraction/"
    output_path = root + Path(thar_master_path).stem + "_1d_orders.fits"
    return output_path
    





def createEtalonPeakFittingOutputPath(optimalExtracted1DspectrumPaths):
    """
    Create the path of the FITS files containing the etalon peak fitting results.

    Input:
        - optimalExtracted1DspectrumPaths: List of paths to the FITS files containing the optimally extracted 1D spectra 
    """

    output_paths = []
    root = "/ProcessedData/WaveCalibration/"

    for path in optimalExtracted1DspectrumPaths:
        fileStem = Path(path).stem.replace("_1d_orders", "_etalon_peak_fitparameters")
        output_paths.append(root + fileStem + ".fits")

    return output_paths









def create_yaml_configfile(yaml_input, rootFolderRawData, rootFolderProcessedData, outputfile="config.yaml"):
    """
    Create the dvc input file that can be used to run dvc.

    Parameters:
        yaml_input: Dict with the paths to the raw calibration images and the
                    science images.
        rootFolderRawData: Full path to the base folder of the raw images. 
        rootFolderProcessedData: Full path to the base folder of the outputfiles. 
        outputfile: Name of the path to the pipeline config file. If none is give it
                   uses the default "config.yaml" name.

    Example:
        create_dvc_inputfile("inputfile.yaml", outputfile="config.yaml")
        -> This will create a "config.yaml" file.
    """

    masterBiasPath = createCalibrationOutputPath(yaml_input["rawBiasImages"])
    masterDarkPath = createCalibrationOutputPath(yaml_input["rawDarkImages"])
    masterFlatPath = createCalibrationOutputPath(yaml_input["rawFlatImages"])
    masterThArPath = createCalibrationOutputPath(yaml_input["rawThArImages"])
    biasAndDarkSubtractedSciencePaths = createScienceOutputPath(yaml_input["rawScienceImages"])

    smoothMasterPath            = createMaskOutputPathPrefix(masterFlatPath) + "_smoothed_master_flat.fits"
    maskBoundariesPath          = createMaskOutputPathPrefix(masterFlatPath) + "_2d_mask_boundaries.fits"
    maskVisualisationImagePath  = createMaskOutputPathPrefix(masterFlatPath) + "_2d_mask_image.fits"
    oneDimScienceOrdersPaths    = createScienceOptimalOrderExtractionOutputPath(biasAndDarkSubtractedSciencePaths)
    oneDimOrderExtractedMasterFlatPath = createOneDimOrderExtractedMasterFlatPath(masterFlatPath)
    oneDimThArOrdersPath        = createThArOptimalOrderExtractionOutputPath(masterThArPath)
    etalonPeakFitParametersPath = createEtalonPeakFittingOutputPath(oneDimScienceOrdersPaths)

    config_yaml = {"Configuration":
                  { 
                    "rootFolderRawData": rootFolderRawData,
                    "rootFolderProcessedData": rootFolderProcessedData,
                  },

                  "MasterBiasImage":
                  { 
                    "inputPath": yaml_input["rawBiasImages"],
                    "outputPath": masterBiasPath
                  },

                  "MasterDarkImage":
                  { 
                    "inputPath": yaml_input["rawDarkImages"],
                    "outputPath": masterDarkPath
                  },

                  "MasterFlatImage":
                  { 
                    "inputPath": yaml_input["rawFlatImages"],
                    "outputPath": masterFlatPath
                  },

                  "MasterThArImage":
                  { 
                    "inputPath": yaml_input["rawThArImages"],
                    "outputPath" : masterThArPath
                  }, 

                  "BiasAndDarkSubtractedScienceImage":
                  {
                    "inputPath" : yaml_input["rawScienceImages"],
                    "outputPath" : biasAndDarkSubtractedSciencePaths
                  },

                  "TwoDimensionalOrderMaskTracing":
                  { 
                    "inputPath": masterFlatPath,
                    "outputPath": maskBoundariesPath,
                    "outputPathSmoothedMasterFlat": smoothMasterPath
                  },

                  "TwoDimensionalOrderMaskVisualisation":
                  {
                      "inputPath": maskBoundariesPath,
                      "outputPath": maskVisualisationImagePath,
                  },

                  "OptimalOrderExtraction":
                  { 
                    "inputPath": biasAndDarkSubtractedSciencePaths + [masterThArPath],
                    "outputPath": oneDimScienceOrdersPaths + [oneDimThArOrdersPath]
                  },

                  "OneDimensionalOrderExtractedMasterFlat":
                  {
                    "inputPathMasterFlat": masterFlatPath,
                    "inputPathMask": maskBoundariesPath,
                    "outputPath": oneDimOrderExtractedMasterFlatPath,
                  },

                  "EtalonPeakFitting":
                  { 
                    "inputPath": oneDimScienceOrdersPaths,
                    "outputPath": etalonPeakFitParametersPath
                  }
                 }


    file=open(outputfile, "w")
    yaml.dump(config_yaml,file, sort_keys=False)
    file.close()







"""
Example usage:

    $ python configure.py inputfile.yaml config.yaml /Users/joris/MARVELpipeline/Data/  /Users/joris/MARVELpipeline/Data/

The raw data will then be looked for in /Users/joris/MARVELpipeline/Data/RawData  and the output files of the 
pipeline will be written in /Users/joris/MARVELpipeline/Data/ProcessedData. The paths in inputfile.yaml and
config.yaml are relative to these root paths. 
"""
if __name__ == "__main__":

    if len(sys.argv) != 5:
        print("Usage: python configure.py <inputfile.yaml> <outputfile.yaml> <root folder raw images> <root folder output files>")
        exit(0)
    else:
        rootFolderRawData = sys.argv[3]
        rootFolderProcessedData = sys.argv[4]
        yamlInput = verifyInputfiles(sys.argv[1], rootFolderRawData)
        yamlOutputFile = sys.argv[2]

    create_yaml_configfile(yamlInput, rootFolderRawData, rootFolderProcessedData, outputfile=yamlOutputFile)




