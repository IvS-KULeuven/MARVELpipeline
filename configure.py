import errno
import yaml
import hashlib
import os
import sys
from pathlib import Path



def isCorrectPath(path):
    return os.path.exists(path)







def readInInputFile(filename):
    """
    Read in the user defined input file. This function checks that the
    input file is sensible and returns a dict with the values from the file.

    Parameters:
        filename: Name of the input file we should parse.

    Output:
        dict with the values of the input file.
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
    expectedKeys = {'rawScienceImage', 'rawFlatImage', 'rawBiasImage'}
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

            stripped_list = [ path for path in yaml_input[keyValue] if isCorrectPath(path) ]

            if len(stripped_list) == 0:
                raise Exception(f"Files for {keyValue} do not exisits.")

            elif len(stripped_list) != len(yaml_input[keyValue]):
                print(f"WARNING: Not every file specified in {keyValue} exists.\nWe continue the simulation with the files\n\n{stripped_list}")
                yaml_input[keyValue] = stripped_list

        elif isinstance(yaml_input[keyValue], str):

            if not isCorrectPath(yaml_input[keyValue]):
                raise Exception(f"File {yaml_input[keyValue]} for {keyValue} does not exists.")

            else:
                yaml_input[keyValue] = [ yaml_input[keyValue] ]

    return yaml_input
                









            
        
def create_dvc_file(yaml_input, dvc_param="params.yaml"):
    """
    Create the dvc input file that can be used to run dvc.

    Parameters:
        yaml_input: Dict with the paths to the raw calibration images and the
                    science images.
        dvc_param: Name of the path to the dvc input file. If none is give it
                   uses the default params.yaml value.
    """

    home = os.getcwd() + "/"

    masterBias        = createCalibrationOutputPath(yaml_input["rawBiasImage"])
    masterFlat        = createCalibrationOutputPath(yaml_input["rawFlatImage"])
    calibratedScience = createScienceOutputPath(yaml_input["rawScienceImage"])

    smoothMaster = createMaskOutputPathPrefix(masterFlat) + "_smoothed_master_flat.fits"
    maskPath     = createMaskOutputPathPrefix(masterFlat) + "_2d_mask.fits"
    scienceMask  = createScienceMaskOutputPath(calibratedScience)
    oneDOrder    = createOptimalOrderExtractionOutputPath(calibratedScience)
    

    param_yaml = {"rawBiasImage" :
                    { "path" : yaml_input["rawBiasImage"],
                      "outputpath" : masterBias
                    },
                  "rawFlatImage" :
                    { "path" : yaml_input["rawFlatImage"],
                      "outputpath" : masterFlat
                    },
                  "rawScienceImage" :
                    { "path" : yaml_input["rawScienceImage"],
                      "outputpath" : calibratedScience
                    },
                  "orderImage" :
                    { "smoothMasterFlat" : smoothMaster,
                      "maskOutputpath" : maskPath,
                      "scienceOutputpath" : scienceMask
                    },
                  "optimalOrderExtraction" :
                    { "outputpath" : oneDOrder
                    },
                  "configuration" :
                    { "rootFolder" : home
                    }
                 }


    file=open(dvc_param,"w")
    yaml.dump(param_yaml,file)
    file.close()

    











def createCalibrationOutputPath(inputPaths):
    """
    Function that generates one output path given a list of calibration (either bias or flat) files.

    Parameters:
        inputPaths: List of paths that correspond to raw bias/flat images

    Output:
        output path for a master bias/flat image. This should be unique for different
        inputs.
    """

    # Combine the sorted list in one string to generate a unique hash
    
    inputPaths.sort()

    # We use the date of the earliest file

    fileName = inputPaths[0].split("/")[-1]
    dateTime = fileName.split("_")[0]
    date = dateTime.split("T")[0]
    
    # Specify the type of file

    imageType = fileName.split("_")[1]
    if imageType == "BBBBB":
        outputType = "MasterBias"
    elif imageType == "FFFFF":
        outputType = "MasterFlat"
    elif imageType == "DDDDD":
        outputType = "MasterDark"

    # Get the root directory

    outputPath = "Data/ProcessedData/" + outputType + "/" + date + "_" + outputType + ".fits"

    return outputPath







def createScienceOutputPath(inputPaths):
    """
    Function that generates one output path given a list of science images.

    Parameters:
        inputPaths: List of paths that correspond to raw science images.

    Output:
        list of paths to be used to save bias subtracted science images
    """

    root = "Data/ProcessedData/BiasSubtractedScience/"
    output_paths = []

    for path in inputPaths:
        fileStem = Path(path).stem
        output_paths.append(root + fileStem + "_bias_subtracted.fits")

    return output_paths








def createMaskOutputPathPrefix(masterFlatPath):
    """
    Create the first part of the path for the mask and the smooth master flat images

    Parameters:
        masterFlatPath: path of the master flat file.

    """

    # Root of the file

    root = "Data/ProcessedData/ExtractedOrders/Mask/"

    # Get date of the files

    fileStem = Path(masterFlatPath).stem
    dateTime = fileStem.split("_")[0]
    output_path = root + dateTime

    return output_path

    








def createScienceMaskOutputPath(scienceImagePaths):
    """
    Create a list with the paths to the science mask.

    Parameters:
        scienceImagePaths: List with (absolute or relative) paths to bias subtracted science images.
    """
    output_paths = []
    root = "Data/ProcessedData/ExtractedOrders/Science/"
    for path in scienceImagePaths:
        fileStem = Path(path).stem.replace("_bias_subtracted", "")
        output_paths.append(root + fileStem + "_2d_orders.fits")

    return output_paths
    





def createOptimalOrderExtractionOutputPath(scienceImagePaths):
    """
    Create the path for optimal order extraction output.

    Parameters:
        scienceImagePaths: List of paths to the bias subtracted science images. 
    """

    output_paths = [] 
    root = "Data/ProcessedData/OptimalExtraction/"
    for path in scienceImagePaths:
        fileStem = Path(path).stem.replace("_bias_subtracted", "")
        output_paths.append(root + fileStem + "_1d_orders.fits")

    return output_paths
    
    
    




def createEtalonPeakFittingOutputPath(optimalExtracted1DspectrumPaths):
    """
    Create the path of the FITS files containing the etalon peak fitting results.

    Input:
        - optimalExtracted1DspectrumPaths: List of paths to the FITS files containing the optimally extracted 1D spectra 
    """

    output_paths = []
    root = "Data/ProcessedData/WaveCalibration/"

    for path in optimalExtracted1DspectrumPaths:
        fileStem = Path(path).stem
        output_paths.append(root + fileStem + "etalon_peak_fitparameters.fits")

    return output_paths









if __name__ == "__main__":

    if len(sys.argv) == 1:
        raise TypeError("Expected 1 arguments, got 0.\nPlease provide input file argument as:\n\tpython marvel.py inputfile.yaml\nor\n\tpython marvel.py inputfile.yaml outputfile.yaml")
    else:
        input = readInInputFile(sys.argv[1])

    if len(sys.argv) == 2:
        create_dvc_file(input)
    else:
        create_dvc_file(input, sys.argv[2])
        




