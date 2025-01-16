
import time
from datetime import datetime
import yaml
import argparse as ap
import textwrap
import subprocess
from masterBias import MasterBias 
from masterDark import MasterDark
from masterFlat import MasterFlat
from masterThAr import MasterThAr
from biasAndDarkCorrection import BiasAndDarkCorrectedScienceFrames
from mask_fits_image import MaskImageCreator



def masterBias(config):

    rootFolderRawData  = config["Configuration"]["rootFolderRawData"]
    rootFolderProcessedData = config["Configuration"]["rootFolderProcessedData"]
    relativeBiasImagePaths = config["MasterBiasImage"]["inputPath"]
    absoluteBiasImagePaths = [ rootFolderRawData + path for path in relativeBiasImagePaths]
    absoluteMasterBiasImagePath = rootFolderProcessedData + config["MasterBiasImage"]["outputPath"]
    masterB = MasterBias(BiasImages=absoluteBiasImagePaths)
    masterB.run(absoluteMasterBiasImagePath)




def masterDark(config):

    rootFolderRawData  = config["Configuration"]["rootFolderRawData"]
    rootFolderProcessedData = config["Configuration"]["rootFolderProcessedData"]
    master_dark_params = config["MasterDarkImage"]
    master_bias_params = config["MasterBiasImage"]
    raw_dark_paths = [ rootFolderRawData + path for path in master_dark_params["inputPath"] ]
    master_bias_path = rootFolderProcessedData + master_bias_params["outputPath"]
    masterDark = MasterDark()
    masterDark.run(raw_dark_paths, master_bias_path, rootFolderProcessedData + master_dark_params["outputPath"])




def masterFlat(config):

    rootFolderRawData  = config["Configuration"]["rootFolderRawData"]
    rootFolderProcessedData = config["Configuration"]["rootFolderProcessedData"]
    flat_params = config["MasterFlatImage"]
    bias_params = config["MasterBiasImage"]
    dark_params = config["MasterDarkImage"]
    raw_flat_paths = [ rootFolderRawData + path for path in flat_params["inputPath"] ]
    master_bias_path = rootFolderProcessedData + bias_params["outputPath"]
    master_dark_path = rootFolderProcessedData + dark_params["outputPath"]
    master_flat_path = rootFolderProcessedData + flat_params["outputPath"] 
    masterFlat = MasterFlat()
    masterFlat.run(raw_flat_paths, master_bias_path, master_dark_path, outputFilePath=master_flat_path)




def masterThAr(config):

    rootFolderRawData  = config["Configuration"]["rootFolderRawData"]
    rootFolderProcessedData = config["Configuration"]["rootFolderProcessedData"]
    masterThAr_params = config["MasterThArImage"]
    masterbias_params = config["MasterBiasImage"]
    masterdark_params = config["MasterDarkImage"]
    raw_ThAr_paths = [ rootFolderRawData + path for path in masterThAr_params["inputPath"] ]
    master_bias_path = rootFolderProcessedData + masterbias_params["outputPath"]
    master_dark_path = rootFolderProcessedData + masterdark_params["outputPath"]
    ThAr_output_path = rootFolderProcessedData + masterThAr_params["outputPath"]
    masterThAr = MasterThAr()
    masterThAr.run(raw_ThAr_paths, master_bias_path, master_dark_path, ThAr_output_path)




def biasAndDarkCorrection(config):

    rootFolderRawData  = config["Configuration"]["rootFolderRawData"]
    rootFolderProcessedData = config["Configuration"]["rootFolderProcessedData"]
    science_params = config["BiasAndDarkSubtractedScienceImage"]
    bias_params = config["MasterBiasImage"]
    dark_params = config["MasterDarkImage"]
    raw_science_paths = [rootFolderRawData + path for path in science_params["inputPath"]]
    master_bias_path = rootFolderProcessedData + bias_params["outputPath"]
    master_dark_path = rootFolderProcessedData + dark_params["outputPath"]
    output_science_paths = [rootFolderProcessedData + path for path in science_params["outputPath"]]
    calibration = BiasAndDarkCorrectedScienceFrames()
    calibration.run(raw_science_paths, master_bias_path, master_dark_path, outputPaths=output_science_paths)





def maskDetermination(config):
    subprocess.run(["./pipeline/target/release/2dmaskdetermination", "--configpath", config['origin']])

    # The following serves to past the masks into a fits image that can be visualized for quality control
    rootFolderProcessedData = config["Configuration"]["rootFolderProcessedData"]
    mask_visualization_params = config["TwoDimensionalOrderMaskVisualisation"]
    inputPath = rootFolderProcessedData + mask_visualization_params["inputPath"]
    outputPath = rootFolderProcessedData + mask_visualization_params["outputPath"]
    maskImageCreator = MaskImageCreator()
    maskImageCreator.run(inputPath, outputPath)




def oneDimOrderExtraction(config):
    capture = subprocess.run(["./pipeline/target/release/1dflatrelativeorderextraction", "--configpath", config['origin']], \
                             capture_output=True, text=True)
    # Uncomment the following lines to debug
    # print(capture.stdout)
    # print(capture.stderr)



def etalonPeakFitting(config):
    capture = subprocess.run(["./pipeline/target/release/fitetalonpeaks", "--configpath", config['origin']],  \
                             capture_output=True, text=True)
    # Uncomment the following lines to debug
    print(capture.stdout)
    print(capture.stderr)








if __name__ == "__main__":

    # These are the different steps of the MARVEL data reduction pipeline 

    pipeline_steps = {  1: [masterBias,            'Compute master bias image'],
                        2: [masterDark,            'Compute master dark image'],
                        3: [masterFlat,            'Compute master flat image'],
                        4: [masterThAr,            'Compute master ThAr image'],
                        5: [biasAndDarkCorrection, 'Correct CCD images for bias and dark'],
                        6: [maskDetermination,     'Determine 2D mask of the orders'],
                        7: [oneDimOrderExtraction, 'Extract 1D orders'],
                        8: [etalonPeakFitting,     'Fit etalon peaks']}


    # Create some help text when the --help or -h option was specified

    description = "MARVEL Data Reduction Pipeline Steps: \n\n"
    for istep in range(1, len(pipeline_steps)+1):
        description += f" {istep}) {pipeline_steps[istep][1]}\n"

    epilog = textwrap.dedent("""
    Examples:
       $ python marvelpipe.py config.yaml              - Run all pipeline steps on the files specified in config.yaml
       $ python marvelpipe.py -f 3 config.yaml         - Run pipeline steps 3, 4, ... until the end
       $ python marvelpipe.py -f 3 -l 6 config.yaml    - Run pipeline steps 3, 4, 5, and 6
    """)


    # The following class is a hacky way to use multiple formatters in argparse

    class CustomFormatter(ap.ArgumentDefaultsHelpFormatter, ap.RawDescriptionHelpFormatter):
        pass


    # Set up the command line argument parser

    parser = ap.ArgumentParser(prog='marvelpipe.py', formatter_class=CustomFormatter, description=description, epilog=epilog)
    parser.add_argument('configfile', help='YAML configuration file')  
    parser.add_argument('-f', '--first', default=1, type=int, required=False, help='First pipeline step to execute')
    parser.add_argument('-l', '--last', default=len(pipeline_steps), type=int, required=False, help='Last pipeline step to execute')

    args = parser.parse_args()
    if args.first < 0:
        args.first = 0
    if args.last > len(pipeline_steps):
        args.last = len(pipeline_steps)


    # Load the yaml file with the pipeline configuration into a dictionary

    config = yaml.safe_load(open(args.configfile))
    config['origin'] = args.configfile


    # Report the number of raw files found

    num_raw_bias = len(config['MasterBiasImage']['inputPath'])
    num_raw_dark = len(config['MasterDarkImage']['inputPath'])
    num_raw_flat = len(config['MasterFlatImage']['inputPath'])
    num_raw_thar = len(config['MasterThArImage']['inputPath'])
    num_raw_science = len(config['BiasAndDarkSubtractedScienceImage']['inputPath'])
   
   
    # Run the different pipeline steps

    print("MARVEL Data Reduction Pipeline - {0}".format(datetime.now()))
    print(f"\n{args.configfile} contains the paths to: ")
    print(f" - {num_raw_bias} raw bias images")
    print(f" - {num_raw_dark} raw dark images")
    print(f" - {num_raw_flat} raw flat images")
    print(f" - {num_raw_thar} raw ThAr images")
    print(f" - {num_raw_science} raw science images")
    print("\nRunning pipeline:")

    for istep in range(args.first, args.last+1):
        pipeline_step, description = pipeline_steps[istep]

        print("Step {0}: {1}".format(istep, description).ljust(50, ' '), end=' ')
        t1 = time.time()
        pipeline_step(config)
        t2 = time.time()
        print("[{:.1f} s]".format(t2-t1))


    # That's it!

    print("Ready.")




