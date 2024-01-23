# This is an example script on how to run the different MARVEL
# pipeline components on a given set of data input files.

import os
from database import DatabaseFromLocalFile
from pipeline import MasterBias, MasterFlat, BiasCorrectedScienceFrames
from orderMaskExtraction import OrderMaskExtraction
from orderExtraction import OrderExtraction
from optimalExtraction import OptimalExtraction
from wavelengthCalibration import WavelengthCalibration


# All information of pipeline component results will be stored both in
# the FITS header as well as in a database, so that searching becomes
# easier. In this example we will use a local ascii file as a database
# proxy.

db = DatabaseFromLocalFile("pipelineDatabase.txt")
db.save()
print("")


# Create a Master Bias Image
raw_bias_root = "Data/RawData/CalibrationImages/Bias"
raw_bias_paths = [raw_bias_root + '/20230504T135640_BBBBB_H_0000.fits',
                  raw_bias_root + '/20230504T135642_BBBBB_H_0000.fits',
                  raw_bias_root + '/20230504T135644_BBBBB_H_0000.fits',
                  raw_bias_root + '/20230504T135648_BBBBB_H_0000.fits',
                  raw_bias_root + '/20230504T135652_BBBBB_H_0000.fits']

masterBias = MasterBias(db, BiasImages=raw_bias_paths)
masterBias.run("masterBias.fits")
print("")


# Create a Master Flat Image
raw_flat_root = 'Data/RawData/CalibrationImages/Flat'
raw_flat_paths = [raw_flat_root + '/20230504T142007_FFFFF_H_0001.fits',
                  raw_flat_root + '/20230504T142402_FFFFF_H_0001.fits',
                  raw_flat_root + '/20230504T142759_FFFFF_H_0001.fits',
                  raw_flat_root + '/20230504T143157_FFFFF_H_0001.fits',
                  raw_flat_root + '/20230504T143554_FFFFF_H_0001.fits']

master_bias_path = "Data/ProcessedData/MasterBias/masterBias.fits"
masterFlat = MasterFlat(db, FlatImages=raw_flat_paths,
                        BiasImages=master_bias_path)
masterFlat.run("masterFlat.fits")


# Correct science images for the bias. We call the output
# a "calibrated science image".
raw_science_root = 'Data/RawData/ScienceFrames'
raw_science_paths = [raw_science_root + "/20230507T075752_ESSSS_H_0600.fits",
                     raw_science_root + "/20230507T080243_ESSSS_H_0600.fits"]
master_bias_path = "Data/ProcessedData/MasterBias/masterBias.fits"

for raw_science_path in raw_science_paths:
    calibration = BiasCorrectedScienceFrames(db,
                                             ScienceImages=raw_science_path,
                                             BiasImages=master_bias_path)
    basename, extension = os.path.splitext(os.path.basename(raw_science_path))
    outputPath = basename + "_BC" + extension
    calibration.run(outputPath)

print(" ")


# Using a master flat field, determine the pixel masks for each order.

masterflat_path = "Data/ProcessedData/MasterFlat/masterFlat.fits"
maskExtractor = OrderMaskExtraction(db, debug=1, FlatImages=masterflat_path)
maskExtractor.run("orderMask.fits")
print(" ")


# Using the pixel mask, extract the orders ("EO") of a
# (bias corrected) science image
bias_cor_sc_root = "Data/ProcessedData/BiasCorrectedScience"
science_image_paths = [bias_cor_sc_root + "/20230507T075752_ESSSS_H_0600_BC.fits",
                       bias_cor_sc_root + "/20230507T080243_ESSSS_H_0600_BC.fits"]

order_mask_path = "Data/ProcessedData/ExtractedOrders/orderMask.fits"
for science_image_path in science_image_paths:
    orderExtractor = OrderExtraction(db, debug=1,
                                     ExtractedOrders=order_mask_path,
                                     ScienceImages=science_image_path)
    basename, extension = os.path.splitext(
        os.path.basename(science_image_path))
    outputPath = basename[:-3] + "_EO" + extension
    orderExtractor.run(outputPath)
print("")


# Extract a 1D spectrum for each order, not yet wavelength calibrated but corrected for
# the instrumental response curve. We call this an "optimal science extraction" ("OP").

extracted_science_paths = ["Data/ProcessedData/ExtractedOrders/20230507T075752_ESSSS_H_0600_EO.fits",
                           "Data/ProcessedData/ExtractedOrders/20230507T080243_ESSSS_H_0600_EO.fits"]
extracted_flat_path = "Data/ProcessedData/ExtractedOrders/orderMask.fits"
for extracted_science_path in extracted_science_paths:
    optimalScience = OptimalExtraction(db, debug=1, ExtractedOrders=[
        extracted_science_path, extracted_flat_path])
    basename, extension = os.path.splitext(os.path.basename(extracted_science_path))
    outputPath = basename[:-3] + "_OP" + extension
    optimalScience.run(outputPath)
print("")


# Do the wavelength calibration

optimal_extracted_science_paths = ["Data/ProcessedData/OptimalExtraction/20230507T075752_ESSSS_H_0600_OP.fits",
                                   "Data/ProcessedData/OptimalExtraction/20230507T080243_ESSSS_H_0600_OP.fits"]
for optimal_extracted_science_path in optimal_extracted_science_paths:
    wavelength_calibration = WavelengthCalibration(db, debug=1, OptimalExtracted=optimal_extracted_science_path)
    basename, extension = os.path.splitext(os.path.basename(optimal_extracted_science_path))
    outputPath = basename[:-3] + "_WC" + extension
    wavelength_calibration.run(outputPath)
print("")


# Save all the information of the pipeline process to the database.

db.save()
