# This is an example script on how to run the different MARVEL pipeline components
# on a given set of data input files.

from database import DatabaseFromLocalFile
from pipeline import MasterBias, MasterFlat, CalibratedScienceFrames 
from orderMaskExtraction import OrderMaskExtraction 
from orderExtraction import OrderExtraction
from optimalExtraction import OptimalExtraction



# All information of pipeline component results will be stored both in the FITS header
# as well as in a database, so that searching becomes easier. In this example we will
# use a local ascii file as a database proxy. 

db = DatabaseFromLocalFile("pipelineDatabase.txt")
db.save()
print("")


# Create a Master Bias Image

raw_bias_paths = ['Data/RawData/CalibrationImages/Bias/bias_0001.fits',
                  'Data/RawData/CalibrationImages/Bias/bias_0002.fits',
                  'Data/RawData/CalibrationImages/Bias/bias_0003.fits',
                  'Data/RawData/CalibrationImages/Bias/bias_0004.fits',
                  'Data/RawData/CalibrationImages/Bias/bias_0005.fits',
                  'Data/RawData/CalibrationImages/Bias/bias_0006.fits',
                  'Data/RawData/CalibrationImages/Bias/bias_0007.fits',
                  'Data/RawData/CalibrationImages/Bias/bias_0008.fits']

masterBias = MasterBias(db, BiasImages=raw_bias_paths)
masterBias.run("masterBias.fits")
print("")



# Create a Master Flat Image

raw_flat_paths = ['Data/RawData/CalibrationImages/Flat/flat_0001.fits',
                  'Data/RawData/CalibrationImages/Flat/flat_0002.fits',
                  'Data/RawData/CalibrationImages/Flat/flat_0003.fits',
                  'Data/RawData/CalibrationImages/Flat/flat_0004.fits',
                  'Data/RawData/CalibrationImages/Flat/flat_0005.fits']
master_bias_path = "Data/ProcessedData/MasterBias/masterBias.fits"

masterFlat = MasterFlat(db, FlatImages=raw_flat_paths, BiasImages=master_bias_path)
masterFlat.run("masterFlat.fits")



# Correct a science image for the bias. We call the output a "calibrated science image". 

raw_science_path = "Data/RawData/ScienceFrames/science_0001.fits"
master_bias_path = "Data/ProcessedData/MasterBias/masterBias.fits"

calibration = CalibratedScienceFrames(db, ScienceImages=raw_science_path,
                                          BiasImages=master_bias_path)
calibration.run("biasCorrectedScience.fits")
print(" ")



# Using a master flat field, determine the pixel masks for each order.  

masterflat_path = "Data/ProcessedData/MasterFlat/masterFlat.fits"
maskExtractor = OrderMaskExtraction(db, debug=1, FlatImages=masterflat_path)
maskExtractor.run("orderMask.fits")
print(" ")



# Using the pixel mask, extract the orders of a science image

science_image_path = "Data/ProcessedData/CalibratedScience/biasCorrectedScience.fits"
order_mask_path    = "Data/ProcessedData/ExtractedOrders/orderMask.fits"
orderExtractor = OrderExtraction(db, debug=1, ExtractedOrders=order_mask_path,
                                              ScienceImages=science_image_path)
orderExtractor.run("extractedOrdersScience.fits")
print("")



# Extract a 1D spectrum for each order, not yet wavelength calibrated but corrected for
# the instrumental response curve. We call this an "optimal science extraction".

extractedSciencePath = "Data/ProcessedData/ExtractedOrders/extractedOrdersScience.fits"
extractedFlatPath = "Data/ProcessedData/ExtractedOrders/orderMask.fits"
optimalScience = OptimalExtraction(db, debug=1, ExtractedOrders=[extractedSciencePath, extractedFlatPath])
optimalScience.run("optimal_extracted_science_flux.fits")
print("")


# Save all the information of the pipeline process to the database. 

db.save()
