# We use this module for opperations that will be commonly used in different pipelineModules




from astropy.io import fits
import numpy as np
from numba import njit
from pymongo import MongoClient
import matplotlib.pyplot as plt
import os

client = MongoClient()
db = client["databaseMARVEL"]




def saveFigure(figure):
    ...

def getImages(paths):
    return np.array([getImage(path) for path in paths])


def getImage(path):

    hdul = fits.open(path)
    return hdul[0].data

@njit()
def getGain(path):
    # TODO: This should be read from the fits file [e/ADU]
    return 9.4

@njit()
def getExposureTime(path):
    # TODO: This should be read from the fits file [s]
    return 5

@njit()
def getDarkRate(path):
    # TODO: This should be read from the fits file [e/pix/s]
    return 0.1




def addToDataBase(dict, overWrite=False):
    # We should have proper error catching for:
    # 1. We want to add somthing that already exist
    # 2. image is not in the dictornary
    # 3. Check that dict is in the right format 

    images     = {"Master Dark Image": "DarkImages", "Master Bias Image": "BiasImages", "Master Flat Image": "FlatImages", "Calibrated Science Image" : "ScienceImages", "Extracted Flat Orders" : "ExtractedOrders", "Extracted Science Orders" : "ExtractedOrders", "Optimal Extracted" : "OptimalExtracted"}
    typeImage  = dict["type"]

    collection = db[images[typeImage]]
    isInCollection = np.all([x == dict for x in collection.find({"_id" : dict["_id"]})])
    
    if not isInCollection:
        collection.insert_one(dict)

    elif overWrite:
        collection.delete_many({"_id": dict["_id"]})
        collection.insert_one(dict)
    else:
        print("Document is already in database")
    
    
        
def getPositionsOfOrders():
    p = {}

    positions = [ 291.32511331,  308.13650784,  324.87281272,  341.54605083,  358.18362429,
                  559.03889327,  575.74417217,  592.37975593,  608.98547766,  625.55506675,
                  819.5088899,   836.13313052,  852.68787197,  869.24913251,  885.74779306,
                  1073.06227972, 1089.60280963, 1106.13686539, 1122.60538217, 1139.09384015,
                  1319.98270952, 1336.45545127, 1352.93951101, 1369.3564903,  1385.81517755,
                  1560.52638844, 1576.94800396, 1593.38875425, 1609.77244001, 1626.17138969,
                  1794.95252285, 1811.31919088, 1827.71412935, 1844.06152865, 1860.3943134,
                  2023.45313782, 2039.82477921, 2056.14846741, 2072.43885688, 2088.79336825,
                  2246.30687052, 2262.58037622, 2278.91678608, 2295.18077681, 2311.44475605,
                  2463.68339985, 2479.94342109, 2496.17060916, 2512.38585697, 2528.68716431,
                  2675.74486456, 2691.94827731, 2708.16345344, 2724.38677553, 2740.6114645,
                  2882.67699209, 2898.9007074,  2915.07395269, 2931.23289217, 2947.41607895,
                  3084.72970505, 3100.85724657, 3116.97723705, 3133.13494815, 3149.32390321,
                  3281.97914122, 3298.08205422, 3314.1629092,  3330.26008735, 3346.38386363,
                  3474.56185421, 3490.65820808, 3506.7381852,  3522.81868878, 3538.9249549,
                  3662.70552415, 3678.77373716, 3694.80639259, 3710.85801125, 3726.93442405,
                  3846.54666403, 3862.54197817, 3878.53952323, 3894.55011113, 3910.58375253,
                  4026.18486199, 4042.14818234, 4058.11165534, 4074.08690924, 4090.09592736,
                  4201.70833144, 4217.62104735, 4233.5177607,  4249.46852503, 4265.45257841,
                  4373.25239469, 4389.14943174, 4405.05712164, 4420.97360096, 4436.92802654,
                  4540.98100505, 4556.84751486, 4572.72199251, 4588.60908269, 4604.53436854,
                  4705.01444015, 4720.8082634,  4736.62459244, 4752.4880146,  4768.39067691,
                  4865.3307715,  4881.14375678, 4896.95297952, 4912.77877029, 4928.64694428,
                  5022.21751492, 5037.92289016, 5053.70551834, 5069.51401269, 5085.36091964,
                  5175.6070596,  5191.31337982, 5207.06401577, 5222.83264602, 5238.63959574,
                  5325.64776186, 5341.38160841, 5357.07321854, 5372.76975286, 5388.57386883,
                  5472.46482913, 5488.14002709, 5503.83473178, 5519.53010318, 5535.24595125,
                  5616.11354213, 5631.7572275,  5647.41483106, 5663.09664765, 5678.81028736,
                  5756.65016822, 5772.32157135, 5787.91544015, 5803.58346472, 5819.27445489,
                  5894.17489377, 5909.85494946, 5925.44617871, 5941.0574191,  5956.7291757,
                  6028.81900749, 6044.43669469, 6060.05276183, 6075.58537882, 6091.25112,
                  6160.58307642, 6176.20623045, 6191.74527856, 6207.33482006, 6222.88397367,
                  6289.57948682, 6305.11476372, 6320.70853578, 6336.24307511, 6351.81291397,
                  6415.84893126, 6431.3871914,  6446.93767603, 6462.45441488, 6477.99468634,
                  6539.44229845, 6555.0226688,  6570.48173089, 6586.02912417, 6601.49415723,
                  6660.5037988,  6676.05439693, 6691.50603125, 6707.01622514, 6722.47797038,
                  6779.07210805, 6794.54825914, 6810.02235917, 6825.47093751, 6840.93318094,
                  6895.15593203, 6910.61277381, 6926.0512185,  6941.49618333, 6956.88721257,
                  7008.83670378, 7024.29865431, 7039.70323453, 7055.09173619, 7070.53196713,
                  7120.201617,   7135.64025293, 7151.00020766, 7166.4203243,  7181.78669318,
                  7229.30196155, 7244.71354667, 7260.05258432, 7275.44006805, 7290.78721593,
                  7336.17227977, 7351.57518821, 7366.91246493, 7382.21054671, 7397.57631266,
                  7440.93512107, 7456.23743733, 7471.60090711, 7486.90999963, 7502.18336319,
                  7543.57421032, 7558.88699438, 7574.16904855, 7589.45578236, 7604.76631993,
                  7644.16651062, 7659.43420465, 7674.75918478, 7690.00417865, 7705.26048406,
                  7742.83274924, 7758.08162013, 7773.30337396, 7788.55792306, 7803.84975333,
                  7839.53370014, 7854.78855974, 7869.99912179, 7885.21651343, 7900.45763509,
                  7934.37316052, 7949.61367344, 7964.83028168, 7980.05625639, 7995.27477325,
                  8027.43341335, 8042.67661918, 8057.87311515, 8073.08064617, 8088.28450734,
                  8118.83958934, 8133.99108617, 8149.17299994, 8164.36248544, 8179.60391042,
                  8208.48028621, 8223.67987352, 8238.87918076, 8254.05588821, 8269.22080374,
                  8296.60738944, 8311.80081436, 8326.96972418, 8342.12795053, 8357.27254234,
                  8383.20080873, 8398.35059084, 8413.5610607,  8428.72084649, 8443.85596185,
                  8468.40043133, 8483.58530475, 8498.72627701, 8513.85696893, 8528.96797479,
                  8552.23059812, 8567.38141239, 8582.53340123, 8597.65561739, 8612.74596354,
                  8634.86271825, 8649.97247074, 8665.07142738, 8680.15389707, 8695.22124246,
                  8716.24857006, 8731.38082028, 8746.4840728,  8761.5555527,  8776.6298448,
                  8796.66171709, 8811.74327634, 8826.80184025, 8841.85725287, 8856.92769545,
                  8876.04108954, 8891.08956775, 8906.13550036, 8921.20256057, 8936.30465761,
                  8954.64942694, 8969.69611155, 8984.76011015, 8999.8164109,  9014.86796454,
                  9032.56039406, 9047.6197992,  9062.6787886,  9077.72916908, 9092.75885384,
                  9109.93892114, 9124.99229017, 9140.03326311, 9155.06329463, 9170.07176993,
                  9187.01183413, 9202.04895811, 9217.07383914, 9232.09589269, 9247.07872144,
                  9263.95951302, 9278.97679476, 9293.97823002, 9308.96462467, 9323.94690983,
                  9340.9710289,  9355.96298507, 9370.94862026, 9385.92533157, 9400.89531516,
                  9418.35583235, 9433.33273106, 9448.30956075, 9463.27443336, 9478.23426121]


    fibers = []
    orders = []

    for i in np.arange(np.size(positions)):
        order = int((i)/5)+1
        fiber = int((i+1)%5)
        if (fiber == 0):
            fiber=5

        fibers.append(fiber)
        orders.append(order)

    p = np.array([fibers, orders, positions])

    return p








def checkInputExtractedData(path, order, fiber):

    """
    Tis function checks that for the functions getExtractedPosition and getExtractedFlux the input
    that is specified is sensible. If this is the case the function returns the specified table, if
    not the function return None.
    """

    # Check that path exist
    if not os.path.isfile(path):
        print("Error: path does not exist")
        return

    # Check that type of fits file is Extracted Flux
    hdul = fits.open(path)

    fileType = hdul[0].header["type"]

    if not (("Extracted" in fileType) and ("Orders" in fileType)):
        print("Error: filetype {} is not a type of Extracted Orders".format(fileType))
        return

    # Try to find order/fiber:
    # First try find correct table at expect location
    idx = (order-1)*5 + fiber

    if ((hdul[idx].header["order"] == order) and (hdul[idx].header["fiber"] == fiber)):
        table = hdul[idx]
    else:
    # If the correct table is not at the expected location, we loop over every order
    # and check if the correct table among them.
        locationFound = False
        for idx in np.arange(np.size(hdul)):
            try:
                correctLocation = ((hdul[idx].header["order"] == order) and (hdul[idx].header["fiber"] == fiber))
            except:
                continue
            if correctLocation:
                locationFound = True
                table = hdul[idx]
                break
        if not locationFound:
            print("order {o} and fiber {f} not found in file {file}.".format(o=order, f=fiber, file=path))
            return
    return table


def getExtractedPosition(path, order, fiber):

    # First we check that the input that is provided is sensible
    table = checkInputExtractedData(path, order, fiber)

    # If it is sensible, we return the Positions
    if table is None:
        return
    else:
        return np.array([ (x, y) for x, y in zip(table.data['X'], table.data['Y'])], dtype=np.int16)




def getExtractedFlux(path, order, fiber):

    # First we check that the input that is provided is sensible
    table = checkInputExtractedData(path, order, fiber)

    # If it is sensible, we return the Flux
    if table is None:
        return
    else:
        return (table.data["Flux"]).astype(np.float64)

def getFibersAndOrders(path):
    """
    This function returns the fibers and orders that are in the extracted data. 
    INPUT: path of the extracted flux fits file 
    OUTPUT: list of fibers, list of fibers 
    """
    # Check that path exist
    if not os.path.isfile(path):
        print("Error: path does not exist")
        return

    # Check that type of fits file is Extracted Flux
    hdul = fits.open(path)

    fileType = hdul[0].header["type"]

    if not (("Extracted" in fileType) and ("Orders" in fileType)):
        print("Error: filetype {} is not a type of Extracted Orders".format(fileType))
        return

    orders = hdul[0].header["orders"]
    orders = [int(i) for i in orders[1:-1].split(", ")]

    fibers = hdul[0].header["fibers"]
    fibers = [int(i) for i in fibers[1:-1].split(", ")]

    return fibers, orders



def getAllExtractedInfo(path):
    """ 
    This function extracts all the info (positions and flux) for all the orders from 
    extracted images. 
    INPUT: path the the extracted images.
    OUPUT: dictionary with keys the orders and values a dictorionary with fibers/(positions, flux) 
           as key/values. 
    REMARK: The idea behind this function is that we only open the file one time and extract all the 
            information at once instead of opening and closing the same file for different orders. 
    """

    fibers, orders = getFibersAndOrders(path)

    # Check that path exist
    if not os.path.isfile(path):
        print("Error: path does not exist")
        return

    # Check that type of fits file is Extracted Flux
    hdul = fits.open(path)

    fileType = hdul[0].header["type"]

    if not (("Extracted" in fileType) and ("Orders" in fileType)):
        print("Error: filetype {} is not a type of Extracted Orders".format(fileType))
        return



    def getTable(order, fiber):
        # Try to find order/fiber:
        # First try find correct table at expect location
        idx = (order-1)*5 + fiber

        if ((hdul[idx].header["order"] == order) and (hdul[idx].header["fiber"] == fiber)):
            table = hdul[idx]

        else:
            # If the correct table is not at the expected location, we loop over every order
            # and check if the correct table among them.
            locationFound = False
            for idx in np.arange(np.size(hdul)):
                try:
                    correctLocation = ((hdul[idx].header["order"] == order) and (hdul[idx].header["fiber"] == fiber))
                except:
                    continue
                if correctLocation:
                    locationFound = True
                    table = hdul[idx]
                    break
            if not locationFound:
                print("order {o} and fiber {f} not found in file {file}.".format(o=order, f=fiber, file=path))
                return
        return table

    fluxes = {}

    for o in orders:
        for f in fibers:
            table = getTable(o, f)
            if table is None:
                print("Error in file {}. Not correct format.".format(path))
                return 
            flux = (table.data["Flux"]).astype(np.float64)
            xPos = (table.data['X']).astype(np.int16)
            yPos = (table.data['Y']).astype(np.int16)

            if o in fluxes:
                fluxes[o].update({f : (xPos, yPos, flux)})
            else:
                fluxes[o] = {f : (xPos, yPos, flux)}

    return fluxes

    




if __name__ == "__main__":
    path = os.getcwd() + "/Data/ProcessedData/ExtractedOrders/extracted_flat_orders.fits"
    getAllExtractedInfo(path)


    # path_flat = os.getcwd() + "/Data/ProcessedData/MasterFlat/master_flat.fits"
    # full_image = getImage(path_flat)

    # positions = getExtractedPosition(path, 20, 4)
    # flux      = getExtractedFlux(path, 20, 4)

    # xx, yy = zip(*positions)
    # ofst = np.min(xx)
    # mask = [ True if x in np.arange(4440 + ofst, 4460 + ofst) else False for x in xx]

    # x_mask = np.array(xx)[mask]
    # y_mask = np.array(yy)[mask]

    # x0 = np.min(x_mask)
    # y0 = np.min(y_mask)

    
        
    # buffr = 20

    # print(np.max(x_mask))
    # print(x0)
    # print("==========")
    # print(np.max(y_mask))
    # print(y0)
    # print("==========")
    # print(np.min(yy), np.max(yy))

    # image = np.zeros((np.max(x_mask)-x0+1, np.max(y_mask)-y0+1))
    # mask_image = np.empty((np.max(x_mask)-x0+1, np.max(y_mask)-y0+1+2*buffr), dtype=bool)
    # mask_image.fill(False)

    # buffer1 = np.zeros((np.max(x_mask)-x0-1, buffr))
    # buffer2 = np.zeros((np.max(x_mask)-x0-1, buffr))
    
    # for x, y in zip(x_mask, y_mask):
    #     mask_image[x-x0, y-y0+buffr] = True

    # image = full_image[x0:np.max(x_mask)+1, y0-buffr:np.max(y_mask)+buffr+1]
    # extracted_img = image[mask_image]

    
    # x_max = np.max(y_mask)+2*buffr + 1 - y0 
    # cm = plt.cm.get_cmap("RdYlBu")
    
    # plt.plot(np.transpose(image))
    
    # plt.fill_between(np.arange(buffr), 30000, 0, color='grey')
    # plt.fill_between(np.arange(x_max-buffr, x_max), 30000, 0, color='grey')
    # plt.show()

    # plt.imshow(np.transpose(image))
    # plt.show()

        
    # a = int(len(image[0][buffr:-buffr])/2)+1


    # plt.plot([np.mean(x[buffr:-buffr]) for x in image])
    # plt.show()
