# We use this module for opperations that will be commonly used in different pipelineModules




from astropy.io import fits
import numpy as np
from numba import njit
from pymongo import MongoClient
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

    images     = {"Master Dark Image": "DarkImages", "Master Bias Image": "BiasImages", "Master Flat Image": "FlatImages", "Calibrated Science Image" : "ScienceImages", "Extracted Flat Orders" : "ExtractedOrders"}
    typeImage  = dict["type"]
    collection = db[images[typeImage]]

    isInCollection = np.all([x == dict for x in collection.find({"_id" : dict["_id"]})])

    if not isInCollection:
        collection.inster_one(dict)

    elif overWrite:
        collection.delete_many({"_id": dict["_id"]})
        collection.insert_one(dict)
    else:
        print("Document is already in database")
    
    
        
def getPositionsOfOrders():
    p = {}
    positions = [ 276.47034757,  293.29151457,  310.02540273,  326.69533521,  343.33412994,
                  544.81733451,  561.52993913,  578.16816779,  594.77397072,  611.34838329,
                  805.89583395,  822.52396273,  839.08604478,  855.64288321,  872.15142763,
                  1060.03078413, 1076.57288292, 1093.10222727, 1109.57997859, 1126.06574561,
                  1307.49950152, 1323.98236364, 1340.45814225, 1356.88917131, 1373.33886427,
                  1548.57285103, 1565.00483684, 1581.43936976, 1597.82852469, 1614.2276996,
                  1783.50697697, 1799.88069926, 1816.27169165, 1832.62626626, 1848.97029067,
                  2012.4966333,  2028.86548177, 2045.19293021, 2061.49400573, 2077.84073522,
                  2235.81530904, 2252.09805389, 2268.42266958, 2284.69600658, 2300.96612213,
                  2453.63432501, 2469.89739046, 2486.13183756, 2502.35643693, 2518.64524964,
                  2666.12911951, 2682.34012918, 2698.55284353, 2714.77374492, 2731.00008524,
                  2873.48056255, 2889.69724295, 2905.8721813,  2922.03283324, 2938.22340887,
                  3075.92275674, 3092.05623548, 3108.18211758, 3124.33547446, 3140.5228562,
                  3273.55744873, 3289.65985127, 3305.74587484, 3321.84634763, 3337.97406045,
                  3466.51604359, 3482.60660853, 3498.68261805, 3514.76529123, 3530.86948563,
                  3655.01399839, 3671.0756758,  3687.10863157, 3703.16053125, 3719.23410192,
                  3839.18975558, 3855.18577666, 3871.18428103, 3887.19504057, 3903.2303084,
                  4019.15206183, 4035.11446118, 4051.07855914, 4067.05520655, 4083.0631696,
                  4194.99493977, 4210.90911916, 4226.81361559, 4242.76181204, 4258.74556228,
                  4366.85265715, 4382.74579881, 4398.65075966, 4414.56654042, 4430.52156323,
                  4534.87393071, 4550.73673709, 4566.61157353, 4582.49843329, 4598.42311403,
                  4699.18752196, 4714.98358156, 4730.80033616, 4746.66528973, 4762.56499905,
                  4859.78583132, 4875.59348888, 4891.40140799, 4907.22761904, 4923.09825375,
                  5016.93459551, 5032.64549592, 5048.42244384, 5064.23014409, 5080.07679083,
                  5170.58206666, 5186.28637915, 5202.03602523, 5217.80440809, 5233.6164422,
                  5320.87026978, 5336.60357033, 5352.29338351, 5367.99415323, 5383.802532,
                  5467.93062264, 5483.6007161,  5499.29484865, 5514.9938392,  5530.71448992,
                  5611.8116421,  5627.45238863, 5643.10804203, 5658.79428705, 5674.51066553,
                  5752.57157529, 5768.24083532, 5783.83436648, 5799.5038061,  5815.20186385,
                  5890.31257984, 5905.99104177, 5921.57994397, 5937.19628705, 5952.87176223,
                  6025.1655429,  6040.78130841, 6056.39616071, 6071.93147197, 6087.60220031,
                  6157.13141789, 6172.75301585, 6188.29018591, 6203.88335827, 6219.43799964,
                  6286.32141102, 6301.85536944, 6317.44826399, 6332.98769397, 6348.56123187,
                  6412.77928345, 6428.31489145, 6443.86678921, 6459.38586516, 6474.93177754,
                  6536.55262145, 6552.13287731, 6567.59370386, 6583.14381099, 6598.61357654,
                  6657.78972455, 6673.34019849, 6688.79435223, 6704.30936524, 6719.77348689,
                  6776.52748189, 6792.00255091, 6807.48250133, 6822.93459153, 6838.39781376,
                  6892.7750038,  6908.23332428, 6923.67533133, 6939.12268669, 6954.51858596,
                  7006.61457325, 7022.07888762, 7037.48679237, 7052.87933593, 7068.32488768,
                  7118.13196512, 7133.57385512, 7148.93637816, 7164.36320259, 7179.73217468,
                  7227.38126648, 7242.79661041, 7258.13947084, 7273.53196985, 7288.88080687,
                  7334.39581257, 7349.80247333, 7365.14435328, 7380.4447697,  7395.81451559,
                  7439.30037827, 7454.60445026, 7469.97195087, 7485.28355235, 7500.55852004,
                  7542.07353484, 7557.39063506, 7572.67393223, 7587.9620371,  7603.27723802,
                  7642.79516329, 7658.06814181, 7673.39479634, 7688.64171006, 7703.89948757,
                  7741.58635502, 7756.83924692, 7772.06161255, 7787.32189989, 7802.61517515,
                  7838.40913602, 7853.66721849, 7868.88051201, 7884.10069845, 7899.34726358,
                  7933.36761135, 7948.61201168, 7963.82911585, 7979.06078845, 7994.28182282,
                  8026.54161756, 8041.78600585, 8056.98730925, 8072.1994811,  8087.4049997,
                  8118.05818416, 8133.21278398, 8148.39855214, 8163.58961481, 8178.83446714,
                  8207.80506012, 8223.0082302,  8238.21146775, 8253.38975396, 8268.55656278,
                  8296.03337131, 8311.23239918, 8326.40316335, 8341.56359151, 8356.70845558,
                  8382.72735235, 8397.87742691, 8413.09251235, 8428.25387445, 8443.38975797,
                  8468.01982291, 8483.20877179, 8498.35094508, 8513.48292909, 8528.59482817,
                  8551.94211933, 8567.09489802, 8582.24583715, 8597.3696133,  8612.4598388,
                  8634.66100577, 8649.77149577, 8664.86860776, 8679.95026871, 8695.01702119,
                  8716.1257905,  8731.26146683, 8746.36431981, 8761.43444258, 8776.51094225,
                  8796.6188383,  8811.70113776, 8826.7575624,  8841.81430576, 8856.88565171,
                  8876.06831061, 8891.1168787,  8906.16325439, 8921.23341483, 8936.33660725,
                  8954.74859894, 8969.79446379, 8984.85981612, 8999.91775128, 9014.97026853,
                  9032.724891,   9047.78518792, 9062.84472106, 9077.89539614, 9092.92734815,
                  9110.16435628, 9125.21649696, 9140.2554798,  9155.28750395, 9170.29639876,
                  9187.29014879, 9202.32578198, 9217.34812559, 9232.37145401, 9247.35445645,
                  9264.2826628,  9279.29955006, 9294.30297299, 9309.29109291, 9324.27277666,
                  9341.33600799, 9356.32889281, 9371.31776225, 9386.29471818, 9401.26896406,
                  9418.76224184, 9433.74115057, 9448.72054669, 9463.68950458, 9478.64371499]


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
        return table.data["Flux"]




if __name__ == "__main__":
    path = os.getcwd() + "/Data/ProcessedData/ExtractedOrders/extracted_flat_orders.fits"

    print(type(getExtractedPosition(path, 3, 2)[:5]))
    print(type(getExtractedFlux(path, 3, 2)[:5]))
