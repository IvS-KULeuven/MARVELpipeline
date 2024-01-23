import yaml
import time
from orderExtraction import OrderExtraction
from database        import DatabaseFromLocalFile


class OrderScienceExtraction(OrderExtraction):
    """
    Class that performs the extraction of the orders given an science image and an order mask.
    """
    def __init__(self, database=None, debug=0, **scienceAndMaskHash):
        super().__init__(database=database, debug=debug, **scienceAndMaskHash)
    





if __name__ == "__main__":

    t1 = time.time()
    o_params = yaml.safe_load(open("params.yaml"))["orderMaskImage"]
    s_params = yaml.safe_load(open("params.yaml"))["rawScienceImage"]
    e_params = yaml.safe_load(open("params.yaml"))["orderScienceExtraction"]

    databaseName = "pipelineDatabase.txt"
    db = DatabaseFromLocalFile(databaseName)

    # Extracted Science orders
    calibrated_science_path = s_params["outpath"]
    order_mask_path = o_params["outpath"]

    scienceExtractor = OrderScienceExtraction(db, debug=0, ExtractedOrders= order_mask_path, ScienceImages=calibrated_science_path)
    scienceExtractor.run(e_params["outpath"])
    db.save()

    t2 = time.time()
    print(f"This took: {t2-t1}")
