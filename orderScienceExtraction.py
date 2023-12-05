import yaml
from orderExtraction import OrderExtraction
from database        import DatabaseFromLocalFile


class OrderScienceExtraction(OrderExtraction):
    """
    Class that performs the extraction of the orders given an science image and an order mask.
    """
    def __init__(self, database=None, debug=0, **scienceAndMaskHash):
        super().__init__(database=database, debug=debug, **scienceAndMaskHash)
    





if __name__ == "__main__":
    o_params = yaml.safe_load(open("params.yaml"))["orderMaskImage"]
    s_params = yaml.safe_load(open("params.yaml"))["rawScienceImage"]
    e_params = yaml.safe_load(open("params.yaml"))["orderScienceExtraction"]

    databaseName = "pipelineDatabase.txt"
    print("Creating a local database file with the name: ", databaseName)

    db = DatabaseFromLocalFile(databaseName)
    print(" ")
    print("X")
    # Extracted Science orderders
    calibrated_science_path = s_params["outpath"]
    order_mask_path = o_params["outpath"]

    scienceExtractor = OrderScienceExtraction(db, debug=1, ExtractedOrders= order_mask_path, ScienceImages=calibrated_science_path)

    scienceExtractor.run(e_params["outpath"])
