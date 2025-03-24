import numpy as np
from pathlib import Path
from astropy.io import fits
from astropy.io.fits import Header
import tools



class MasterFlatOneDimOrders:
    """
    Class that extracts the 1D flat fields order per order, for each fiber. 
    The product of this class is for monitoring purposes, and is not used anywhere else in the pipeline.
    """

    def __init__(self):
        pass



    def run(self, master_flat_path, twodim_mask_path, output_path):
        """
        Extract the one-dimensional flatfield curve for each order and each fiber
        This flatfield curve still contains the blaze function.

        Input:
            master_flat_path: path to the 2D master flatfield FITS file
            twodim_mask_path: path to the FITS file containing the mask boundaries of each order/fiber,
                              as derived from the master flat.
            output_path: path to save the result to. 

        Output:
            One dimensional flatfield curve for each order and each fiber written to a FITS file
        """

        # Get the master flat

        masterFlat = tools.getImage(master_flat_path)

        # Get the "ridge" of maximum flatfield values for each order/fiber

        mask_fits_file = fits.open(twodim_mask_path)
        ridge_masks = mask_fits_file[0].data 

        # Loop over all masks (one per order / fiber) and save the 1D flatfield profiles in a FITS file

        hdu_list = fits.HDUList([fits.PrimaryHDU()])       # Compulsory but empty primary HDU

        for imask in range(len(ridge_masks)):

            fiber = 5 - imask % 5                          # Fiber 1 is the calibration fiber, 2,3,4,5 are the science fibers
            order = 98 -  int(np.floor(imask / 5.0))       # Ranging from 98 (blue) down to 33 (red). 

            current_mask = ridge_masks[imask]
            xpixel = []
            ridgeValues = []                        # Ridge of the 2D flatfield profile of one order/fiber
            for irow in range(len(current_mask)):
                if current_mask[irow] != 0:
                    xpixel.append(irow)
                    ridgeValues.append(masterFlat[irow, current_mask[irow]])
            xpixel = np.array(xpixel)
            ridgeValues = np.array(ridgeValues)

            coldefs = fits.ColDefs([fits.Column(name='x_pixel',  format='D', array=xpixel),             \
                                    fits.Column(name='spectrum', format='D', array=ridgeValues) ])
            header = Header({'EXTNAME': f"order: {order}, fiber: {fiber}"})
            new_table_hdu = fits.BinTableHDU.from_columns(columns=coldefs, header=header)
            hdu_list.append(new_table_hdu)

        hdu_list.writeto(output_path, overwrite=True)

        # That's it!

        return 




