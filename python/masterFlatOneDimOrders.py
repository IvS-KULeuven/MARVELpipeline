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
        ridges = mask_fits_file[0].data                  # Column values of the ridge of a particular order
        lower_edges = mask_fits_file[1].data             # Lower column values of the edge of each order
        upper_edges = mask_fits_file[2].data             # Upper column vlaues of the edge of each order

        # Loop over all masks (one per order / fiber) and save the 1D flatfield profiles in a FITS file

        hdu_list = fits.HDUList([fits.PrimaryHDU()])       # Compulsory but empty primary HDU

        for imask in range(len(ridges)):

            fiber = 5 - imask % 5                          # Fiber 1 is the calibration fiber, 2,3,4,5 are the science fibers
            order = 98 -  int(np.floor(imask / 5.0))       # Ranging from 98 (blue) down to 33 (red). 

            current_ridge = ridges[imask]                  # Column value of the maximum of the cross-section of current order
            current_lower_edge = lower_edges[imask]        # Lower column values of the cross-section of current order
            current_upper_edge = upper_edges[imask]        # Upper column values of the cross-section of current order
            xpixel = []                                    # Along order pixel coordinate
            max_values = []                                # Max value of the order cross section at each row (the order ridge)
            sum_values = []                                # Mean value of the order cross section at each row 
            for irow in range(len(current_ridge)):
                if current_ridge[irow] != 0:
                    xpixel.append(irow)
                    max_values.append(masterFlat[irow, current_ridge[irow]])
                    sum = masterFlat[irow, current_lower_edge[irow]:(current_upper_edge[irow]+1)].sum()
                    sum_values.append(sum)

            xpixel = np.array(xpixel)
            max_values = np.array(max_values)
            sum_values = np.array(sum_values)

            coldefs = fits.ColDefs([fits.Column(name='x_pixel',  format='D', array=xpixel),             \
                                    fits.Column(name='ridge', format='D', array=max_values),            \
                                    fits.Column(name='integrated', format='D', array=sum_values) ])

            header = Header({'EXTNAME': f"order: {order}, fiber: {fiber}"})
            new_table_hdu = fits.BinTableHDU.from_columns(columns=coldefs, header=header)
            hdu_list.append(new_table_hdu)

        hdu_list.writeto(output_path, overwrite=True)

        # That's it!

        return 




