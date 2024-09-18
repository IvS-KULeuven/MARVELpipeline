use std::fs;
use std::path::Path;
use std::time::Instant;
use fitsio::FitsFile;
use fitsio::images::{ImageType, ImageDescription};
use ndarray::{Array2, ArrayD, Axis};
use ndarray_ndimage::{gaussian_filter, BorderMode};
use configuration::parse_file;


type CCDImageType = ArrayD<i32>;




fn main() {

    let now = Instant::now();

    // Load and parse the param.yaml file to get the paths from which we
    // will load the files and to which we will save the output files.

    let config: serde_yaml::Value = parse_file();

    // Get all the paths from which we will read/write

    let order_tracing_config = &config["TwoDimensionalOrderMaskTracing"];
    let configurations = &config["Configuration"];

    let project_root = configurations.get("rootFolder").unwrap().as_str().unwrap();
    let master_flat_path = project_root.to_owned() + order_tracing_config.get("inputPath").unwrap().as_str().unwrap();
    let orders_mask_path = project_root.to_owned() + order_tracing_config.get("outputPath").unwrap().as_str().unwrap();
    let smoothed_master_flat_path = project_root.to_owned() + order_tracing_config.get("outputPathSmoothedMasterFlat").unwrap().as_str().unwrap();


    // Open the master flat file

    let mut fitsfile = FitsFile::open(master_flat_path).unwrap();
    let hdu = fitsfile.primary_hdu().unwrap();
    let master_flat: CCDImageType = hdu.read_image(&mut fitsfile).unwrap();

    let stdev_bias = hdu.read_key::<f32>(&mut fitsfile, "STD_BIAS").unwrap() as u32;
    let num_rows_ccd = hdu.read_key::<f32>(&mut fitsfile, "ROWS").unwrap() as usize;
    let num_cols_ccd = hdu.read_key::<f32>(&mut fitsfile, "COLS").unwrap() as usize;


    // Smooth the master flat


    let master_flat = master_flat.mapv(|elem| elem as f64);
    let stdev = 0.7;
    let order = 0;                                                                                             // convolve with Gaussian 
    let truncate = 3;
    let smoothed_master_flat = gaussian_filter(&master_flat, stdev, order, BorderMode::Reflect, truncate);
    let smoothed_master_flat = smoothed_master_flat.mapv(|elem| elem as u32);

    if Path::new(&smoothed_master_flat_path).exists() {                              // Remove previous file from an older run
        fs::remove_file(&smoothed_master_flat_path).unwrap();
    }

    let description = ImageDescription {
        data_type: ImageType::UnsignedShort,
        dimensions: &[num_rows_ccd, num_cols_ccd],
    };

    let smoothed_master_flat_path = Path::new(&smoothed_master_flat_path);
    if smoothed_master_flat_path.exists() {
        fs::remove_file(&smoothed_master_flat_path).unwrap();
    } else if !smoothed_master_flat_path.parent().unwrap().is_dir() {
        fs::create_dir_all(smoothed_master_flat_path.parent().unwrap()).unwrap();
    }

    let mut fitsfile = FitsFile::create(smoothed_master_flat_path)
        .with_custom_primary(&description)
        .open()
        .unwrap();
    let hdu = fitsfile.primary_hdu().unwrap();
    hdu.write_image(&mut fitsfile, smoothed_master_flat.as_slice().unwrap()).unwrap();



    // Find the maxima of all orders in all fibers along the middle row

    let middle = num_rows_ccd/2;
    let middle_row = smoothed_master_flat.index_axis(Axis(0), middle);
    let mut icol_maxima: Vec::<usize> = Vec::new();

    for icol in 1..num_cols_ccd-1 {
        if (middle_row[icol] > 10 * stdev_bias) && (middle_row[icol] >= middle_row[icol-1]) && (middle_row[icol] > middle_row[icol+1]) { 
                icol_maxima.push(icol);
        }
    }


    // Creating the order masks


    // Loop over all orders/fibers

    let num_spectra = icol_maxima.len();                                             // total number of orders over all fibers
    let mut lower_boundaries = Array2::<u32>::zeros((num_spectra, num_rows_ccd));
    let mut upper_boundaries = Array2::<u32>::zeros((num_spectra, num_rows_ccd));
    let mut max_ridge = Array2::<u32>::zeros((num_spectra, num_rows_ccd));

    for ispectrum in 0..num_spectra {

        // Trace the current order along all rows, first the rows below the middle, then above the middle.

        let mut icol_max_of_current_row = 0;

        // First trace the rows below the middle row

        for irow in (0..=middle).rev() {

            // Find the maximum of the current order for the current row

            if irow == middle {
                icol_max_of_current_row = icol_maxima[ispectrum];
            } else {
                // The curvature of the order is so mild that we only need to look in the neighbouring
                // pixels of the maximum of the previous row

                let mut max = smoothed_master_flat[[irow, icol_max_of_current_row]];
                
                let icol_min_range: usize =
                    if icol_max_of_current_row != 0 {
                        icol_max_of_current_row-1} else {0};

                let icol_max_range: usize =
                    if icol_max_of_current_row != num_cols_ccd {
                        icol_max_of_current_row+1} else {num_cols_ccd};

                
                for icol in [icol_min_range, icol_max_range] {
                    
                    // The orders can hit the left edge of the CCD as well as the top/bottom edge, depending on the
                    // order, so we need to check if icol is not negative. Since icol has type uint, this means checking
                    // whether it didn't wrap around and became a huge value.

                    if (icol < num_cols_ccd) && (smoothed_master_flat[[irow, icol]] > max) {
                            max = smoothed_master_flat[[irow, icol]];
                            icol_max_of_current_row = icol;
                    }
                }
            }

            // If the maximum is high enough above the noise level search for the order boundaries
            // If not, then we stop the order tracing on this end.

            if smoothed_master_flat[[irow, icol_max_of_current_row]] >= 10 * stdev_bias {

                max_ridge[[ispectrum, irow]] = icol_max_of_current_row as u32;

                // Search for the left boundary of the order 

                let mut icol_left_boundary = icol_max_of_current_row;
                while (smoothed_master_flat[[irow, icol_left_boundary]] > 5 * stdev_bias) & (icol_left_boundary > 0) {
                    icol_left_boundary -= 1;
                }

                // Search for the right boundary of the order 

                let mut icol_right_boundary = icol_max_of_current_row;
                while (smoothed_master_flat[[irow, icol_right_boundary]] > 5 * stdev_bias) & (icol_right_boundary < num_cols_ccd-1) {
                    icol_right_boundary += 1;
                }

                // Keep the column boundaries for the current order for this particular row

                lower_boundaries[[ispectrum,irow]] = icol_left_boundary as u32;
                upper_boundaries[[ispectrum,irow]] = icol_right_boundary as u32;

            } else {
                // We stop the order tracing on this end of the order right here.
                break;
            }


        }

        // Then trace the rows beyond the middle row

        // We restart from the row next to the middle one again, so reset the column of maximum flux to the one of the middle row.

        icol_max_of_current_row = icol_maxima[ispectrum];

        for irow in (middle+1)..num_rows_ccd {

            // Find the maximum of the current order for the current row
            // The curvature of the order is so mild that we only need to look in the neighbouring
            // pixels of the maximum of the previous row

            let mut max = smoothed_master_flat[[irow, icol_max_of_current_row]];

            let icol_min_range: usize =
                if icol_max_of_current_row != 0 {
                    icol_max_of_current_row-1} else {0};

            let icol_max_range: usize =
                if icol_max_of_current_row != num_cols_ccd {
                    icol_max_of_current_row+1} else {num_cols_ccd};

            for icol in [icol_min_range, icol_max_range] {
                // The orders can hit the left edge of the CCD as well as the top/bottom edge, depending on the order.
                // So we need to check if icol is not negative. Since icol has type uint, this means checking
                // whether it didn't wrap around and became a huge value.
                if (icol < num_cols_ccd) && (smoothed_master_flat[[irow, icol]] > max) {
                        max = smoothed_master_flat[[irow, icol]];
                        icol_max_of_current_row = icol;
                }
            }

            // If the maximum is high enough above the noise level search for the order boundaries
            
            if smoothed_master_flat[[irow, icol_max_of_current_row]] >= 10 * stdev_bias {

                max_ridge[[ispectrum, irow]] = icol_max_of_current_row as u32;

                // Search for the left boundary of the order 

                let mut icol_left_boundary = icol_max_of_current_row;
                while (smoothed_master_flat[[irow, icol_left_boundary]] > 5 * stdev_bias) & (icol_left_boundary > 0) {
                    icol_left_boundary -= 1;
                }

                // Search for the right boundary of the order 

                let mut icol_right_boundary = icol_max_of_current_row;
                while (smoothed_master_flat[[irow, icol_right_boundary]] > 5 * stdev_bias) & (icol_right_boundary < num_cols_ccd-1) {
                    icol_right_boundary += 1;
                }

                // Keep the column boundaries for the current order for this particular row

                lower_boundaries[[ispectrum,irow]] = icol_left_boundary as u32;
                upper_boundaries[[ispectrum,irow]] = icol_right_boundary as u32;

            } else {
                // We stop the order tracing on this end of the order right here.
                break;
            }
        }
    }


    // Set up the output directory

    let orders_path = Path::new(&orders_mask_path);

    if orders_path.exists() {                                         // Remove previous file from an older run
        fs::remove_file(orders_path).unwrap();
    } else if !orders_path.parent().unwrap().is_dir() {                       // Create directory if needed.
        let _ = fs::create_dir_all(orders_path.parent().unwrap());
    }

    // Save the order mask info to a FITS file

    let description = ImageDescription {
        data_type: ImageType::UnsignedShort,
        dimensions: &[num_spectra, num_rows_ccd],
    };
    let mut fitsfile = FitsFile::create(orders_mask_path)
        .with_custom_primary(&description)
        .open()
        .unwrap();

    let hdu = fitsfile.primary_hdu().unwrap();
    hdu.write_image(&mut fitsfile, max_ridge.as_slice().unwrap()).unwrap();

    let hdu = fitsfile.create_image("EXTNAME".to_string(), &description).unwrap();
    hdu.write_image(&mut fitsfile, lower_boundaries.as_slice().unwrap()).unwrap();

    let hdu = fitsfile.create_image("EXTNAME2".to_string(), &description).unwrap();
    hdu.write_image(&mut fitsfile, upper_boundaries.as_slice().unwrap()).unwrap();

    println!("[{:.1?}]", now.elapsed());

}

