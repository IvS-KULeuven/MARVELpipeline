use std::path::Path;
use std::fs;
use std::time::Instant;
use fitsio::FitsFile;
use fitsio::tables::{ColumnDescription, ColumnDataType};
use ndarray::{ArrayD, Array1, Axis, s};
use itertools::izip;
use configuration::parse_file;

type CCDImageType = ArrayD<u32>;



fn main() {
    
    let now = Instant::now();

    // Load and parse the param.yaml file to get the paths from which we
    // will load the files and to which we will save the output files.

    let config: serde_yaml::Value = parse_file();

    // Get all the paths from which we will read/write

    let configurations = &config["Configuration"];
    let two_dim_order_mask_config = &config["TwoDimensionalOrderMaskTracing"];
    let optimal_order_extraction_config = &config["OptimalOrderExtraction"];

    let project_root = configurations.get("rootFolder").unwrap().as_str().unwrap();
    let master_flat_path= project_root.to_owned() + two_dim_order_mask_config.get("inputPath").unwrap().as_str().unwrap();
    let order_mask_path = project_root.to_owned() + two_dim_order_mask_config.get("outputPath").unwrap().as_str().unwrap();
    let science_paths = optimal_order_extraction_config.get("inputPath").unwrap().as_sequence().unwrap();    
    let output_directories = optimal_order_extraction_config.get("outputPath").unwrap().as_sequence().unwrap();

    // Open the order mask file. This file contains three images, one with the index of the maximum position,
    // one with the begin column index and one with the end column index of each order.

    let mut fitsfile = FitsFile::open(order_mask_path).unwrap();
    let max_hdu = fitsfile.hdu(0).unwrap();
    let lower_hdu = fitsfile.hdu(1).unwrap();
    let upper_hdu = fitsfile.hdu(2).unwrap();

    let max_ridge: CCDImageType = max_hdu.read_image( &mut fitsfile).unwrap();
    let lower_boundaries: CCDImageType = lower_hdu.read_image(&mut fitsfile).unwrap();
    let upper_boundaries: CCDImageType = upper_hdu.read_image(&mut fitsfile).unwrap();

    let num_spectra: usize = max_ridge.len_of(Axis(0));
    let num_rows_ccd: usize = max_ridge.len_of(Axis(1));

    // Open the master flat file.

    let mut fitsfile = FitsFile::open(master_flat_path).unwrap();
    let hdu = fitsfile.primary_hdu().unwrap();
    let master_flat: ArrayD<f64>  = hdu.read_image(&mut fitsfile).unwrap();
    let stdev_bias = hdu.read_key::<f32>(&mut fitsfile, "STD_BIAS").unwrap() as u32;


    for (science_path, output_path) in izip!(science_paths, output_directories) {
        
        let science_path = project_root.to_owned() + science_path.as_str().unwrap();
        let science_path = Path::new(&science_path);

        // Open the bias subtracted science image.

        let mut fitsfile = FitsFile::open(science_path).unwrap();
        let mut hdu = fitsfile.primary_hdu().unwrap();
        let science_image: CCDImageType = hdu.read_image(&mut fitsfile).unwrap();

        // Do the flat-relative optimal extraction using the method of Zechmeister et al. (2014)

        let mut mean_spectrum = Array1::<f64>::zeros(num_rows_ccd);    // For a given order, the mean spectrum averaged over all science fibers
        let mut num_fibers = Array1::<u8>::zeros(num_rows_ccd);        // Given an order, for each pixel, how many fibers contributed to the mean spectrum?


        // Set up the output directory

        let output_path = project_root.to_owned() + output_path.as_str().unwrap();
        let output_path = Path::new(&output_path);

        if !output_path.parent().unwrap().is_dir(){                // Create directory if needed.
            fs::create_dir_all(output_path.parent().unwrap()).unwrap();}

        if output_path.exists() {                                  // Remove previous file form an older run
            fs::remove_file(output_path).unwrap();}

        let mut fitsfile  = FitsFile::create(output_path)
            .open()
            .unwrap();

        // We set up the description of the table that will store all the data for one order

        let fiber_x_description = ColumnDescription::new("x_pixel")
            .with_type(ColumnDataType::Int)
            .create()
            .unwrap();

        let fiber_science_description = ColumnDescription::new("spectrum")
            .with_type(ColumnDataType::Float)
            .create()
            .unwrap();


        let descriptions = [fiber_x_description, fiber_science_description];
      
        for ispectrum in 0..num_spectra {
            let fiber = 5 - ispectrum % 5;                                 // Fiber 1 is the calibration fiber, 2,3,4,5 are the science fibers
            let order = 98 -  (ispectrum as f64 / 5.0).floor() as u32;     // Ranging from 98 (blue) down to 33 (red). Checked with Balmer lines.

            // When we access fiber 5 we start a new order with new table

            hdu = fitsfile.create_table(format!("order: {}, fiber: {}", order, fiber).to_string(), &descriptions).unwrap();


            let mut spectrum: Vec::<f64> = Vec::new();
            let mut xpixel: Vec::<u32> = Vec::new(); // Pixel coordinates on the x-axis on the optimal extracted spectrum

            let mut i=0;
            for irow in 0..num_rows_ccd {
                if (lower_boundaries[[ispectrum, irow]] == 0) & (upper_boundaries[[ispectrum, irow]] == 0) {
                    continue;
                } else {
                    i = i+1;
                    let begin = lower_boundaries[[ispectrum, irow]] as isize;
                    let end   = upper_boundaries[[ispectrum, irow]] as isize + 1;
                    
                    let slice = s![irow as isize, begin..end];
                    let science_cross_order_slice = science_image.slice(slice).mapv(|x| x as f64);
                    let flat_cross_order_slice = master_flat.slice(slice);
                
                    // The following numbers came from Nich
                
                    let readout_noise = 5.0;                                   // [e-/pix]
                    let gain = 3.0;                                            // [e-/ADU] 
                    let readout_noise = readout_noise / gain;                  // [ADU/pix]

                    // Note: if a skybackground would have been subtracted, we would also need to take into account its photon noise
                    //       in the weights.
                    // Fiber 5 is the etalon emission spectrum. No 4σ thresholding there, because then we wouldn't have the close-to-zero flux 
                    // in between the etalon lines.

                    let weights = if fiber == 1 {
                        science_cross_order_slice.mapv(|x| 1.0 / (readout_noise*readout_noise + x)) 
                    } else {
                        science_cross_order_slice.mapv(|x| 
                            if x > (4*stdev_bias) as f64 { 
                                1.0 / (readout_noise*readout_noise + x)
                            } else { 
                                0.0 
                            })
                    };

                    // The following is a bit of a quirk of ndarray:
                    //        owned * view
                    // or     &view * &view
                    // seem to work but not
                    //        view * owned
                    // or     view * view

                    spectrum.push( (science_cross_order_slice * flat_cross_order_slice * &weights).sum()
                        / (&flat_cross_order_slice * &flat_cross_order_slice * &weights).sum() );
                    xpixel.push(irow as u32);
                }

            }
                
            hdu.write_col(&mut fitsfile, format!("x_pixel"),  &xpixel).unwrap();
            hdu.write_col(&mut fitsfile, format!("spectrum"), &spectrum).unwrap();

            // Compute the mean spectrum per order over all 4 science fibers.
            // Note that the spectra of the same order but from different fibers may not have all the same length,
            // this can differ a few pixels, depending on the mask.
        
            if fiber == 1 {                                                  
                // Starting with another wavelength calibration fiber signals a new order.
                // First process the mean spectrum of the prevous order and save it.

                // Crop the mean spectrum to those pixels for which all 4 science fibers are available.

                let fiber_x_description = ColumnDescription::new("x_pixel")
                    .with_type(ColumnDataType::Int)
                    .create()
                    .unwrap();

                let fiber_science_description = ColumnDescription::new("spectrum")
                    .with_type(ColumnDataType::Float)
                    .create()
                    .unwrap();

                let descriptions = [fiber_x_description, fiber_science_description];

                hdu = fitsfile.create_table(format!("order: {}, mean", order).to_string(), &descriptions).unwrap();

                let mut cropped_mean_spectrum: Vec::<f64> = Vec::new();
                let mut cropped_xpixel: Vec::<f64> = Vec::new();
                
                for irow in 0..num_rows_ccd {
                    if num_fibers[irow] == 4 {
                        cropped_mean_spectrum.push(mean_spectrum[irow] / 4.0);
                        cropped_xpixel.push(irow as f64);
                    }
                }

                // Add the mean spectrum to the output fits file

                hdu.write_col(&mut fitsfile, "x_pixel", &cropped_xpixel).unwrap();
                hdu.write_col( &mut fitsfile, "spectrum", &cropped_mean_spectrum).unwrap();

                // Reset the arrays to accommodate the current order. Fiber == 1 is the wavelength calibration order, 
                // so we can skip this one.

                mean_spectrum.iter_mut().for_each(|x| *x = 0.0);
                num_fibers.iter_mut().for_each(|x| *x = 0);

            } else {

                // We're dealing with a science order: contribute to the mean

                for (n, irow) in xpixel.iter().enumerate() {
                    mean_spectrum[*irow as usize] += spectrum[n];
                    num_fibers[*irow as usize] += 1;
                }
            }
        }
        
        

    }
    println!("[{:.1?}]", now.elapsed());
}
