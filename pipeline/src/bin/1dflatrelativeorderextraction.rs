use std::path::Path;
use std::time::Instant;
use fitsio::FitsFile;
use ndarray::{ArrayD, Array1, Axis, s};
use std::fs;
mod configuration;

type CCDImageType = ArrayD<u32>;








fn main() {
    
    let now = Instant::now();

    // Load and parse the param.yaml file to get the paths from which we
    // will load the files and to which we will save the output files.

    let config: serde_yaml::Value = configuration::configuration::parse_file();

    // Get all the paths from which we will read/write

    let mask_image_paths = &config["orderImage"];
    let science_image_paths = &config["rawScienceImage"];
    let configurations = &config["configuration"];
    let flat_images_path = &config["rawFlatImage"];
    let optimal_image_paths = &config["optimalOrderExtraction"];

    let project_root = configurations.get("rootFolder").unwrap().as_str().unwrap();
    let order_mask_path = project_root.to_owned() + mask_image_paths.get("maskOutpath").unwrap().as_str().unwrap();
    let master_flat_path= project_root.to_owned() + flat_images_path.get("outpath").unwrap().as_str().unwrap();
    let c_science_paths = science_image_paths.get("outpath").unwrap().as_sequence().unwrap();
    let output_directory = project_root.to_owned() +optimal_image_paths.get("outpath").unwrap().as_str().unwrap();

    // Open the order mask file. This file contains three images, one with the index of the maximum position,
    //one with the lower bound index and one with the upper bound index.

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

    
    for s_path in c_science_paths {

        let science_path = project_root.to_owned() + s_path.as_str().unwrap();
        let science = Path::new(&science_path);
        let file_stem = science.file_stem().unwrap().to_str().unwrap();
        let exposure_id = &file_stem[0..15];

        // Set up the output directory

        let output_directory = format!("{output_directory}/{exposure_id}");
        let path = Path::new(&output_directory);

        if !path.is_dir() {
            fs::create_dir_all(path).unwrap();}

        // Open the calibrated science file.

        let mut fitsfile = FitsFile::open(science).unwrap();
        let hdu = fitsfile.primary_hdu().unwrap();
        let science_image: CCDImageType = hdu.read_image(&mut fitsfile).unwrap();

        // Do the flat-relative optimal extraction using the method of Zechmeister et al. (2014)

        let mut mean_spectrum = Array1::<f64>::zeros(num_rows_ccd);  // For a given order, the mean spectrum averaged over all science fibers
        let mut num_fibers = Array1::<u8>::zeros(num_rows_ccd);        // Given an order, for each pixel, how many fibers contributed to the mean spectrum?

        
        for ispectrum in 0..num_spectra {
            let fiber = 5 - ispectrum % 5;                                 // Fiber 1 is the calibration fiber, 2,3,4,5 are the science fibers
            let order = 98 -  (ispectrum as f64 / 5.0).floor() as u32;      // Ranging from 98 (blue) down to 33 (red). Checked with Balmer lines.

            let mut spectrum: Vec::<f64> = Vec::with_capacity(num_rows_ccd);
            let mut xpixel: Vec::<u32> = Vec::with_capacity(num_rows_ccd); // Pixel coordinates on the x-axis on the optimal extracted spectrum

            for irow in 0..num_rows_ccd {
                if (lower_boundaries[[ispectrum, irow]] == 0) & (upper_boundaries[[ispectrum, irow]] == 0) {
                    continue;
                } else {
                    let begin = lower_boundaries[[ispectrum, irow]] as isize;
                    let end   = upper_boundaries[[ispectrum, irow]] as isize + 1;
                    let slice = s![irow as isize, begin..end];
                    let science_cross_order_slice = science_image.slice(slice).mapv(|x| x as f64);
                    let flat_cross_order_slice = master_flat.slice(slice);
                
                    // The following numbers came from Nich.
                
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

            // Save the orders to a csv file using the proper order and fiber number in the filename
            
            let output_path = format!("{output_directory}/{exposure_id}_order{order}_fiber{fiber}.csv");

            let mut outputfile = csv::Writer::from_path(output_path).unwrap();
            outputfile.write_record(["pixel","flat_relative_flux"]).unwrap();
            for (pix, spec) in xpixel.iter().zip(spectrum.iter()) { 
                outputfile.write_record(&[*pix as f64, *spec].map(|elem| elem.to_string()) ).unwrap();
            }
            outputfile.flush().unwrap();


            // Compute the mean spectrum per order over all 4 science fibers.
            // Note that the spectra of the same order but from different fibers may not have all the same length,
            // this can differ a few pixels, depending on the mask.
        

            if fiber == 1 {                                                  
                // Starting with another wavelength calibration fiber signals a new order.
                // First process the mean spectrum of the prevous order and save it.

                // Crop the mean spectrum to those pixels for which all 4 science fibers are available.

                let mut cropped_mean_spectrum: Vec::<f64> = Vec::with_capacity(num_rows_ccd);
                let mut cropped_xpixel: Vec::<f64> = Vec::with_capacity(num_rows_ccd);
                for irow in 0..num_rows_ccd {
                    if num_fibers[irow] == 4 {
                        cropped_mean_spectrum.push(mean_spectrum[irow] / 4.0);
                        cropped_xpixel.push(irow as f64);
                    }
                }

                // Save the mean spectrum of the previous order 

                let previous_order = order - 1;
                let output_path = format!("{output_directory}/{exposure_id}_order{previous_order}_mean.csv");

                let mut outputfile = csv::Writer::from_path(output_path).unwrap();
                outputfile.write_record(["pixel","flat_relative_flux"]).unwrap();
                for (pix, spec) in cropped_xpixel.iter().zip(cropped_mean_spectrum.iter()) { 
                    outputfile.write_record(&[*pix, *spec].map(|elem| elem.to_string()) ).unwrap();
                }
                outputfile.flush().unwrap();

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
