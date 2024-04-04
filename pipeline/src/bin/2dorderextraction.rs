use std::path::Path;
use std::fs;
use std::time::Instant;
use fitsio::FitsFile;
use fitsio::images::{ImageType, ImageDescription};
use ndarray::{Array2, ArrayD, Axis, s};
use itertools::izip;

use configuration::parse_file;


type CCDImageType = ArrayD<u32>;



fn main() {

    let now = Instant::now();

    // Load and parse the param.yaml file to get the paths from which we
    // will load the files and to which we will save the output files.

    let config: serde_yaml::Value = parse_file();

    // Get all the paths from which we will read/write

    let mask_image_paths = &config["orderImage"];
    let science_image_paths = &config["rawScienceImage"];
    let configurations = &config["configuration"];

    let project_root = configurations.get("rootFolder").unwrap().as_str().unwrap();
    let order_mask_path = project_root.to_owned() + mask_image_paths.get("maskOutputpath").unwrap().as_str().unwrap();
    let c_science_paths = science_image_paths.get("outputpath").unwrap().as_sequence().unwrap();
    let extracted_orders_paths = mask_image_paths.get("scienceOutputpath").unwrap().as_sequence().unwrap();


    // Open the order mask file. This file contains three images, one with the index of the maximum position,
    //one with the lower bound index and one with the upper bound index.

    let mut fitsfile = FitsFile::open(order_mask_path).unwrap();
    let max_hdu = fitsfile.hdu(0).unwrap();
    let lower_hdu = fitsfile.hdu(1).unwrap();
    let upper_hdu = fitsfile.hdu(2).unwrap();

    let max_ridge: CCDImageType = max_hdu.read_image( &mut fitsfile).unwrap();
    let lower_boundaries: CCDImageType = lower_hdu.read_image(&mut fitsfile).unwrap();
    let upper_boundaries: CCDImageType = upper_hdu.read_image(&mut fitsfile).unwrap();



    // Loop over the science images

    for (c_science_path, e_science_path) in izip!(c_science_paths, extracted_orders_paths) {
        let science_path = project_root.to_owned() + c_science_path.as_str().unwrap();
        let output_path = project_root.to_owned() + e_science_path.as_str().unwrap();

        // Open the calibrated science images

        let mut fitsfile = FitsFile::open(science_path).unwrap();
        let hdu = fitsfile.primary_hdu().unwrap();
        let raw_science_image: CCDImageType = hdu.read_image(&mut fitsfile).unwrap();

    
        // extract its orders 

        let num_spectra: usize = max_ridge.len_of(Axis(0));
        let num_rows_ccd: usize = max_ridge.len_of(Axis(1));
        let mut image_science = Array2::<u32>::zeros((num_rows_ccd, num_rows_ccd));

        for ispectrum in 0..num_spectra {
            for irow in 0..num_rows_ccd {
                if (lower_boundaries[[ispectrum, irow]] == 0) & (upper_boundaries[[ispectrum, irow]] == 0) {
                    continue;
                } else {
                    let begin = lower_boundaries[[ispectrum, irow]] as isize;
                    let end   = upper_boundaries[[ispectrum, irow]] as isize + 1;
                    let slice = s![irow as isize, begin..end];
                
                    image_science.slice_mut(slice).assign(&raw_science_image.slice(slice));
                }
            }
        }


        // Paste the extracted orders on a 2D image and save it to a FITS file just for verification

        let path = Path::new(&output_path);
        if path.exists() {                                     // Remove previous filefrom an older run
            fs::remove_file(&output_path).unwrap();
        }  else if !path.parent().unwrap().is_dir() {         // Make sure parent directory exists
            fs::create_dir_all(path.parent().unwrap()).unwrap();}

        let description = ImageDescription {
            data_type: ImageType::UnsignedShort,
            dimensions: &[num_rows_ccd, num_rows_ccd],
        };
        let mut fitsfile = FitsFile::create(output_path)
            .with_custom_primary(&description)
            .open()
            .unwrap();
        let hdu = fitsfile.primary_hdu().unwrap();
        hdu.write_image(&mut fitsfile, image_science.as_slice().unwrap()).unwrap();
    }

    println!("[{:.1?}]", now.elapsed());

}
