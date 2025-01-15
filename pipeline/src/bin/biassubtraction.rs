use std::fs;
use std::path::Path;
use std::time::Instant;
use fitsio::FitsFile;
use fitsio::images::{ImageType, ImageDescription};
use ndarray::{ArrayD, Axis};
use itertools::izip;
use clap::Parser;

use configuration::parse_file;

type CCDImageType = ArrayD<u32>;


// Create a struct for the command line arguments

#[derive(Parser, Debug)]
struct Args {
    #[arg(short, long)]
    config_path: String,
}



fn main() {

    let now = Instant::now();

    // Load and parse the 'param.yaml' file to get the paths from which we
    // will load the files and to which we will save the output files.

    let args = Args::parse();
    let config: serde_yaml::Value = parse_file(&args.config_path);

    // Get all the paths rom which we will read/write

    let science_image_paths = &config["BiasSubtractedScienceImage"];
    let bias_image_paths = &config["MasterBiasImage"];
    let configurations = &config["Configuration"];

    let project_root = configurations.get("rootFolder").unwrap().as_str().unwrap();
    let master_bias_path = project_root.to_owned() + bias_image_paths.get("outputPath").unwrap().as_str().unwrap();
    let raw_science_paths = science_image_paths.get("inputPath").unwrap().as_sequence().unwrap();
    let cal_science_paths = science_image_paths.get("outputPath").unwrap().as_sequence().unwrap();
    

    // Open the master bias file

    let mut fitsfile = FitsFile::open(master_bias_path).unwrap();
    let hdu = fitsfile.primary_hdu().unwrap();
    let master_bias: ArrayD<i32> = hdu.read_image(&mut fitsfile).unwrap();


    // Substract the master bias from every science file

    for (tail_r_path, tail_c_path)  in izip!(raw_science_paths, cal_science_paths) {
        let science_path = project_root.to_owned() + tail_r_path.as_str().unwrap();
        let output_path = project_root.to_owned() + tail_c_path.as_str().unwrap();

        // open the science file

        let mut fitsfile = FitsFile::open(science_path).unwrap();
        let hdu = fitsfile.primary_hdu().unwrap();

        let science_image: ArrayD<i32> = hdu.read_image(&mut fitsfile).unwrap();

        // substract the master bias        

        let bias_substracted_science: ArrayD<i32>  = science_image - master_bias.clone();

        // if image has a negative value, then we add an offset

        let bias_substracted_science: CCDImageType = bias_substracted_science.mapv(
            |elem| { if elem < 0 {0_32} else {elem as u32}});

        // save image to fits file

        let path = Path::new(&output_path);
        if path.exists() {                                  // Remove previous file form an older run
            fs::remove_file(path).unwrap();
        } else if !path.parent().unwrap().is_dir() {        // Make sure parent directory exists
            fs::create_dir(path.parent().unwrap()).unwrap();
        }        


        let num_cols_ccd = bias_substracted_science.len_of(Axis(0));
        let num_rows_ccd =  bias_substracted_science.len_of(Axis(1));
            
        let description = ImageDescription {
            data_type: ImageType::UnsignedShort,
            dimensions: &[num_rows_ccd, num_cols_ccd],
        };
   
        let mut fitsfile = FitsFile::create(output_path)
            .with_custom_primary(&description)
            .open()
            .unwrap();
        let hdu = fitsfile.primary_hdu().unwrap();
        hdu.write_image(&mut fitsfile, bias_substracted_science.as_slice().unwrap()).unwrap();

    }
        
    println!("[{:.1?}]", now.elapsed());
}
