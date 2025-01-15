use std::fs;
use std::path::Path;
use itertools::izip;
use std::time::Instant;
use fitsio::FitsFile;
use fitsio::tables::{ColumnDescription, ColumnDataType};
use nalgebra::DVector;
use varpro::prelude::*;
use varpro::solvers::levmar::{LevMarProblemBuilder, LevMarSolver};
use clap::Parser;
use configuration::parse_file;


fn gaussian(x: &DVector<f64>, mu: f64, sigmasq: f64) -> DVector<f64> {
    x.map(|z| (-0.5 * ((z-mu).powi(2)/sigmasq)).exp())
}


fn gaussian_dmu(x: &DVector<f64>, mu: f64, sigmasq: f64) -> DVector<f64> {
    x.map(|z| (z-mu)/sigmasq * (-0.5 * ((z-mu).powi(2)/sigmasq)).exp())
}


fn gaussian_dsigmasq(x: &DVector<f64>, mu: f64, sigmasq: f64) -> DVector<f64> {
    x.map(|z| 0.5 * (z-mu).powi(2)/sigmasq.powi(2) * (-0.5 * ((z-mu).powi(2)/sigmasq)).exp())
}



pub fn max(data: &[f64]) -> f64 {
    let init_max: &f64 = &data[0];
    data.iter().fold(*init_max, |maxvalue, &value| { if maxvalue > value { maxvalue } else { value } })
}


pub fn min(data: &[f64]) -> f64 {
    let init_min: &f64 = &data[0];
    data.iter().fold(*init_min, |minvalue, &value| { if minvalue < value { minvalue } else { value } })
}


pub fn mean(y: &[f64]) -> f64 {
    y.iter().sum::<f64>() / (y.len() as f64)                          // weights must be normalized
}
 

pub fn stdev(y: &[f64], mean: f64) -> f64 {
    (y.iter().map(|x| (x-mean).powi(2)).sum::<f64>() / (y.len() as f64 - 1.0)).sqrt()
}




// Create a struct for the command line arguments

#[derive(Parser, Debug)]
struct Args {
    #[arg(short, long)]
    config_path: String,
}





fn main() {

    // Starting timing this pipeline module so that we can report how long it takes

    let now = Instant::now();

    // Load and parse the param.yaml file to get the paths of the FITS files we need to process. 

    let args = Args::parse();
    let config: serde_yaml::Value = parse_file(&args.config_path);

    // Get all the paths from which we will read/write

    let general_config = &config["Configuration"];
    let root_folder_processed_data = general_config.get("rootFolderProcessedData").unwrap().as_str().unwrap();

    let etalon_config = &config["EtalonPeakFitting"];
    let input_paths = etalon_config.get("inputPath").unwrap().as_sequence().unwrap();
    let output_paths = etalon_config.get("outputPath").unwrap().as_sequence().unwrap();

    // Loop over all science frames, the output of the peak fitting will be written to a separate FITS file.

    for (input_path, output_path) in izip!(input_paths, output_paths) {

        // Construct the absolute path of the input FITS file containing the 1D optimal extracted orders

        let input_path = root_folder_processed_data.to_owned() + input_path.as_str().unwrap();

        // Construct the absolute path of the output FITS file that will contain the etalon peak fit parameters 
        // If the output file already exists, delete the file.

        let output_path = root_folder_processed_data.to_owned() + output_path.as_str().unwrap();
        let output_path = Path::new(&output_path);
        if output_path.exists() {                                           // Remove previous file from an older run
            fs::remove_file(output_path).unwrap();
        }

        // Open the output FITS file

        let mut output_fits_file  = FitsFile::create(output_path).open().unwrap();

        // Set the descriptions of each of the columns in the output FITS table

        let offset_description = ColumnDescription::new("offset")
            .with_type(ColumnDataType::Float)
            .create()
            .unwrap();

        let mu_description = ColumnDescription::new("mu")
            .with_type(ColumnDataType::Float)
            .create()
            .unwrap();

        let sigma_description = ColumnDescription::new("sigma")
            .with_type(ColumnDataType::Float)
            .create()
            .unwrap();

        let height_description = ColumnDescription::new("height")
            .with_type(ColumnDataType::Float)
            .create()
            .unwrap();

        let offset_err_description = ColumnDescription::new("offset_err")
            .with_type(ColumnDataType::Float)
            .create()
            .unwrap();

        let mu_err_description = ColumnDescription::new("mu_err")
            .with_type(ColumnDataType::Float)
            .create()
            .unwrap();

        let sigma_err_description = ColumnDescription::new("sigma_err")
            .with_type(ColumnDataType::Float)
            .create()
            .unwrap();

        let height_err_description = ColumnDescription::new("height_err")
            .with_type(ColumnDataType::Float)
            .create()
            .unwrap();

        let descriptions = [offset_description, mu_description, sigma_description, height_description,
                            offset_err_description, mu_err_description, sigma_err_description, height_err_description];

        // Open the input fits file. Each order has a separate HDU. Loop over all orders (HDUs), reading
        // only the 1st fiber which is the etalon fiber. Skip the first (# 0) HDU because this is the primary
        // HDU which does not contain a spectrum.

        let mut input_fits_file = FitsFile::open(input_path).unwrap();
        let num_hdus = input_fits_file.iter().count();

        for ihdu in 1..num_hdus {

            let hdu = input_fits_file.hdu(ihdu).unwrap();

            // Read the order information. Example: the key 'EXTNAME' could contain the string "order: 98, fiber: 2".

            let order_fiber_id: String = hdu.read_key(&mut input_fits_file, "EXTNAME").unwrap();

            // Split this string into the order and the fiber part

            let order_fiber_id: Vec<&str> = order_fiber_id.split(',').collect();

            // If this does not concern fiber 1, then simply skip this HDU

            if order_fiber_id[1] != " fiber: 1" {
                continue;
            }

            // Read the relevant columns. Convert the pixel coordinates from integers to floats,
            // as they will be used in a fitting algorithm.

            let xpixel: Vec<i32> = hdu.read_col(&mut input_fits_file, "x_pixel").unwrap();
            let xpixel: Vec<f64> = xpixel.iter().map(|&value| value as f64).collect();
            let flux:   Vec<f64> = hdu.read_col(&mut input_fits_file, "spectrum").unwrap();

            // Find all the peaks of the etalon emission lines. We need to take care not to erronously identify 
            // an emission line as two lines in the rare cases that:
            //     a) the top of the line is a plateau of 2 points 
            //     b) the top of the line are 2 neighboring peaks (3 points: up-down-up) 
            // In the case of a plateau we select the left most point.
            // We start from n=6 because we're not interested in half peaks at the edges of the etalon spectrum. 
            // The thresholding with respect to the minimum flux in the local neighborhood is to avoid the small 
            // peaks in the continuum in between the emission lines. We can't take a global minimum because the 
            // etalon spectrum may not be perfectly normalized and can be curved upwards at the edges.

            let mut ipeaks: Vec::<usize> = Vec::new();
            for n in 6..xpixel.len()-6 {
                if (flux[n-2] < flux[n]) && (flux[n-1] < flux[n]) && (flux[n] >= flux[n+1]) && (flux[n] >= flux[n+2])
                   && (flux[n] > 0.1) {
                    ipeaks.push(n);
                }
            } 

            // The following arrays will be written in a FITS table 

            let mut offset:     Vec<f64> = vec![];
            let mut mu:         Vec<f64> = vec![]; 
            let mut sigma:      Vec<f64> = vec![];
            let mut height:     Vec<f64> = vec![];
            let mut offset_err: Vec<f64> = vec![];
            let mut mu_err:     Vec<f64> = vec![]; 
            let mut sigma_err:  Vec<f64> = vec![];
            let mut height_err: Vec<f64> = vec![];

            // Fit each of the peaks with a Gaussian profile

            for i in 0..ipeaks.len() {

                // The fitting module expects nalgebra matrices rather than built-in vectors, so convert.
                // Only select the part relevant for the current emission line to fit.

                let range_of_one_peak = ipeaks[i]-5..=ipeaks[i]+5;
                let x_onepeak: DVector<f64> = DVector::from_vec(xpixel[range_of_one_peak.clone()].to_vec());
                let flux_onepeak: DVector<f64> = DVector::from_vec(flux[range_of_one_peak].to_vec());

                // Fit a constant background + a single Gaussian to the emission line
                // The constant background is necessary as the flat-relative etalon spectrum may be 
                // curved upwards at the edges.

                let error_message = format!("Problem fitting peak # {} in {}", i, order_fiber_id[0]);

                let model = SeparableModelBuilder::<f64>::new(["mu", "sigmasq"])
                    .invariant_function(|x| DVector::from_element(x.nrows(), 1.0))
                    .function(&["mu", "sigmasq"], gaussian)
                    .partial_deriv("mu", gaussian_dmu)
                    .partial_deriv("sigmasq", gaussian_dsigmasq)
                    .independent_variable(x_onepeak)
                    .initial_parameters(vec![xpixel[ipeaks[i]], 1.0])
                    .build()
                    .expect(&error_message);

                let problem = LevMarProblemBuilder::new(model)
                    .observations(flux_onepeak)
                    .build()
                    .expect(&error_message);

                let (fit_result, fit_statistics) = LevMarSolver::new()
                    .fit_with_statistics(problem)
                    .expect("Problem with fitting the etalon lines");

                let nonlinear_params = fit_result.nonlinear_parameters();
                let linear_params    = fit_result.linear_coefficients().unwrap();
                let variance_nonlinear_params = fit_statistics.nonlinear_parameters_variance();
                let variance_linear_params = fit_statistics.linear_coefficients_variance();

                // Verify if the gaussian fit led to a sane result. If so, save the fit parameters.

                let fitted_offset     = linear_params[0];
                let fitted_mu         = nonlinear_params[0];
                let fitted_sigma      = nonlinear_params[1].sqrt();
                let fitted_height     = linear_params[1];
                let fitted_offset_err = variance_linear_params[0].sqrt();
                let fitted_mu_err     = variance_nonlinear_params[0].sqrt();
                let fitted_sigma_err  = variance_nonlinear_params[1].sqrt();
                let fitted_height_err = variance_linear_params[1].sqrt();

                // Only include the line in the list, if there are no counter indications that the fit succeeded. 

                if (fitted_sigma < 3.0) && (fitted_height > 0.0) && (fitted_mu_err != 0.0) && (fitted_sigma_err != 0.0)
                    && (fitted_mu_err < 0.2 * fitted_mu) && (fitted_sigma_err < 0.2 * fitted_sigma) {
                    offset.push(fitted_offset);
                    mu.push(fitted_mu);
                    sigma.push(fitted_sigma);
                    height.push(fitted_height);

                    offset_err.push(fitted_offset_err);
                    mu_err.push(fitted_mu_err);
                    sigma_err.push(fitted_sigma_err);
                    height_err.push(fitted_height_err);
                }   
            }

            let fits_table = output_fits_file.create_table(order_fiber_id[0], &descriptions).unwrap();
            fits_table.write_col(&mut output_fits_file, "offset", &offset).unwrap();
            fits_table.write_col(&mut output_fits_file, "mu", &mu).unwrap();
            fits_table.write_col(&mut output_fits_file, "sigma", &sigma).unwrap();
            fits_table.write_col(&mut output_fits_file, "height", &height).unwrap();
            fits_table.write_col(&mut output_fits_file, "offset_err", &offset_err).unwrap();
            fits_table.write_col(&mut output_fits_file, "mu_err", &mu_err).unwrap();
            fits_table.write_col(&mut output_fits_file, "sigma_err", &sigma_err).unwrap();
            fits_table.write_col(&mut output_fits_file, "height_err", &height_err).unwrap();

        }
    }

    // Print out how long this pipeline module took to execute

    println!("[{:.1?}]", now.elapsed());
}


