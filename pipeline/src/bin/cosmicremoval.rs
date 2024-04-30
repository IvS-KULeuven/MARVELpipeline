use std::{fs, vec};
use fitsio::FitsFile;
use itertools::Itertools;
use std::time::Instant;

use std::path::Path;
use ndarray::{Array, Array1, Array2, ArrayD, Axis, s, concatenate, array};
use ndarray_stats::QuantileExt;

use fitsio::images::{ImageType, ImageDescription};







type CCDImageType = ArrayD<u32>;

fn main() {

    let now = Instant::now();

    // We get the paths to all the nessecairy files

    let mask_path = "/home/driess/Projects/MARVELpipeline/Data/ProcessedData/ExtractedOrders/Mask/20230504_2d_mask.fits";
    let science_path = "/home/driess/Projects/MARVELpipeline/Data/ProcessedData/BiasSubtractedScience/20230614T185728_ESSSS_H_0600_bias_subtracted.fits";
    let flat_path = "/home/driess/Projects/MARVELpipeline/Data/ProcessedData/MasterFlat/20230504_MasterFlat.fits";

    // Open up extracted order mask file

    let mut fitsfile = FitsFile::open(mask_path).unwrap();
    let max_hdu = fitsfile.hdu(0).unwrap();
    let lower_hdu = fitsfile.hdu(1).unwrap();
    let upper_hdu = fitsfile.hdu(2).unwrap();

    let max_ridge: CCDImageType = max_hdu.read_image(&mut fitsfile).unwrap();
    let lower_boundaries: CCDImageType = lower_hdu.read_image(&mut fitsfile).unwrap();
    let upper_boundaries: CCDImageType = upper_hdu.read_image(&mut fitsfile).unwrap();

    let num_spectra: usize = max_ridge.len_of(Axis(0));
    let num_rows_ccd: usize = max_ridge.len_of(Axis(1));

    // Open up the science frame

    let mut fitsfile = FitsFile::open(science_path).unwrap();
    let hdu = fitsfile.primary_hdu().unwrap();

    let science_image: CCDImageType = hdu.read_image(&mut fitsfile).unwrap();
    let num_columns_ccd: usize = science_image.len_of(Axis(0));

    // Open up the master flat frame

    let mut fitsfile = FitsFile::open(flat_path).unwrap();
    let hdu = fitsfile.primary_hdu().unwrap();

    let flat_image: CCDImageType = hdu.read_image( &mut fitsfile).unwrap();

    // Set up the parameters to get the noise

    let readout_noise = 5.0;                                   // [e-/pix]
    let gain = 3.0;                                            // [e-/ADU] 
    let readout_noise = readout_noise / gain;                  // [ADU/pix]

    // Set up the paths to the output files.

    let unfiltered_difference_path = Path::new("/home/driess/Projects/MARVELpipeline/Data/ProcessedData/CosmicMask/unfiltered_difference.fits");
    let median_filtered_difference_path = Path::new("/home/driess/Projects/MARVELpipeline/Data/ProcessedData/CosmicMask/median_filtered_difference.fits");    
    let mask_path = Path::new("/home/driess/Projects/MARVELpipeline/Data/ProcessedData/CosmicMask/mask.fits");


    if !mask_path.parent().unwrap().is_dir() {
        fs::create_dir_all(mask_path.parent().unwrap()).unwrap();}

    // Initialize the image with the difference between the image and scaled flat image.

    let mut difference_expected_actual = Array2::<f64>::zeros((num_rows_ccd, num_columns_ccd));
    let mut mask = Array2::<usize>::ones((num_rows_ccd, num_columns_ccd));
    let mut profile = Array2::<f64>::zeros((num_rows_ccd, num_columns_ccd));


    // Create a image that reflext the difference between the observed spectra and the expected spectra from the flat image.

    for ispectrum in 0..num_spectra {

        // We neglect the calibration fiber 

        if ispectrum % 5 == 4 {
            continue;
        }
        
        for irow in 0..num_rows_ccd {
            
            let begin = (lower_boundaries[[ispectrum, irow]]) as isize;
            let end = (upper_boundaries[[ispectrum, irow]]) as isize + 1;
            let slice = s![irow as isize, begin..end];

            // We skip rows that don't extract pixels

            if (end - begin) == 1 {
                continue;
            }

            // We get the flat and science spectrum in the cross-order direction and
            // the sum of the flux in that direction.

            let flat_cross_order_slice = flat_image
                .slice(slice)
                .mapv(|x| x as f64);
            let box_extracted_flat_stripe = flat_cross_order_slice.sum();

            let science_cross_order_slice = science_image
                .slice(slice)
                .mapv(|x| x as f64);
            let stand_spec = science_cross_order_slice.sum();

            // The profile describes the shape of the flat image in the cross-order direction.

            profile.slice_mut(slice).assign(&(flat_cross_order_slice / box_extracted_flat_stripe));

            // Calculate the first guess of the spectrum (expected) and compare this
            // to the observed spectrum.
            
            let expected = profile.slice(slice).to_owned() * stand_spec;
            let actual = &science_cross_order_slice.clone() * &(mask.slice(slice).mapv(|x| x as f64));

            difference_expected_actual.slice_mut(slice).assign(&(actual-expected));
        }
    }


    
    if unfiltered_difference_path.exists() {
        fs::remove_file(unfiltered_difference_path).unwrap();
    }
    let description = ImageDescription {
        data_type: ImageType::Float,
        dimensions: &[num_rows_ccd, num_columns_ccd]
    };
    let mut fitsfile = FitsFile::create(unfiltered_difference_path)
        .with_custom_primary(&description)
        .open()
        .unwrap();
    let hdu = fitsfile.primary_hdu().unwrap();
    hdu.write_image( &mut fitsfile, difference_expected_actual.clone().as_slice().unwrap()).unwrap();



    
    // We apply a median filter on the difference between the expected and measured spectrum. 
    
    let median_filtered_difference = median_filter(&difference_expected_actual, 201).mapv(|x| x.abs());


    if median_filtered_difference_path.exists() {
        fs::remove_file(median_filtered_difference_path).unwrap();
    }

    let mut fitsfile = FitsFile::create(    median_filtered_difference_path)
        .with_custom_primary(&description)
        .open()
        .unwrap();
    let hdu = fitsfile.primary_hdu().unwrap();
    hdu.write_image( &mut fitsfile, median_filtered_difference.clone().as_slice().unwrap()).unwrap();



    // We keep track of the lines that still have pixels that we want to exclude

    let mut lines_to_consider = (0..num_spectra).collect_vec();
    while lines_to_consider.len() > 0 {

        let mut next_lines_to_consider : Vec<usize> = Vec::new();

        for ispectrum in lines_to_consider.clone().into_iter() {

            // We neglect the calibration fiber 

            if ispectrum % 5 == 4 {
                continue;
            }

            // Initialize the values and indices of those lines we can exclude.

            let mut reject_value : Array1<f64> = Array1::<f64>::zeros(num_rows_ccd);
            let mut reject_index : Array1<isize> = Array1::<isize>::zeros( num_rows_ccd);


            for irow in 0..num_rows_ccd {
                
                let begin = (lower_boundaries[[ispectrum, irow]]) as isize;
                let end = (upper_boundaries[[ispectrum, irow]]) as isize + 1;
                let slice = s![irow as isize, begin..end];                       

                if (end - begin) == 1 {
                    continue;
                }
                    

                let science_cross_order_slice = science_image.slice(slice)
                    .mapv(|x| x as f64);
                let profile_row = profile.slice(slice).to_owned();
                let stand_spec = science_cross_order_slice.sum();

                // Calculate the first guess of the difference between stripe and scaled flat stripe
                let mask_line = mask.slice(slice)
                    .to_owned().mapv(|x| x as f64);
                let expected = &profile_row * &(stand_spec * &mask_line) ;
                let data_var = (stand_spec * profile_row).mapv(|x| x.abs() / gain + readout_noise) ;
                let noise_rev = data_var.mapv(|x| 1./ x.sqrt() );


                let actual = &mask_line * &science_cross_order_slice.clone();
                let diff = actual - expected;
            
                let reject = diff.mapv(|x| x.abs())
                    - median_filtered_difference
                        .slice(slice)
                        .mapv(|x| x.abs());
                let reject = reject * noise_rev;
                let reject = reject.mapv(|x| {
                    if x >= 50. {x} else {0.}});

                // Get the max values of this slice and keep the coordinates

                reject_value[[irow]] = *reject.max().unwrap();
                reject_index[[irow]] = begin + reject.argmax().unwrap() as isize;
            }

            let values_max = reject_value.max().unwrap();
            let row_max = reject_value.argmax().unwrap();
            let col_max = reject_index[[row_max]] as usize;

            // If max == 0 remove the line from those that are considered,
            // else we flag the value in the mask

            if *values_max != 0. {

                // We exclude this pixel from the map.
                
                mask[[row_max, col_max]] = 0;
                next_lines_to_consider.push(ispectrum);
            }
        }
        lines_to_consider = next_lines_to_consider;
    }

    
    

    let mask = mask.mapv(|x| x as f64);


    if mask_path.exists() {
        fs::remove_file(mask_path).unwrap();
    }

    let mut fitsfile = FitsFile::create(mask_path)
        .with_custom_primary(&description)
        .open()
        .unwrap();
    let hdu = fitsfile.primary_hdu().unwrap();
    hdu.write_image( &mut fitsfile, mask.clone().as_slice().unwrap()).unwrap();

    println!("[{:.1?}]", now.elapsed());
    
}








// We apply a simple median filter on the input array where
// we reflect the array when we are out of the boundaries.
// We assume that size is odd.
// Median filter over one column with a isize

fn median_filter(input : &Array2<f64>, size : isize) -> Array2<f64>
{

    let shape = input.shape();
    let num_rows = shape[0];
    let num_columns = shape[1];
    let mut output : Array2<f64> = Array::zeros((num_rows, num_columns));

    let half_size = (((size - 1) / 2) + 1) as usize;

    // We apply the filter on one column.

    for column in 0..num_columns {

        // Select the pixels that we consider for row 0
        
        let slice1 = s![0..half_size, column];
        let slice2 = s![0..half_size-1;-1, column];
        let selected_values = concatenate![Axis(0), input.slice(slice2), input.slice(slice1)];

        // We sort these values so that the median values of these selected values are found. 

        let ( mut selected_values, _idxs) = sort(selected_values);
        let mut median = selected_values[half_size-1];       

        for row in 0..num_rows {

            output[[row, column]] = median;

            // We drop a value from the selected_values.
            let to_drop : f64;
            if row < half_size-1 {
                to_drop = input[[half_size-2-row,column]];
            }
            else {
                to_drop = input[[row+1-half_size, column]];
            }

            // We add a value from the selected_values (such that the array is still sorted).
            let to_add: f64;
            if row >= num_rows - half_size {
                to_add = input[[2*num_rows - row - half_size-1,column]];
            } else {
                to_add = input[[ row + half_size, column]];
            }

            selected_values = drop(selected_values.clone(), to_drop).unwrap();
            selected_values = insert(&selected_values, to_add);

            median = selected_values[half_size-1];

        }
    }

    output
}





// This function sorts the input array using the merge sort algorithm.
// It returns a sorted array and an array with the indices of where the values
// of the sorted can be found in the original array.
// eg if: (values, indices) = sort(input), then: value[i] = input[ indices[i] ] forall i in range(values)

fn sort(input : Array1<f64>) -> (Array1<f64>, Array1<i64>)
{
    let size = input.shape()[0];

    if size == 1 {
        
        return (input, Array1::zeros(1));
    }

    let half_size = size / 2;
    let (input_lower, lower_indices) = sort(input.slice(s![0..half_size]).into_owned());
    let (input_upper, upper_indices) = sort(input.slice(s![half_size..size]).into_owned());
    merge(input_lower, input_upper, lower_indices, upper_indices)

}



// We merge two sorted lists with corresponding indices.
// The sorting should respect the order of these indices.

fn merge(input_lower : Array1<f64>, input_upper : Array1<f64>, lower_indices : Array1<i64>,
    upper_indices : Array1<i64>) -> (Array1<f64>, Array1<i64>)
{
    // We set up the output arrays.
    let combined_sizes = input_lower.shape()[0] + input_upper.shape()[0];
    let mut output_array: Array1<f64> = Array::zeros(combined_sizes);
    let mut output_indices: Array1<i64> = Array::zeros( combined_sizes);

    // Set up the index of both arrays.
    let mut idx1 = 0;
    let mut idx2 = 0;

    for i in 0..combined_sizes {
        
        // If we have reached the end of a input, we simply add the rest of the
        // other input array to the output array.
        if idx1 == input_lower.shape()[0] {
            output_array[i] = input_upper[idx2];
            output_indices[i] = input_lower.shape()[0] as i64 + upper_indices[idx2];
            idx2 = idx2+1;
        } else if idx2 == input_upper.shape()[0] {
            output_array[i] = input_lower[idx1];
            output_indices[i] = lower_indices[idx1];
            idx1 = idx1+1;
        } else {            
            if input_lower[idx1] <= input_upper[idx2] {
                output_array[i] = input_lower[idx1];
                output_indices[i] = lower_indices[idx1];
                idx1 = idx1+1;
            }
            else {
                output_array[i] = input_upper[idx2];
                output_indices[i] = input_lower.shape()[0] as i64+upper_indices[idx2];
                idx2 = idx2+1;
            }
        }
    }
    (output_array, output_indices)
}



// This function drops a value from an input array and returns an
// array like the input array without the value. 
fn drop(input_array: Array1<f64>, value: f64) -> Option<Array1<f64>> {

    let length_array = input_array.shape()[0];

    // If array is empty, no value can be dropped.
    
    if length_array == 0 { return None;}

    // If array only has one element, that element should be the value.

    if length_array == 1 {
        if input_array[0] == value {
            return Some(array![]);
        }
        else {
            return None;
        }
    }

    let middle_index = length_array / 2;

    // The to-be-dropped value is the middle value of the array.

    if input_array[middle_index] == value {
        let lower_slice = s![0..middle_index];
        let upper_slice = s![middle_index+1..length_array];
        let dropped_array = concatenate![Axis(0), input_array.slice(lower_slice), input_array.slice(upper_slice)];
        return Some(dropped_array);
    }

    // If the to-be-dropped value is larger then the middle value.
    // This means that the value is on the right side of the sorted array.

    else if input_array[middle_index] < value {
        let lower_slice = s![0..middle_index+1];
        let upper_slice = s![middle_index+1..length_array];
        let dropped_part = drop(input_array.slice(upper_slice).to_owned(), value);
        let other_part = input_array.slice(lower_slice);

        match dropped_part {
            None => return None,
            Some(i) => return Some(concatenate![Axis(0), other_part, i]),}

    }

    // If the to-be-dropped is smaller then the middle value.
    // This means that the value is on the left side of the sorted array.

    else {
        let lower_slice = s![0..middle_index];
        let upper_slice = s![middle_index..length_array];
        let dropped_part = drop(input_array.slice(lower_slice).to_owned(), value);
        let other_part = input_array.slice(upper_slice);

        match dropped_part {
            None => return None,
            Some(i) => return Some(concatenate![Axis(0), i, other_part]),
        };
    }
}



// Inerst a value in a sorted list
fn insert(input_array: &Array1<f64>, value: f64) -> Array1<f64> {
    let length_array = input_array.shape()[0];
    if length_array == 0 {
        return array![value];
    }

    let middle_index = length_array / 2;

    // If the value to be added is the same as the value on the middle index.
    if input_array[middle_index] == value {
        let lower_slice = s![0..middle_index+1];
        let upper_slice = s![middle_index+1..length_array];

        return concatenate![Axis(0),
            input_array.slice(lower_slice),
            array![value],
            input_array.slice( upper_slice)];
    // If the value to be added is larger than the value on the middle index.
    // The value should be added on the right side of the array.
    } else if input_array[middle_index] < value {
        let lower_slice = s![0..middle_index+1];
        let upper_slice = s![middle_index+1..length_array];

        return concatenate![Axis(0),
            input_array.slice(lower_slice),
            insert(&input_array.slice( upper_slice).to_owned(), value)];
    // If the value to be added is smaller than the value on the middle index.
    // The value should be added on the left side of the array.
    } else {
        let lower_slice = s![0..middle_index];
        let upper_slice = s![middle_index..length_array];
        return concatenate![Axis(0),
            insert(&input_array.slice(lower_slice).to_owned(),value),
            input_array.slice( upper_slice)];       
    }    
}
