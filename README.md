# MARVEL pipeline
Data processing and radial velocity estimation pipeline of the MARVEL spectrograph


## Running the different components of the pipeline

The pipeline is made out of different components that can run one after the other. Every time a component finishes it will create a FITS file with the intermediate results saved into it. This file can then serve as input for another component. An example of how these components rely on the results of previous components is illustrated in the following graph:

![Components needed to create an "Optimal Science Extraction" image and how these are in turn dependent on other components.](./Docs/Images/my_output_file.png "Optimal Extraction")


## Dependencies

### Python and Rust 

The pipeline components are written in [Python](https://www.python.org/) or [rust](https://foundation.rust-lang.org/). 
In the following section, we explain how to install the Python packages using poetry. 
To [install rust](https://www.rust-lang.org/tools/install) we recommend 
following the official installation guide given on its website. 

### Python dependencies

Before starting to install the packages, we recommend using a software tool to manage different environments like [conda](https://docs.conda.io/projects/conda/en/stable/commands/create.html) to create an environment for marvel. Using conda, use the command:

```
conda create -n marvelpipeline [python=<python version>]
```

where the `<python-version>` should be higher than Python 3.9. We can then activate the environment with

```
conda activate marvelpipeline
```

To install the needed Python packages we recommend installing the package dependency manager [poetry](https://python-poetry.org/) using [the official installation](https://python-poetry.org/docs/). This tool allows us to install all the packages with one command in a reproducible way. We can install the dependencies with the command:

```
poetry install
```

This command will look for the `pyproject.toml` file in the directory, and install all the dependencies that are specified in that file.


## Compiling the pipeline

As the computationally heavy parts of the pipeline are written in Rust, you need to build the rust code before it can be run. This process a made easy by simply running the following bash shell script (to be executed in the home folder of the Marvel pipeline repository):

```
./compile
```


## Configuration files

#### The `inputfile.yaml` file

In `inputfile.yaml` you can specify which CCD image FITS files (raw bias images, raw dark images, raw flat images, and raw science images) should be processed. The paths specified are relative to a root path that is given later.


#### The `config.yaml` file

`config.yaml` is the inputfile that the pipeline will use to process the MARVEL images. It contains both the paths to the input FITS files as
well as the paths to the ouput file names. `config.yaml` can be auto-generated from `inputfile.yaml` using the
`configure.py` script:

```
python configure.py inputfile.yaml config.yaml <root folder to raw data> <root folder to processed data>
```

For example:
```
python configure.py inputfile.yaml config.yaml /home/marvel/rawdata/  /home/marvel/processeddata/

```

All folder locations of raw MARVEL data files mentioned in the config.yaml file are relative to the specified root folder of raw data.
Similarly, all folder locations of a pipeline product mentioned in the config.yaml file are relative to the specified root folder of
processed data.



## Running the pipeline

The pipeline can be run on the command line using `marvelpipe.py`. It has a convenient help function:

```
python marvelpipe.py -h
```

listing the different pipeline steps and the options:

```
usage: marvelpipe.py [-h] [-f FIRST] [-l LAST] configfile

MARVEL Data Reduction Pipeline Steps:

 1) Compute master bias image
 2) Compute master dark image
 3) Compute master flat image
 4) Compute master ThAr image
 5) Correct CCD images for bias and dark
 6) Determine 2D mask of the orders
 7) Extract 1D orders

positional arguments:
  configfile                YAML configuration file

options:
  -h, --help                Show this help message and exit
  -f FIRST, --first FIRST   First pipeline step to execute (default: 1)
  -l LAST, --last LAST      Last pipeline step to execute (default: 7)

Examples:
   $ python marvelpipe.py config.yaml              - Run all pipeline steps on the files specified in config.yaml
   $ python marvelpipe.py -f 3 config.yaml         - Run pipeline steps 3, 4, ... until the end
   $ python marvelpipe.py -f 3 -l 6 config.yaml    - Run pipeline steps 3, 4, 5, and 6

```

If you want to run the entire pipeline, simply run:

```
python marvelpipe.py config.yaml
```

If some of the pipeline steps were computed before, and you don't want to recompute them, you can specify
the first and/or the last pipeline step you want to execute. For example:

```
python marvelpipe.py -f 3 config.yaml
```
skips the first two steps and starts with the 3rd step, while:
```
python marvelpipe.py -f 3 -l 6 config.yaml

```
only executes steps 3, 4, 5, and 6.


## The pipeline output folders

The results of the pipeline are stored in the `ProcessedData/` directory inside the root folder you specified. It contains the following subdirectories:

- `MasterBias/`: contains the master bias FITS file: a median of several raw bias images.
- `MasterFlat/`: contains the master flat FITS file: a median of several raw flat images, bias subtracted.
- `BiasSubtractedScience/`: contains 2D science image FITS files, of which the master bias was subtracted.
- `ExtractedOrders/`:  this contains 2 subdirectories
	- `Mask/`:
		- `*_2d_mask.fits`: contains for each row of the CCD image the column where the order begins, the column where the order ends, and the column where the order has the maximum flux. This file is used to extract the 1d orders of each science frame.
		- `*_smoothed_master_flat.fits`: a smoothed version of the master flatfield image. This is the image that is used to determine the position of the orders.
	- `Science/`: for each bias subtracted science frame, the `*_2d_orders.fits` file contains a CCD image where all pixels outside the orders are set to zero. This file is only used for debugging and visualization purposes, and is not further used by the pipeline.
- `OptimalExtraction/`: for each bias subtracted science frame, the `*_1d_orders.fits` file contains the 1D flat-relative extracted orders. These 1D orders are the ratio of the stellar spectrum divided by the flatfield lamp spectrum. As both the stellar spectrum and the flatfield spectrum experience the same blaze function, the ratio should be free of the blaze function. For each order a flat-normalized flux as a function of pixel coordinate is stored. The orders are not merged.
- `WaveCalibration/`: TBD
 

## The pipeline output files

TBW


