# GSI-Spectral-Analyser
A simple Python script (Spectral-Analyser.py) to show the spectral absorption in the IR spectrum within the band windows of a given instrument.

The requirements to run this function are:

- Instrument band csv file
- Spectral data (from http://www.spectraplot.com/absorption) within a folder, where file names represent gas species

The spectral data in the IR and NIR is available on DropBox in `James/GSI Data/Spectroscropy`.

## Instrument data

The instrument band csv file must have columns representing:

`band id (or name), wavelength (nm)`

Accompanied by a spectral resolution argument in the command line when the program is launched.

Or if the band window is given in the csv (not just the center):

`band id (or name), wavelength (nm) start, wavelength (nm) end`

In which case the spectral resolution argument is not required.


## Creating a virtual environment

An `environment.yml` file has been created for use in conda to generate a working virtual environment.
To build the virtual environment run the following command from within the project folder:

`conda env create -f environment.yml`

Following by the following command to enter the virtual environment:

`activate spectralAnalyser`

## Launching the script manually

The script is launched at the command line using:

`python Spectral-Analyser.py --spectral_data_path="./IR and NIR/" --instr_bands_file="Orbita Bands OHS.csv" --spectral_res="2.5" --gas="CH4"`

Note that dependancies will need to be installed



Upon running the script, an `output` folder will be created which will contain imagery and an absorbance matrix as a csv

