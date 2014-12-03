this folder contains the earthquake (EQ) epicentre databases

*.mat files are the raw data read into a MATLAB structure
*_prob.mat files are probabilistic sets (see eq_global_probabilistic)

isc-gem-cat.csv from www.isc.ac.uk/iscgem (main page) and directly from www.isc.ac.uk/iscgem/download.php and there http://colossus.iris.washington.edu/iscgem/download/isc-gem-cat.zip
The ISC-GEM Global Instrumental Earthquake Catalogue (1900-2009) is the result of a special effort to adapt and substantially extend and improve currently existing bulletin data to serve the requirements of the specific user group who assess and model seismic hazard and risk. See also www.globalquakemodel.org/gem/terms/licensing

signigeq_Oct_13_2014.xlsx from http://www.ngdc.noaa.gov/nndc/struts/form?t=101650&s=1&d=1,
see the code eq_signigeq_read. signigeq_CLEAN is the version climada reads, since the columns YEAR, MONTH, DAY, HOUR, MINUTE, SECOND, INTENSITY and LATITUDE and LONGITUDE need to be defined as Numeric and the file then saved as signigeq_CLEAN(with extension .xlsx, depends on the machine whether shown (Windws) or not (Mac))

centennial_Y2K.CAT.txt from http://earthquake.usgs.gov/data/centennial/centennial_Y2K.CAT, see also centennial.pdf and centennial_README.rtf in docs and the code eq_centennial_read
