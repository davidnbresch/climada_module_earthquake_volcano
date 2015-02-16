this folder contains the earthquake (EQ) epicentre databases

*.mat files are the raw data read into a MATLAB structure
*_prob.mat files are probabilistic sets (see eq_global_probabilistic)

isc-gem-cat.csv from www.isc.ac.uk/iscgem (main page) and directly from www.isc.ac.uk/iscgem/download.php and there http://colossus.iris.washington.edu/iscgem/download/isc-gem-cat.zip
The ISC-GEM Global Instrumental Earthquake Catalogue (1900-2009) is the result of a special effort to adapt and substantially extend and improve currently existing bulletin data to serve the requirements of the specific user group who assess and model seismic hazard and risk. See also www.globalquakemodel.org/gem/terms/licensing
Note: In order to generate an EQ hazard event set for Switzerland, some historic EQs from the ECOS-09 (the Earthquake Catalog of the Swiss Seismological Service) have been added manually to isc-gem-cat.csv.
See http://www.seismo.ethz.ch/prod/catalog/index for details on that catalog.
The following entries of the ECOS-09 catalog have been used:
date,	                lat,	 lon,   mw,  intens., region
1356/10/18 21:--:--, 47.47, 7.6,   6.6,  IX,     	Basel
1811/06/07 13:45:00, 46.85, 9.53,  3.9,  V, 	 Chur
1855/07/26 13:20:00, 46.23, 7.82,  5.3,  VII,    	Stalden-Visp
1917/06/20 23:09:00, 47.61, 9.1,   4.1,  V,      Engwilen
1946/05/30 00:35:00, 46.3,  7.52,  5.0,  VI,     Sierre
1961/08/09 13:--:--, 46.81, 10.35, 4.9,  V,      Scuol
1964/02/17 12:20:00, 46.88, 8.27,  4.8,  VII,    Sarnen
1971/09/29 07:18:52, 47.0,  9.0,   4.9,  VI,     Glarus
1992/05/00 06:44:40, 47.145,9.518, 4.34, IV-V,   Buchs

When redownloading the ISC-GEM Catalog, these entries are thus not contained in isc-gem-cat.csv by default - re-attach them manually in case you want to do an EQ risk assessment for Switzerland. All the data needed are given above, just put in blanks for parameters that are specified in the ISC-GEM catalog but not in ECOS-09, and use a value between 0km and 10km for 'depth' (most entries in ECOS-09 do not have a depth specified, but usually shallow EQs can be assumed for Switzerland). Epicentral intensity and region can be ignored (they are not among the parameters listed in the ISC-GEM Catalog). 


signigeq_Oct_13_2014.xlsx from http://www.ngdc.noaa.gov/nndc/struts/form?t=101650&s=1&d=1,
see the code eq_signigeq_read. signigeq_CLEAN is the version climada reads, since the columns YEAR, MONTH, DAY, HOUR, MINUTE, SECOND, INTENSITY and LATITUDE and LONGITUDE need to be defined as Numeric and the file then saved as signigeq_CLEAN(with extension .xlsx, depends on the machine whether shown (Windws) or not (Mac))

centennial_Y2K.CAT.txt from http://earthquake.usgs.gov/data/centennial/centennial_Y2K.CAT, see also centennial.pdf and centennial_README.rtf in docs and the code eq_centennial_read
