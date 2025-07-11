Change Log for SPEC_PLOTS
=========================

v1.37.1 - 2025 July 2 
----------------------
* Updated readspec() for JWST to support new time series observations
  (TSO) where all exposures are in a single FITS extension called
  "EXTRACT1D". Previously each exposure was stored in its own
  extension, all called "EXTRACT1D".

v1.37.0 - 2025 June 17 
----------------------
* Switched to using `tostring_argb()` in `calc_covering_fraction`
  for compatibility with matplotlib > 3.10.0
* Removed a misleading comment in readspec to reflect that JWST are
  not always stored in the first extension, nor do they need to be for
  `spec_plots` to work.
* All error states end with a sys.error(1) error code now.
* Updated numpy.fromstring() to use numpy.frombuffer() for
  compatibility with numpy > v2.2.0
  
  
v1.36.1 - 2025 Feb. 14 
----------------------
* Fixed incorrect units being shown for JWST NIRSPEC 1D files  
* NIRCam x1d files now supported  
  
v1.36.0 - 2024 Feb. 9 
-----------------
* Added basic support HST Hubble Advanced Spectral Products (HASP)  
 
v1.35.1 - 2023 Aug. 4
-----------------
* Added support for flux uncertainties being in FLUX_ERROR column
  instead of ERROR column for JWST x1d FITS files.
* Fixed readability of y-axis tick labels by adjusting font size,
  angle, and padding.

v1.35.0 - 2023 Jul. 31
-----------------
* Added required config file for ReadTheDocs.
* Updated package dependencies to newer versions.

v1.34.6 - 2018 Mar. 15
-----------------
* Changed default backend back to Agg.

v1.34.5 - 2018 Mar. 15
-----------------
* Output file names when verbose turned on now accurate.
* Changed output file names, size no longer in the file name.
* Fixed multi-panel large spectra so labels weren't overlapping with other text.
* Changed matplotlib backend to be TkAGG instead of Agg.

v1.34.4 - 2018 Feb. 02
-----------------
* Fixed Python 2.x issue with output_type argparse option.
* Can now specify more than one output type from the command line.

v1.34.3 - 2017 Dec. 08
-----------------
* Added binary FITS table as an output option for COS spectra.

v1.34.2 - 2017 Mar. 15
-----------------
* Removed useless argument checks in main function.
* First build that includes a conda package.

v1.34.1 - 2016 Nov. 28
-----------------
* Added basic support for NIRSPEC and NIRISS.
* Updated ERR keyword to ERROR keyword for JWST instruments.
* Automated version numbers in API doc.
* Renamed some modules that were named after MIRI to be JWST (more generic).

v1.34 - 2016 Oct.
-----------------
* Added dual-support for Python 2.7 and Python 3.5.
* Added basic support for JWST MIRI 1D spectra.
