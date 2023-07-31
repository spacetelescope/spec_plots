Change Log for SPEC_PLOTS
=========================

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
