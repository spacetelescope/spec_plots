############################################################################################
## MAKE_HST_SPEC_PREVIEWS
##
## Given an HST spectrum file name, this script will read in the data and generate preview plots of the spectra.  The plots are generated in different dimensions (large, medium, small, thumbnail) and, depending on the instrument/configuration, are plotted in different formats to maximize readability and usability.
##
############################################################################################



############################################################################################
## Place import commands and logging options.
############################################################################################
import sys
import os
import logging
import pyfits
import specutils_cos # Local package.
import specutils_stis # Local package.
from optparse import OptionParser

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
############################################################################################



############################################################################################
## This class defines a generic Exception to use for errors raised in MAKE_HST_SPEC_PREVIEWS and specific to this module.  It simply returns the given value when raising the exception, e.g., raise HSTSpecPrevError("Print this string") -> __main__.MyError: 'Print this string.'
############################################################################################
class HSTSpecPrevError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
############################################################################################



############################################################################################
## This function sets up command-line options and arguments through the OptionParser class.
############################################################################################
def setup_input_options(parser):
    parser.add_option("-f", action="store", type="string", dest="input_file", default="", help="[Required] Full path to input file (HST spectrum) for which to generate preview plots.  Include the file name in the path.")
    parser.add_option("-o", action="store", type="string", dest="output_path", default="", help="[Required] Full path to output file (HST spectrum) for which to generate preview plots.  Do not inclue file name in path.")
    parser.add_option("-v", action="store_true", dest="verbose", default=False, help='[Optional] Turn on verbose messages/logging.  Default = "False".')
    parser.add_option("-t", action="store", dest="output_type", default="png", help='[Optional] Specify where plots should be output.  Options are "png" or "screen".  Default = "png".')
############################################################################################



############################################################################################
## This function checks input arguments satisfy some minimum requirements.
############################################################################################
def check_input_options(parser,opts):
    ## Make sure the input file is specified on input (non-empty string).  Catch as an OptionsParser error in this case.
    if not opts.input_file.strip():
        parser.error("File name must be specified.")

    ## Make sure input file exists on disk, and is readable.  Catch as an IOError in this case.
    try:
        with open(opts.input_file) as ifile:
            pass
    except IOError as input_error:
        print "File could not be read, or does not exist: ", opts.input_file

    ## Make sure the output path is specified on input (non-empty string).  Catch as an OptionsParser error in this case.
    if not opts.output_path.strip():
        parser.error("Output path must be specified.")

    ## Make sure the output type is understood.
    if opts.output_type != "png" and opts.output_type != "screen":
        parser.error('Invalid choice of output type (-o option).  Valid options are "png" or "screen".')

    ## Make sure the output_type string is trimmed and lowercase.
    opts.output_type = opts.output_type.strip().lower()
############################################################################################



############################################################################################
## This function determines the instrument name from a FITS file primary header.
############################################################################################
def get_instrument_name(input_file):
    hdulist = pyfits.open(input_file)
    ## Make sure the INSTRUME keyword exists in the primary header, otherwise, catch as a PYFITS KeyError in this case.
    try:
        this_instrument = hdulist[0].header['INSTRUME']
    except KeyError as header_error:
        ## Make sure we close the FITS file before exiting with an Exception
        hdulist.close()
        print "INSTRUME keyword not found in file's primary header: ", input_file
    ## Close the FITS file.
    hdulist.close()
    return this_instrument.strip().upper()
############################################################################################



############################################################################################
## This is the main routine.
############################################################################################
def main():
    ## Define input options here.
    parser = OptionParser(usage="%prog -f <input file> -o <output path> [-v] [-t <plot type {png,screen}>]", version="%prog 1.0")
    setup_input_options(parser)

    ## Parse input options from the command line.
    opts, args = parser.parse_args()

    ## Check input arguments are valid and sensible.
    check_input_options(parser,opts)

    ## Print file name, if verbose is turned on.
    if opts.verbose:
        print "Input file: " + opts.input_file

    ## Derive output file name from input file name.
    if opts.output_type == "png":
        output_file = os.path.join(opts.output_path,"") + os.path.basename(opts.input_file).split(".fits")[0] + ".png"

    ## Print name of output file, if verbose is turned on and not printing to screen.
    if opts.verbose and opts.output_type != "screen":
        print "Output file: " + output_file
    else:
        print "Output file: Plotting to screen."

    ## Read in the FITS file to determine which instrument it comes from.  Print the name of the instrument found in the header if verbose is turned on.
    this_instrument = get_instrument_name(opts.input_file)
    if opts.verbose:
        print "Instrument: " + this_instrument

    ## Read in the FITS file to extract wavelengths, fluxes, and flux uncertainties, using the local package appropriate for the instrument used in the input file.
    if this_instrument == 'COS':
        ## Get wavelengths, fluxes, flux uncertainties.
        spec_tuple = specutils_cos.readspec(opts.input_file, opts.verbose)
        ## Plot spectra preview plots.
        specutils_cos.plotspec(spec_tuple, opts.output_type, output_file)
    elif this_instrument == 'STIS':
        pass
        ## specutils_stis.readspec(opts.input_file)
        ## specutils_stis.plotspec()
    else:
        raise HSTSpecPrevError('"INSTRUME" keyword not understood: ' + this_instrument)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    logger.setLevel(logging.INFO)
    main()
############################################################################################
