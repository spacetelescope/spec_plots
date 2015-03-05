#!/usr/bin/env python

__version__ = '1.32.0'

"""
.. module:: make_hst_spec_previews

   :synopsis: Given an HST spectrum file name, this script will read in the data and generate preview plots of the spectra.  The plots are generated in different dimensions (large, thumbnail) and, depending on the instrument/configuration, are plotted in different formats to maximize readability and usability.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

import argparse
from astropy.io import fits
import numpy
from os import path
""" These are local modules that are imported. """
from utils import specutils, specutils_cos, specutils_stis

flux_scale_factor_default = 10.
fluxerr_scale_factor_default = 5.
n_consecutive_default = 20
output_path_default = ""
output_type_default = "png"
dpi_val_default = 96.
debug_default = False
full_ylabels_default = False
nooptimize_default = False
verbose_default = False

#--------------------

class HSTSpecPrevError(Exception, object):
    """
    This class defines a generic Exception to use for errors raised in MAKE_HST_SPEC_PREVIEWS.  It simply prints the given string when raising the exception. e.g.,

    .. code-block:: python
     
         raise HSTSpecPrevError("Print this string")
         HSTSpecPrevError: *** MAKE_HST_SPEC_PREVIEWS ERROR: 'Print this string'
    """

    def __init__(self, value):
        """
        Initiate the Exception.

        :param value: The string to print on error.

        :type value: str
        """
        self.value = value

    def __str__(self):
        """
        Overrides the str function for this class.
        """
        return "*** MAKE_HST_SPEC_PREVIEWS ERROR: "+repr(self.value)

#--------------------

def check_input_options(args):
    """
    Check that input arguments satisfy some minimum requirements.

    :param args: Stores arguments and options.

    :type args: argparse.Namespace object.

    :raises: HSTSpecPrevError
    """

    """ Make sure the input file is specified on input (non-empty string), and remove any leading/trailing whitespace.  Catch as an ArgumentsParser error in this case.  If it does exist, check that the input file exists at the time of the command-line parameter checking. """
    if not args.input_file:
        raise HSTSpecPrevError("File name must be specified.")
    else:
        args.input_file = args.input_file.strip()

    """ Make sure the output_type string is trimmed and lowercase. """
    args.output_type = args.output_type.strip().lower()

    """ The DPI value must be greater than zero... """
    if args.dpi_val < 0.:
        raise ValueError("DPI value must be > 0.")

#--------------------

def get_instrument_name(input_file):
    """
    Retrieves the instrument name from a FITS file primary header.

    :param input_file: Name of input FITS file.

    :type input_file: str

    :returns: str -- The name of the instrument based on the FITS header keyword, in uppercase with leading/trailing whitespace removed.

    :raises: KeyError, HSTSpecPrevError
    """

    """ Make sure file exists. """
    if not path.isfile(input_file):
        raise HSTSpecPrevError("Input file not found, looking for file: " + input_file)

    with fits.open(input_file) as hdulist:
        """ Make sure the INSTRUME keyword exists in the primary header, otherwise, catch as a KeyError in this case. """
        try:
            this_instrument = hdulist[0].header['INSTRUME']
        except KeyError:
            print "*** MAKE_HST_SPEC_PREVIEWS ERROR: INSTRUME keyword not found in this file's primary header: ", input_file
            exit(1)
    return this_instrument.strip().upper()

#--------------------

def make_hst_spec_previews(input_file, flux_scale_factor=flux_scale_factor_default, fluxerr_scale_factor=fluxerr_scale_factor_default, n_consecutive=n_consecutive_default, output_path=output_path_default, output_type=output_type_default, dpi_val=dpi_val_default, debug=debug_default, full_ylabels=full_ylabels_default, optimize=not nooptimize_default, verbose=verbose_default):
    """
    Main function in the module.

    :param input_file: The full path and name of the FITS file to create a preview for.
    
    :type input_file: str

    :param flux_scale_factor: The ratio between the flux and the median flux that defines the pass/fail criterion within the edge trim test.

    :type flux_scale_factor: float

    :param fluxerr_scale_factor: The ratio between the flux uncertainty and the median flux uncertainty that defines the pass/fail criterion within the edge trim test.

    :type fluxerr_scale_factor: float

    :param n_consecutive: The number of consecutive data points that must pass the edge trim test to define the start and end of the spectrum for plotting purposes.

    :type n_consecutive: int

    :param output_path: Full path to output plot files.  Do not include the output file name in this path.  If not supplied, plots will be created in the same directory as the input file.

    :type output_path: str

    :param output_type: The file type of the plots to make.

    :type output_type: str

    :param dpi_val: The DPI value of your device's monitor, which will affect the size of the output plots.

    :type dpi_val: float

    :param debug: If True, turns on debug mode, which will plot to the screen and color-code fluxes based on different rejection criteria.

    :type debug: bool

    :param full_ylabels: If True, label the y-axis with full values, including powers of ten in scientific notation.

    :type full_ylabels: bool

    :param optimize: If set to True, will use a slightly optimized version of determining the plot covering fraction.

    :type optimize: bool

    :param verbose: Turn on verbose messages/logging?

    :type verbose: bool

    :raises: HSTSpecPrevError
    """

    """ Print file name, if verbose is turned on. """
    if verbose:
        print "Input file: " + input_file

    """ Derive output file name from input file name. """
    if output_type != "screen":
        output_file = path.join(output_path,"") + path.basename(input_file).split(".fits")[0] + "." + output_type
    else:
        output_file = None

    """ Print name of output file, if verbose is turned on and not plotting to screen. """
    if verbose and output_type != "screen":
        print "Output file: " + output_file
    elif verbose:
        print "Output file: Plotting to screen."

    """ Read in the FITS file to determine which instrument it comes from.  Print the name of the instrument found in the header if verbose is turned on. """
    this_instrument = get_instrument_name(input_file)
    if verbose:
        print "Instrument: " + this_instrument

    """ Read in the FITS files and create plots using the local package appropriate for the instrument used in the input file. """
    if this_instrument == 'COS':
        """ Get wavelengths, fluxes, flux uncertainties. """
        cos_spectrum = specutils_cos.readspec(input_file)

        """ Get a list of segment names sorted such that the bluest segment is first. """
        cos_segment_names = specutils_cos.get_segment_names(cos_spectrum)

        """ Trim the wavelengths < 900 Angstroms for FUVB segment if optical element used is G140L. """
        if "FUVB" in cos_spectrum.segments and cos_spectrum.optical_element == "G140L":
            specutils_cos.extract_subspec(cos_spectrum, "FUVB", min_wl=900.)

        """ Create a stitched spectrum for use when making thumb-size plots. """
        stitched_spectrum = specutils.stitch_components(cos_spectrum, n_consecutive, flux_scale_factor, fluxerr_scale_factor, segment_names=cos_segment_names)

        """ Calculate plot metrics for the each segment. """
        segment_plot_metrics = [specutils.calc_plot_metrics("cos", cos_spectrum.segments[x].wavelengths, cos_spectrum.segments[x].fluxes, cos_spectrum.segments[x].fluxerrs, cos_spectrum.segments[x].dqs, n_consecutive, flux_scale_factor, fluxerr_scale_factor) for x in cos_segment_names]

        """ Make "large-size" plot. """
        specutils_cos.plotspec(cos_spectrum, output_type, output_file, n_consecutive, flux_scale_factor, fluxerr_scale_factor, segment_plot_metrics, dpi_val=dpi_val, output_size=1024, debug=debug, full_ylabels=full_ylabels, title_addendum=stitched_spectrum["title"], optimize=optimize)

        if not debug:
            """ Calculate plot metrics for the stitched spectrum. """
            stitched_plot_metrics = [specutils.calc_plot_metrics("cos", stitched_spectrum["wls"], stitched_spectrum["fls"], stitched_spectrum["flerrs"], stitched_spectrum["dqs"], n_consecutive, flux_scale_factor, fluxerr_scale_factor)]

            """ Make "thumbnail-size" plot, if requested. """
            specutils_cos.plotspec(cos_spectrum, output_type, output_file, n_consecutive, flux_scale_factor, fluxerr_scale_factor, stitched_plot_metrics, dpi_val=dpi_val, output_size=128, stitched_spectrum=stitched_spectrum, optimize=optimize)

    elif this_instrument == 'STIS':
        """ Get wavelengths, fluxes, flux uncertainties. """
        stis_spectrum = specutils_stis.readspec(input_file)

        """ Get the indices to plot.  If the number of associations is <=3 then we will plot all of them, but if >3 then only the first, middle, and last association will be plotted, so it's not necessary to stitch or calculate plot metrics for any of the others. """
        indices_to_plot = specutils_stis.get_association_indices(stis_spectrum.associations)

        """ Create a stitched spectrum for each association, used when making thumb-size plots.  The stitched spectrum joins the orders together (if there is more than one).  The length of the returned list is equal to the number of associations that will be plotted. """
        stitched_spectra = [specutils.stitch_components(x, n_consecutive, flux_scale_factor, fluxerr_scale_factor) for x in numpy.asarray(stis_spectrum.associations)[indices_to_plot]]

        """ Calculate plot metrics for the each association that will be plotted.  The length of the returned list is equal to the number of associations that will be plotted. """
        association_plot_metrics = [specutils.calc_plot_metrics("stis", x["wls"], x["fls"], x["flerrs"], x["dqs"], n_consecutive, flux_scale_factor, fluxerr_scale_factor) for x in stitched_spectra]

        """ Make "large-size" plot. """
        specutils_stis.plotspec(stis_spectrum, indices_to_plot, stitched_spectra, output_type, output_file, n_consecutive, flux_scale_factor, fluxerr_scale_factor, association_plot_metrics, dpi_val=dpi_val, output_size=1024, debug=debug, full_ylabels=full_ylabels, optimize=optimize)

        """ Make "thumbnail-size" plot, if requested.  Notice that in this case we always plot just the first association, by passing only `stitched_spectra[0]`. """
        if not debug:
            specutils_stis.plotspec(stis_spectrum, indices_to_plot, [stitched_spectra[0]], output_type, output_file, n_consecutive, flux_scale_factor, fluxerr_scale_factor, association_plot_metrics, dpi_val=dpi_val, output_size=128, optimize=optimize)

    else:
        raise HSTSpecPrevError('"INSTRUME" keyword not understood: ' + this_instrument)

#--------------------

def setup_args():
    """
    Set up command-line arguments and options.

    :returns: ArgumentParser -- Stores arguments and options.
    """
    parser = argparse.ArgumentParser(description="Create spectroscopic preview plots given an HST spectrum FITS file.")

    parser.add_argument("input_file", action="store", type=str, help="[Required] Full path to input file (HST spectrum) for which to generate preview plots.  Include the file name in the path.")

    parser.add_argument("-d", action="store_true", dest="debug", default=debug_default, help='[Optional] Turn on debug mode, which will plot to the screen and color-code fluxes based on different rejection criteria.  Default = %(default)s.')

    parser.add_argument("--dpival", action="store", type=float, dest="dpi_val", default=dpi_val_default, help="[Optional] Specify the DPI value of your device's monitor, which will affect the size of the output plots.  Default = %(default)s.")

    parser.add_argument("-e", action="store", type=float, dest="fluxerr_scale_factor", default=fluxerr_scale_factor_default, help="[Optional] Specify the ratio between the flux uncertainty and the median flux uncertainty that defines the pass/fail criterion within the edge trim test.  Default = %(default)s.")

    parser.add_argument("-n", action="store", type=int, dest="n_consecutive", default=n_consecutive_default, help="[Optional] Specify the number of consecutive data points that must pass the edge trim test to define the start and end of the spectrum for plotting purposes.  Default = %(default)i.")

    parser.add_argument("--nooptimize", action="store_true", dest="nooptimize", default=nooptimize_default, help='[Optional] If True, do not use an optimized method of calculating plot covering fraction to determine line transparency?  Default = %(default)s.')

    parser.add_argument("-o", action="store", type=str, dest="output_path", default=output_path_default, help="[Optional] Full path to output plot files.  Do not include file name in path.    Default = %(default)s, the same directory as the input file.", metavar='output path')

    parser.add_argument("-s", action="store", type=float, dest="flux_scale_factor", default=flux_scale_factor_default, help="[Optional] Specify the ratio between the flux and the median flux that defines the pass/fail criterion within the edge trim test.  Default = %(default)s.")

    parser.add_argument("-t", action="store", type=str, dest="output_type", default=output_type_default, help='[Optional] Specify the file type of the plots to make.  Default = %(default)s.', choices=['png','PNG','eps', 'EPS', 'screen','SCREEN'], metavar='{png,ps,screen}')

    parser.add_argument("-v", action="store_true", dest="verbose", default=verbose_default, help='[Optional] Turn on verbose messages/logging?  Default = %(default)s.')

    parser.add_argument("-y", action="store_true", dest="full_ylabels", default=full_ylabels_default, help='[Optional] Label the y-axis with full values, including powers of ten in scientific notation?  Default = %(default)s.')

    return parser

#--------------------

if __name__ == "__main__":

    """ Create ArgumentParser object that holds arguments and options. """
    args = setup_args().parse_args()

    """ Check arguments and options. """
    check_input_options(args)

    """ Call main function. """
    make_hst_spec_previews(args.input_file, \
                               flux_scale_factor = args.flux_scale_factor, \
                               fluxerr_scale_factor = args.fluxerr_scale_factor, \
                               n_consecutive = args.n_consecutive, \
                               output_path = args.output_path, \
                               output_type = args.output_type, \
                               dpi_val = args.dpi_val, \
                               debug = args.debug, \
                               full_ylabels = args.full_ylabels, \
                               optimize = not args.nooptimize, \
                               verbose = args.verbose)

#--------------------
