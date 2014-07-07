__version__ = '1.0'

"""
.. module:: specutils_cos
   :synopsis: Contains functions for reading and plotting HST COS spectra.
.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

from numpy import array_equal
import os
import pyfits
import matplotlib.pyplot as pyplot
"""These are local modules that are imported."""
from make_hst_spec_previews import HSTSpecPrevError


def check_segments(segments_list, input_file):
    """
    Checks that the array of "segments" in the COS spectrum header are expected values.  It returns a scalar string representing the band (either "FUV" or "NUV").  If the array is not what's expected for COS, an Exception is raised.
    :param segments_list: List of segment labels.
    :type segments_list: list
    :param input_file: Name of input FITS file.
    :type input_file: str
    :returns: str -- A string representation of the band, either "FUV" or "NUV".
    :raises: HSTSpecPrevError
    .. note::
       This function does not attempt to access the input file, it only requires the name of the file for error reporting purposes.
    """
    list_len = len(segments_list)
    """Segment array must contain either two or three elements."""
    if list_len == 2:
        """Must be ["FUVA", "FUVB"]."""
        if array_equal(segments_list, ["FUVA", "FUVB"]):
            this_band = "FUV"
        else:
            try:
                raise HSTSpecPrevError("The array of SEGMENT strings contains 2 values, but is not equal to [\"FUVA\", \"FUVB\"] in file " + input_file)
            except HSTSpecPrevError as error_string:
                print error_string
                exit(1)
    elif list_len == 3:
        """Must be ["NUVA", "NUVB", "NUVC"]."""
        if array_equal(segments_list, ["NUVA", "NUVB", "NUVC"]):
            this_band = "NUV"
        else:
            try:
                raise HSTSpecPrevError("The array of SEGMENT strings contains 3 values, but is not equal to [\"NUVA\", \"NUVB\", \"NUVC\"] in file " + input_file)
            except HSTSpecPrevError as error_string:
                print error_string
                exit(1)
    else:
        try:
            raise HSTSpecPrevError("The array of SEGMENT strings should contain 2 or 3 values, found " + str(list_len) + " in file " + input_file)
        except HSTSpecPrevError as error_string:
                print error_string
                exit(1)
    return this_band


def readspec(input_file, verbosity=False):
    """
    Reads in a COS spectrum FITS file (x1d, x1dsum, or x1dsum{1,2,34} FITS files) and returns the wavelengths, fluxes, and flux uncertainties for the two (FUV segments) or three (NUV stripes).
    :param input_file: Name of input FITS file.
    :type input_file: str
    :param verbosity: Should verbose output be displayed?  Default is False.
    :type verbosity: bool
    :returns: tuple -- The spectroscopic data {band, wavelength table, flux table, flux uncertainty table}, where:
    band = Scalar string, either "FUV" or "NUV".
    wavelength table = Table of wavelengths for each segment/slice, where a given segment/slice is a row in the table.  FUV has 2 segments, NUV has 3 slices.
    flux table = Table of fluxes for each segment/slice, where a given segment/slice is a row in the table.  FUV has 2 segments, NUV has 3 slices.
    flux uncertainty table = Table of flux uncertainties for each segment/slice, where a given segment/slice is a row in the table.  FUV has 2 segments, NUV has 3 slices.
    :raises: KeyError
    """
    with pyfits.open(input_file) as hdulist:
        """Read the data from the first extension.  For COS, the spectra are always stored as tables in the first FITS extension."""
        cos_tabledata = hdulist[1].data
        """Extract the SEGMENTS.  This is either a 2-element array of ["FUVA", "FUVB"], or a 3-element array of ["NUVA", "NUVB", "NUVC"]."""
        try:
            segment_list = cos_tabledata.field("SEGMENT")
        except KeyError:
            print "*** MAKE_HST_SPEC_PREVIEWS ERROR: SEGMENT column not found in first extension's binary table."
            exit(1)
        band = check_segments(segment_list, input_file)
        """Extract the number of elements (n_wavelengths, n_fluxes, etc.) for each segment.  This will also be either a 2-element array (FUV) or 3-element array (NUV)."""
        try:
            nelems_list = cos_tabledata.field("NELEM")
        except KeyError:
            print "*** MAKE_HST_SPEC_PREVIEWS ERROR: NELEM column not found in first extension's binary table."
            exit(1)
        """Extract wavelength, fluxes, and flux uncertainties for each segment.  These will be either 2xn (FUV) or 3xn (NUV) tables."""
        try:
            wavelength_table = cos_tabledata.field("WAVELENGTH")
        except KeyError:
            print "*** MAKE_HST_SPEC_PREVIEWS ERROR: WAVELENGTH column not found in first extension's binary table."
            exit(1)
        try:
            flux_table = cos_tabledata.field("FLUX")
        except KeyError:
            print "*** MAKE_HST_SPEC_PREVIEWS ERROR: FLUX column not found in first extension's binary table."
            exit(1)
        try:
            fluxerr_table = cos_tabledata.field("ERROR")
        except KeyError:
            print "*** MAKE_HST_SPEC_PREVIEWS ERROR: ERROR column not found in first extension's binary table."
            exit(1)
        return (band, wavelength_table, flux_table, fluxerr_table)


def plotspec(this_spec_tuple, output_type, output_file):
    """
    This function accepts a COS spectrum tuple from the READSPEC function and produces preview plots.
    :param this_spec_tuple: COS spectrum as returned by READSPEC.
    :type this_spec_tuple: tuple
    :param output_type: What kind of output to make.
    :type output_type: str
    :param output_file: Name of output file (including full path).
    :type output_file: str
    :raises: OSError,HSTSpecPrevError
    """
    """Make sure the output path exists, if not, create it."""
    if not os.path.isdir(os.path.dirname(output_file)):
        try:
            os.mkdir(os.path.dirname(output_file))
        except OSError as this_error:
            if this_error.errno == 13: 
                print "*** MAKE_HST_SPEC_PREVIEWS ERROR: Output directory could not be created, "+repr(this_error.strerror)
                exit(1)
            else:
                raise
    """Extract the band from the tuple."""
    this_band = this_spec_tuple[0]
    """Start plot figure."""
    pyplot.figure(1)
    if this_band == "FUV":
        """Then there are two segments, so make a 2-panel plot.
        First plot area (top) is for *FUVB*."""
        segment_1_plotarea = pyplot.subplot(2, 1, 1)
        pyplot.plot(this_spec_tuple[1][0], this_spec_tuple[2][0], 'k')
        """Second plot area (bottom) is for *FUVA*."""
        segment_2_plotarea = pyplot.subplot(2, 1, 2)
        pyplot.plot(this_spec_tuple[1][1], this_spec_tuple[2][1], 'k')
        if output_type == "png":
            pyplot.savefig(output_file, format=output_type)
        elif output_type == "screen":
            pyplot.show()
    elif this_band == "NUV":
        pass
    else:
        try:
            raise HSTSpecPrevError("The band is not understood, expect either \"FUV\" or \"NUV\", found " + this_band + " in file " + input_file)
        except HSTSpecPrevError as error_string:
            print error_string
            exit(1)
