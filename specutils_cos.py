__version__ = '1.0'

"""
.. module:: specutils_cos
   :synopsis: Contains functions for reading and plotting HST COS spectra.
.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

import numpy
import os
import pyfits
import matplotlib.pyplot as pyplot
"""These are local modules that are imported."""
from make_hst_spec_previews import HSTSpecPrevError

class COSSpectrum:
    """
    Defines a COS spectrum, including wavelegnth, flux, and flux errors.  A COS spectrum consists of N segments (N = {2,3}) stored as a dict object.  Each of these dicts contain a COSSegment object that contains the wavelengths, fluxes, flux errors, etc.
    :raises: ValueError
    """
    def __init__(self, band=None, cos_segments=None):
        """
        Create a COSSpectrum object given a band choice (must be "FUV" or "NUV").
        :param band: Which band is this spectrum for ("FUV" or "NUV")?
        :type band: str
        :param cos_segments: [Optional] COSSegment objects to populate the COSSpectrum with.
        :type cos_segments: dict
        """
        if band.strip().upper() == "FUV":
            self.band = band
            if cos_segments is not None:
                if len(cos_segments) == 2:
                    self.segments = {'FUVA':cos_segments['FUVA'], 'FUVB':cos_segments['FUVB']}
                else:
                    raise ValueError("Band is specified as "+band.strip().upper()+", expected 2 COSSegment objects as a list but received "+str(len(cos_segments))+".")
            else:
                self.segments = {'FUVA':COSSegment(), 'FUVB':COSSegment()}
        elif band.strip().upper() == "NUV":
            self.band = band
            pass
        else:
            raise ValueError("Must specify band=\"FUV\" or band=\"NUV\".")


class COSSegment:
    """
    Defines a spectrum from a COS segment.  The data (wavelength, flux, flux errors) are stored as numpy ndarrays.  A scalar int property provides the number of elements in this segment.
    """
    def __init__(self, nelem=None, wavelengths=None, fluxes=None, fluxerrs=None):
        """
        Create a COSSegment object, default to empty values.  Allows user to preallocate space if they desire by setting "nelem" but not providing lists/arrays on input right away.
        :param nelem: Number of elements for this segment's
        """
        if nelem is not None:
            self.nelem = nelem
        else:
            self.nelem = 0
        if wavelengths is not None:
            self.wavelengths = numpy.asarray(wavelengths)
        else:
            self.wavelengths = numpy.zeros(self.nelem)
        if fluxes is not None:
            self.fluxes = numpy.asarray(fluxes)
        else:
            self.fluxes = numpy.zeros(self.nelem)
        if fluxerrs is not None:
            self.fluxerrs = numpy.asarray(fluxerrs)
        else:
            self.fluxerrs = numpy.zeros(self.nelem)


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
        if numpy.array_equal(segments_list, ["FUVA", "FUVB"]):
            this_band = "FUV"
        else:
            try:
                raise HSTSpecPrevError("The array of SEGMENT strings contains 2 values, but is not equal to [\"FUVA\", \"FUVB\"] in file " + input_file)
            except HSTSpecPrevError as error_string:
                print error_string
                exit(1)
    elif list_len == 3:
        """Must be ["NUVA", "NUVB", "NUVC"]."""
        if numpy.array_equal(segments_list, ["NUVA", "NUVB", "NUVC"]):
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
    :returns: COSSpectrum -- The spectroscopic data (wavelength, flux, flux error, etc):
    :raises: KeyError
    """
    with pyfits.open(input_file) as hdulist:
        """Read the data from the first extension.  For COS, the spectra are always stored as tables in the first FITS extension."""
        cos_tabledata = hdulist[1].data
        """Extract the SEGMENTS.  This is either a 2-element array of ["FUVA", "FUVB"], or a 3-element array of ["NUVA", "NUVB", "NUVC"]."""
        try:
            segment_arr = cos_tabledata.field("SEGMENT")
        except KeyError:
            print "*** MAKE_HST_SPEC_PREVIEWS ERROR: SEGMENT column not found in first extension's binary table."
            exit(1)
        band = check_segments(segment_arr, input_file)
        """Extract the number of elements (n_wavelengths, n_fluxes, etc.) for each segment.  This will also be either a 2-element array (FUV) or 3-element array (NUV)."""
        try:
            nelems_arr = cos_tabledata.field("NELEM")
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
        """Create COSSegment objects to populate the COSSpectrum object with"""
        fuva_index = numpy.where(segment_arr == "FUVA")[0][0]
        fuva_cossegment = COSSegment(nelem=nelems_arr[fuva_index], wavelengths=wavelength_table[fuva_index,:], fluxes=flux_table[fuva_index,:], fluxerrs=fluxerr_table[fuva_index,:])
        fuvb_index = numpy.where(segment_arr == "FUVB")[0][0]
        fuvb_cossegment = COSSegment(nelem=nelems_arr[fuvb_index], wavelengths=wavelength_table[fuvb_index,:], fluxes=flux_table[fuvb_index,:], fluxerrs=fluxerr_table[fuvb_index,:])
        """Create COSSpectrum object."""
        return_spec = COSSpectrum(band=band, cos_segments={'FUVA':fuva_cossegment,'FUVB':fuvb_cossegment})
#        return {band, segment_arr, nelems_arr, wavelength_table, flux_table, fluxerr_table}
        return return_spec

def set_plot_xrange(wavelengths,fluxes):
    """
    Given a COS tuple from READSPEC, returns a list of [xmin,xmax] to define an optimal x-axis plot range.
    :param wavelengths: The wavelengths to be plotted.
    :type wavelengths: numpy.ndarray
    :param fluxes: The fluxes to be plotted.
    :type fluxes: numpy.ndarray
    :returns: list -- Two-element list containing the optimal [xmin,xmax] values to define the x-axis plot range.
    """
    """Test if there are any NaN's in the wavelength array.  If so, issue a warning for now..."""
    """ NOTE: Use of the SUM here was reported on stackoverflow to be faster than MIN...it won't matter for the sizes we're dealing with here, but I thought it was a neat trick."""
    if numpy.isnan(numpy.sum(wavelengths)):
        print "***WARNING in SPECUTILS_COS: Wavelength array contains NaN values.  Behavior has not been fully tested in this case."
    """Sort the wavelength and flux arrays.  Don't do it in-place since we aren't (yet) going to sort the fluxes, errors, etc. as well."""
    sorted_indexes = numpy.argsort(wavelengths)
    sorted_wavelengths = wavelengths[sorted_indexes]
    sorted_fluxes = fluxes[sorted_indexes]
    """Find the first element in the array that is NOT 0.0, and the last element in the array that is NOT 0.0."""
    start_index = 0
    end_index = -1
    while sorted_fluxes[start_index] == 0.:
        start_index += 1
    while sorted_fluxes[end_index] == 0.:
        end_index -= 1
    """Return the optimal start and end wavelength values for defining the x-axis plot range."""
    return [wavelengths[start_index],wavelengths[end_index]]
    

def plotspec(cos_spectrum, output_type, output_file):
    """
    Accepts a COS spectrum tuple from the READSPEC function and produces preview plots.
    :param cos_spectrum: COS spectrum as returned by READSPEC.
    :type cos_spectrum: COSSpectrum
    :param output_type: What kind of output to make.
    :type output_type: str
    :param output_file: Name of output file (including full path).
    :type output_file: str
    :raises: OSError,HSTSpecPrevError
    .. note::
       This function assumes a screen resolution of 96 DPI in order to generate plots of the desired sizes.  This is because matplotlib works in units of inches and DPI rather than pixels.
    """
    dpi_val = 96.
    """Make sure the output path exists, if not, create it."""
    if output_type != 'screen':
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
    this_band = cos_spectrum.band
    """Start plot figure."""
    pyplot.figure(figsize=(800./dpi_val,600./dpi_val),dpi=dpi_val)
    if this_band == "FUV":
        """Then there are two segments, so make a 2-panel plot.
        First plot area (top) is for *FUVB*, since that's the bluest segment."""
        segment_1_plotarea = pyplot.subplot(2, 1, 1)
        """Determine optimal x-axis."""
        x_axis_range = set_plot_xrange(cos_spectrum.segments["FUVB"].wavelengths, cos_spectrum.segments["FUVB"].fluxes)
        """Plot the spectrum."""
        pyplot.plot(cos_spectrum.segments["FUVB"].wavelengths, cos_spectrum.segments["FUVB"].fluxes, 'k')
        """Update the x-axis range."""
        pyplot.xlim(x_axis_range)
        """Second plot area (bottom) is for *FUVA, since that's the reddest*."""
        segment_2_plotarea = pyplot.subplot(2, 1, 2)
        """Determine optimal x-axis."""
        x_axis_range = set_plot_xrange(cos_spectrum.segments["FUVA"].wavelengths, cos_spectrum.segments["FUVA"].fluxes)
        """Plot the spectrum."""
        pyplot.plot(cos_spectrum.segments["FUVA"].wavelengths, cos_spectrum.segments["FUVA"].fluxes, 'k')
        """Update the x-axis range."""
        pyplot.xlim(x_axis_range)
        """Display or plot to the desired format."""
        if output_type != "screen":
            pyplot.savefig(output_file, format=output_type, dpi=dpi_val,bbox_inches='tight')
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
