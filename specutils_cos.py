__version__ = '1.05'

"""
.. module:: specutils_cos
   :synopsis: Contains functions for reading and plotting HST COS spectra.
.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

from astropy.io import fits
import math
from matplotlib import rc
import matplotlib.pyplot as pyplot
import numpy
import os
"""These are local modules that are imported."""
from make_hst_spec_previews import HSTSpecPrevError

class COSAvoidRegion:
    """
    Defines a COS avoid region, which is simply a section of wavelength space that should not be included when determining the optimal y-axis plot range.  The object consists of a starting wavelength, ending wavelength, and string description of what that region is.
    :raises: ValueError
    """
    def __init__(self, minwl=None, maxwl=None, description=""):
        if minwl is None:
            raise ValueError("Must specify a minimum wavelength for this COS avoid region.")
        if maxwl is None:
            raise ValueError("Must specify a maximum wavelength for this COS avoid region.")
        if minwl >= maxwl:
            raise ValueError("Minimum wavelength must be less than maximum wavelength for this COS avoid region.  Given min. wavelength = "+str(minwl)+" and max. wavelength = "+str(maxwl)+".")
        self.minwl = minwl
        self.maxwl = maxwl
        self.description = description

class COSSpectrum:
    """
    Defines a COS spectrum, including wavelegnth, flux, and flux errors.  A COS spectrum consists of N segments (N = {2,3}) stored as a dict object.  Each of these dicts contain a COSSegment object that contains the wavelengths, fluxes, flux errors, etc.
    :raises: ValueError
    """
    def __init__(self, band=None, cos_segments=None, orig_file=None):
        """
        Create a COSSpectrum object given a band choice (must be "FUV" or "NUV").
        :param band: Which band is this spectrum for ("FUV" or "NUV")?
        :type band: str
        :param cos_segments: [Optional] COSSegment objects to populate the COSSpectrum with.
        :type cos_segments: dict
        :param orig_file: Original FITS file read to create the spectrum (includes full path).
        :type orig_file: str
        """
        self.orig_file = orig_file
        if band.strip().upper() == "FUV":
            self.band = band
            if cos_segments is not None:
                if len(cos_segments) == 2:
                    self.segments = {'FUVA':cos_segments['FUVA'], 'FUVB':cos_segments['FUVB']}
                elif len(cos_segments) == 1 and 'FUVA' in cos_segments:
                    self.segments = {'FUVA':cos_segments['FUVA']}
                elif len(cos_segments) == 1 and 'FUVB' in cos_segments:
                    self.segments = {'FUVB':cos_segments['FUVB']}
                else:
                    raise ValueError("Band is specified as "+band.strip().upper()+", expected 1 or 2 COSSegment objects as a list but received "+str(len(cos_segments))+".")
            else:
                self.segments = {'FUVA':COSSegment(), 'FUVB':COSSegment()}
        elif band.strip().upper() == "NUV":
            self.band = band
            if cos_segments is not None:
                if len(cos_segments) == 3:
                    self.segments = {'NUVA':cos_segments['NUVA'], 'NUVB':cos_segments['NUVB'], 'NUVC':cos_segments['NUVC']}
                else:
                    raise ValueError("Band is specified as "+band.strip().upper()+", expected 3 COSSegment objects as a list but received "+str(len(cos_segments))+".")
            else:
                self.segments = {'NUVA':COSSegment(), 'NUVB':COSSegment(), 'NUVC':COSSegment()}
        else:
            raise ValueError("Must specify band=\"FUV\" or band=\"NUV\".")


class COSSegment:
    """
    Defines a spectrum from a COS segment.  The data (wavelength, flux, flux errors) are stored as numpy ndarrays.  A scalar int property provides the number of elements in this segment.
    """
    def __init__(self, nelem=None, wavelengths=None, fluxes=None, fluxerrs=None, dqs=None):
        """
        Create a COSSegment object, default to empty values.  Allows user to preallocate space if they desire by setting "nelem" but not providing lists/arrays on input right away.
        :param nelem: Number of elements for this segment's spectrum.
        :type nelem: int
        :param wavelengths: List of wavelengths in this segment's spectrum.
        :type wavelengths: list
        :param fluxes: List of fluxes in this segment's spectrum.
        :type fluxes: list
        :param fluxerrs: List of flux uncertainties in this segment's spectrum.
        :type fluxerrs: list
        :param dqs: List of Data Quality (DQ) flags in this segment's spectrum.
        :type dqs: list
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
        if dqs is not None:
            self.dqs = numpy.asarray(dqs)
        else:
            self.dqs = numpy.zeros(self.nelem)


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
    """Segment array must contain either one, two, or three elements."""
    if list_len == 1:
        """Can happen with FUV data for some reason, so it can be either "FUVA" or "FUVB"."""
        if numpy.array_equal(segments_list, ["FUVA"]) or numpy.array_equal(segments_list, ["FUVB"]):
            this_band = "FUV"
        else:
            try:
                raise HSTSpecPrevError("The array of SEGMENT strings contains 1 value, but is not equal to [\"FUVA\"] or [\"FUVB\"] in file " + input_file)
            except HSTSpecPrevError as error_string:
                print error_string
                exit(1)
    elif list_len == 2:
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
            raise HSTSpecPrevError("The array of SEGMENT strings should contain 1, 2, or 3 values, found " + str(list_len) + " in file " + input_file)
        except HSTSpecPrevError as error_string:
                print error_string
                exit(1)
    return this_band

def generate_cos_avoid_regions():
    """
    Creates a list of COSAvoidRegion objects for use in the plotting routine, specifically designed for COS spectra.
    """
    lya1215_ar = COSAvoidRegion(1214.,1217., "Lyman alpha emission line.")
    return [lya1215_ar]

def plotspec(cos_spectrum, output_type, output_file, output_size=None):
    """
    Accepts a COSSpectrum object from the READSPEC function and produces preview plots.
    :param cos_spectrum: COS spectrum as returned by READSPEC.
    :type cos_spectrum: COSSpectrum
    :param output_type: What kind of output to make?
    :type output_type: str
    :param output_file: Name of output file (including full path).
    :type output_file: str
    :param output_size: Size of plot in pixels (plots are square in dimensions).  Defaults to 1024.
    :param output_size: int
    :raises: OSError,HSTSpecPrevError
    .. note::
       This function assumes a screen resolution of 96 DPI in order to generate plots of the desired sizes.  This is because matplotlib works in units of inches and DPI rather than pixels.
    """
    dpi_val = 96.
    if output_size is not None:
        if not isinstance(output_size, int):
            output_size = int(round(output_size))
    else:
        output_size = 1024
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
    """Extract the band type (FUV, NUV)."""
    this_band = cos_spectrum.band
    segment_names = cos_spectrum.segments.keys()
    """Reverse the list of segment names FOR FUV DATA, becaue the bluest segment is latter in the alphabet. but only for the FUV spectra.  If NUV, then just make sure the segments are sorted alphabetically."""
    if this_band == 'FUV':
        segment_names.sort(reverse=True)
    else:
        segment_names.sort(reverse=False)
    n_segments = len(segment_names)

    if output_size > 128:
        is_bigplot = True
        n_subplots = n_segments
        subplot_segment_names = segment_names
    else:
        is_bigplot = False
        """Change segment_names to be a one-element array, because we need to join all the segments into a single array if making a thumbnail plot.  Make sure it's still a list..."""
        n_subplots = 1
        subplot_segment_names = ["-".join(segment_names)]

    """Start plot figure."""
    this_figure, these_plotareas = pyplot.subplots(nrows=n_subplots, ncols=1,figsize=(output_size/dpi_val,output_size/dpi_val),dpi=dpi_val)

    if not isinstance(these_plotareas, numpy.ndarray):
        these_plotareas = numpy.asarray([these_plotareas])
    if is_bigplot:
        this_figure.subplots_adjust(hspace=0.3,top=0.915)
        this_figure.suptitle(os.path.basename(cos_spectrum.orig_file))
    else:
        this_figure.subplots_adjust(top=0.85,bottom=0.3,left=0.25,right=0.9)

    for i,s in enumerate(subplot_segment_names):
        this_plotarea = these_plotareas[i]
        if is_bigplot:
            all_wls = cos_spectrum.segments[s].wavelengths
            all_fls = cos_spectrum.segments[s].fluxes
            all_dqs = cos_spectrum.segments[s].dqs
        else:
            all_wls, all_fls, all_dqs, title_addendum = stitch_segments(cos_spectrum, segment_names)
        """Determine optimal x-axis."""
        x_axis_range = set_plot_xrange(all_wls, all_fls)
        """Plot the spectrum, but only if valid wavelength ranges for x-axis are returned, otherwise plot a special "Fluxes Are All Zero" plot."""
        if all(numpy.isfinite(x_axis_range)):
            """Create COS avoid regions."""
            avoid_regions = generate_cos_avoid_regions()
            """Determine optimal y-axis, but only provide it with fluxes from the part of the spectrum that will be plotted based on the x-axis trimming."""
            y_axis_range = set_plot_yrange(all_wls, all_fls,avoid_regions=avoid_regions,wl_range=x_axis_range)
            this_plotarea.plot(all_wls, all_fls, 'k')
            """Overplot the x-axis edges that are trimmed to define the y-axis plot range as a shaded area."""
            this_plotarea.axvspan(numpy.nanmin(all_wls), x_axis_range[0],facecolor="lightgrey",alpha=0.5)
            this_plotarea.axvspan(x_axis_range[1], numpy.nanmax(all_wls),facecolor="lightgrey",alpha=0.5)
            """Overplot the avoid regions in a light color as a shaded area."""
            for ar in avoid_regions:
                this_plotarea.axvspan(ar.minwl,ar.maxwl,facecolor="lightgrey",alpha=0.5)
            """Update the x-axis and y-axis range, but only adjust the ranges if this isn't an all-zero flux case."""
            """Note that we change the x-axis range here to be the min. and max. wavelength of this segment, rather than using the truncated version, so that all the plots for a similar instrument setting will have the same starting and ending plot values.  But, we still calculate the trimmed starting and ending wavelengths above for other things, such as defining the y-plot range."""
            x_axis_range = [numpy.nanmin(all_wls),numpy.nanmax(all_wls)]
            this_plotarea.set_xlim(x_axis_range)
            """Change x-axis units to microns if a small plot, because there isn't enough space."""
            if not is_bigplot:
                rc('font', size=10)
                this_plotarea.set_xticklabels(this_plotarea.get_xticks()/10000.,rotation=45.)
                this_plotarea.locator_params(axis="both", nbins=4, steps=[1,2,4,6,8,10])
                this_figure.suptitle(r'$\lambda (\mu$m)', position=(0.83,0.99))
            else:
                """Make sure the font properties go back to normal."""
                pyplot.rcdefaults()
                this_plotarea.set_xlabel("Angstroms")
            this_plotarea.set_ylim(y_axis_range)
        else:
            x_axis_range = [numpy.nanmin(all_wls),numpy.nanmax(all_wls)]
            this_plotarea.set_xlim(x_axis_range)
            this_plotarea.set_axis_bgcolor("lightgrey")
            if not is_bigplot:
                rc('font', size=10)
                this_plotarea.set_xticklabels(this_plotarea.get_xticks()/10000.,rotation=45.)
                this_plotarea.locator_params(axis="both", nbins=4, steps=[1,2,4,6,8,10])
                this_figure.suptitle(r'$\lambda in {\mu}m$', position=(0.83,0.99))
                textsize = "small"
            else:
                """Make sure the font properties go back to normal."""
                pyplot.rcdefaults()
                this_plotarea.set_xlabel("Angstroms")
                textsize = "x-large"
            this_plotarea.text(0.5,0.5,"Fluxes are all 0.",horizontalalignment="center",verticalalignment="center",transform=this_plotarea.transAxes,size=textsize)

    """Display or plot to the desired format."""
    if output_type != "screen":
        """Deconstruct output file to include plot size information."""
        output_splits = os.path.split(output_file)
        file_splits = os.path.splitext(output_splits[1])
        revised_output_file = output_splits[0]+os.path.sep+file_splits[0]+'_{0:04d}'.format(output_size)+file_splits[1]
        this_figure.savefig(revised_output_file, format=output_type, dpi=dpi_val)
    elif output_type == "screen":
        pyplot.show()

def readspec(input_file):
    """
    Reads in a COS spectrum FITS file (x1d, x1dsum, or x1dsum{1,2,3,4}) and returns the wavelengths, fluxes, and flux uncertainties for the two (FUV segments) or three (NUV stripes).
    :param input_file: Name of input FITS file.
    :type input_file: str
    :returns: COSSpectrum -- The spectroscopic data (wavelength, flux, flux error, etc):
    :raises: KeyError
    """
    with fits.open(input_file) as hdulist:
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
        try:
            dq_table = cos_tabledata.field("DQ")
        except KeyError:
            print "*** MAKE_HST_SPEC_PREVIEWS ERROR: DQ column not found in first extension's binary table."
            exit(1)
        """Create COSSegment objects to populate the COSSpectrum object with"""
        if band == 'FUV':
            try:
                fuva_index = numpy.where(segment_arr == "FUVA")[0][0]
            except IndexError:
                """Then there is no FUVA segment, normally this happends if there is only one segment present (in which case it's the other one of the two)."""
                fuva_index = None
            if fuva_index is not None:
                fuva_cossegment = COSSegment(nelem=nelems_arr[fuva_index], wavelengths=wavelength_table[fuva_index,:], fluxes=flux_table[fuva_index,:], fluxerrs=fluxerr_table[fuva_index,:], dqs=dq_table[fuva_index,:])
            try:
                fuvb_index = numpy.where(segment_arr == "FUVB")[0][0]
            except IndexError:
                """Then there is no FUVB segment, normally this happends if there is only one segment present (in which case it's the other one of the two)."""
                fuvb_index = None
            if fuvb_index is not None:
                fuvb_cossegment = COSSegment(nelem=nelems_arr[fuvb_index], wavelengths=wavelength_table[fuvb_index,:], fluxes=flux_table[fuvb_index,:], fluxerrs=fluxerr_table[fuvb_index,:], dqs=dq_table[fuvb_index,:])
        elif band == 'NUV':
            nuva_index = numpy.where(segment_arr == "NUVA")[0][0]
            nuva_cossegment = COSSegment(nelem=nelems_arr[nuva_index], wavelengths=wavelength_table[nuva_index,:], fluxes=flux_table[nuva_index,:], fluxerrs=fluxerr_table[nuva_index,:], dqs=dq_table[nuva_index,:])
            nuvb_index = numpy.where(segment_arr == "NUVB")[0][0]
            nuvb_cossegment = COSSegment(nelem=nelems_arr[nuvb_index], wavelengths=wavelength_table[nuvb_index,:], fluxes=flux_table[nuvb_index,:], fluxerrs=fluxerr_table[nuvb_index,:], dqs=dq_table[nuvb_index,:])
            nuvc_index = numpy.where(segment_arr == "NUVC")[0][0]
            nuvc_cossegment = COSSegment(nelem=nelems_arr[nuvc_index], wavelengths=wavelength_table[nuvc_index,:], fluxes=flux_table[nuvc_index,:], fluxerrs=fluxerr_table[nuvc_index,:], dqs=dq_table[nuvc_index,:])
        """Create COSSpectrum object."""
        if band == 'FUV':
            if fuva_index is not None and fuvb_index is not None:
                return_spec = COSSpectrum(band=band, cos_segments={'FUVA':fuva_cossegment,'FUVB':fuvb_cossegment}, orig_file=input_file)
            elif fuva_index is not None:
                return_spec = COSSpectrum(band=band, cos_segments={'FUVA':fuva_cossegment}, orig_file=input_file)
            else:
                return_spec = COSSpectrum(band=band, cos_segments={'FUVB':fuvb_cossegment}, orig_file=input_file)
        elif band == 'NUV':
            return_spec = COSSpectrum(band=band, cos_segments={'NUVA':nuva_cossegment,'NUVB':nuvb_cossegment,'NUVC':nuvc_cossegment}, orig_file=input_file)
        return return_spec

def _rms(values, offset=0.):
    """
    Calculates the RMS about some offset (default offset is 0.)
    :param values: Array of values to compute the rms of.
    :type values: numpy.ndarray
    :param offset: Optional offset to compute the rms about.  Defaults to 0.0.
    :type offset: float
    :returns: float -- A scalar float containing the rms about the offset.
    """
    return math.sqrt(numpy.nanmean([(x-offset)**2 for x in values]))

def _set_plot_xrange_test(flux_values, median_flux, lt0_fract):
    """
    Defines the test for an invalid part of the spectrum when trimming from the edges along the wavelength (x) axis.
    :param flux_values: An array of fluxes to test.
    :type flux_values: numpy.ndarray
    :param median_flux: A median flux value used in the test.
    :type median_flux: float
    :param lt0_fract: Fraction of fluxes that are less than zero out of the (good) part of the spectrum (0. <= lt0_fract <= 1.).
    :type lt0_fract: float
    :returns: list -- A list of True/False values depening on whether the input flux values pass the test.  Note that if a return value is True, then the flux value is considered PART OF THE SPECTRUM TO TRIM/SKIP OVER.  If median_flux is input as NaN, then this function returns True for all flux_values (i.e., skip all of them since median_flux is not defined).
    """
    if not isinstance(flux_values, numpy.ndarray):
        raise ValueError("The flux values must be passed as a numpy.ndarray object.")
    if numpy.isfinite(median_flux):
        if median_flux > 0. and lt0_fract < 0.1:
            return_var = [x <= 0. or median_flux/x >= 5. for x in flux_values]
        else:
            return_var = [x == 0. or median_flux/x >= 5. for x in flux_values]
    else:
        return_var = [True] * len(flux_values)
    return return_var
    
def set_plot_xrange(wavelengths,fluxes):
    """
    Given an array of wavelengths and fluxes, returns a list of [xmin,xmax] to define an optimal x-axis plot range.
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
    """Find the median flux value, ignoring any NaN values or fluxes that are 0.0."""
    where_finite_and_notzero = numpy.where( (numpy.isfinite(sorted_fluxes)) & (sorted_fluxes != 0.0) )
    n_finite_and_notzero = len(where_finite_and_notzero[0])
    if n_finite_and_notzero > 0:
        median_flux = numpy.median(sorted_fluxes[where_finite_and_notzero])
        n_lessthanzero = len(numpy.where(fluxes[where_finite_and_notzero] < 0.)[0])
        lessthanzero_fract = float(n_lessthanzero) / n_finite_and_notzero
    else:
        median_flux = numpy.nan
        lessthanzero_fract = numpy.nan
    """Find the first element in the array that is NOT considered an "edge effect", and the last element in the array that is NOT considered an "edge effect".  If the input array is all zeroes, then it will find the last and first element, respectively.  Note that the trim does not just stop at the first index that satisfies this requirement, since there can be spikes at the edges that can fool the algorithm.  Instead, it requires that the next "n_consecutive" data points after each trial location also fail the test for "edge effect"."""
    n_consecutive = 20
    start_index = 0
    end_index = -1
    n_fluxes = len(sorted_fluxes)
    done_trimming = False
    while not done_trimming:
        if start_index > n_fluxes-n_consecutive-1:
            done_trimming = True
        elif not numpy.any(_set_plot_xrange_test(sorted_fluxes[start_index:start_index+n_consecutive+1], median_flux, lessthanzero_fract)):
            """Test if next "n_consecutive" points also *fail( the edge effect test, e.g., they are from the *good* part of the spectrum, and if so, then we have found a good location and can break out of the while loop."""
            done_trimming = True
        else:
            start_index += 1
    done_trimming = False
    while not done_trimming:
        if end_index < -1*(n_fluxes-n_consecutive):
            done_trimming = True
        elif end_index != -1 and not numpy.any(_set_plot_xrange_test(sorted_fluxes[end_index-n_consecutive:end_index+1], median_flux, lessthanzero_fract)):
            """Test if next "n_consecutive" points also *fail( the edge effect test, e.g., they are from the *good* part of the spectrum, and if so, then we have found a good location and can break out of the while loop."""
            done_trimming = True
        elif end_index == -1 and not numpy.any(_set_plot_xrange_test(sorted_fluxes[end_index-n_consecutive:], median_flux, lessthanzero_fract)):
            """Also test if next "n_consecutive" points also *fail( the edge effect test, e.g., they are from the *good* part of the spectrum, and if so, then we have found a good location and can break out of the while loop.  This extra test is needed due to the vagaries of how python slicing syntax works with negaive indexes.  Probably could just re-write this entirely to use non-negative indexes, but the logic works either way."""
            done_trimming = True
        else:
            end_index -= 1
    """Return the optimal start and end wavelength values for defining the x-axis plot range.  Note that if the fluxes are all zeroes, then start index will be past end index, so we return NaN values to indicate a special plot should be made in that case.  The odd conditional below checks to make sure the end index (working from the back of the list via negative indexes) stops before reaching the start index (which works from the front using zero-based, positive indexes), otherwise return NaN values because the array is all zeroes."""
    if n_fluxes + end_index > start_index:
        return [sorted_wavelengths[start_index],sorted_wavelengths[end_index]]
    else:
        return [numpy.nan,numpy.nan]

def set_plot_yrange(wavelengths,fluxes,avoid_regions=None,wl_range=None):
    """
    Given an array of wavelengths, fluxes, and avoid regions, returns a list of [ymin,ymax] to define an optimal y-axis plot range.
    :param wavelengths: The wavelengths to be plotted.
    :type wavelengths: numpy.ndarray
    :param fluxes: The fluxes to be plotted.
    :type fluxes: numpy.ndarray
    :param avoid_regions: A list of wavelength ranges to avoid when calculating optimal y-axis plot range.
    :type avoid_regions: list of COSAvoidRegion objects.
    :param wl_range: The min. and max. wavelength that defines the x-axis plot range.  The default is None, in which case the min. and max. if the input wavelength array will be used.
    :type wl_range: list
    :returns: list -- Two-element list containing the optimal [ymin,ymax] values to define the y-axis plot range.
    .. note::
       This function makes use of an internal look-up table of wavelength regions where known contaminating emission lines or other strong UV artifacts can affect the zoom level of the plot.
    """
    if wl_range is None:
        wl_range = [numpy.nanmin(wavelengths), numpy.nanmax(wavelengths)]
    """This list will keep track of which fluxes to retain when defining the y-axis plot range, where setting the value to 1 means keep this flux for consideration."""
    keep_indices = [1] * len(wavelengths)
    if avoid_regions is not None:
        for i,ar in enumerate(avoid_regions):
            if i == 0:
                reject_indices = [i for i in range(len(wavelengths)) if wavelengths[i] >= ar.minwl and wavelengths[i] <= ar.maxwl or wavelengths[i] < wl_range[0] or wavelengths[i] > wl_range[1]]
            else:
                """Don't need to worry about checking wavelengths within bounds after the first avoid region is examined."""
                reject_indices = [i for i in range(len(wavelengths)) if wavelengths[i] >= ar.minwl and wavelengths[i] <= ar.maxwl]
            for j in reject_indices:
                keep_indices[j] = 0
    keep_fluxes = numpy.asarray([f for ii,f in enumerate(fluxes) if keep_indices[ii] == 1 and numpy.isfinite(fluxes[ii])])
    """Don't just take the pure min and max, since weird defects can affect the calculation.  Instead, take the 1th and 99th percentile fluxes within the region to consider."""
    min_flux = numpy.percentile(keep_fluxes,1.)
    max_flux = numpy.percentile(keep_fluxes,99.)
    """Determine a y-buffer based on the difference between the max. and min. flux."""
    ybuffer = 0.1 * (max_flux-min_flux)
    return [min_flux-ybuffer, max_flux+ybuffer]


def stitch_segments(input_spectrum, segment_names):
    """
    Given a COSSpectrum object, stitches each segment into a contiguous array.  Does not do any fitting or adjustments between the orders, but does trim out the edges of each order with flux exactly equal to 0.
    :param input_spectrum: The COS spectrum to stitch.
    :type input_spectrum: COSSpectrum
    :param segment_names: List of segment names.
    :type segment_names: list
    :returns: numpy array, numpy array, numpy array, str -- The stitched wavelengths, fluxes, flux errors, and an informational plot title in the event that all the fluxes were exactly equal to 0.
    """
    all_wls = []
    all_fls = []
    all_dqs = []
    n_segments = len(segment_names)
    all_dq_flags = numpy.zeros(n_segments)
    return_title = ""
    for j in segment_names:
        these_wls = input_spectrum.segments[j].wavelengths
        these_fls = input_spectrum.segments[j].fluxes
        these_dqs = input_spectrum.segments[j].dqs
        """Trim from the edges anything with a DQ flag > 0."""
        start_index = 0
        while start_index < len(these_dqs):
            if these_dqs[start_index] == 0:
                break
            else:
                start_index += 1
        end_index = -1
        while end_index >= -1*len(these_dqs):
            if these_dqs[end_index] == 0:
                break
            else:
                end_index -= 1
        """Only append the parts of this order's spectrum that are not DQ > 0 at the edges."""
        if len(these_dqs) + end_index > start_index:
            if end_index == -1:
                all_wls += list(these_wls[start_index:])
                all_fls += list(these_fls[start_index:])
                all_dqs += list(these_dqs[start_index:])
            else:
                all_wls += list(these_wls[start_index:end_index+1])
                all_fls += list(these_fls[start_index:end_index+1])
                all_dqs += list(these_dqs[start_index:end_index+1])
        else:
            all_dq_flags[j] = 1
            """Then the trimming from the edges passed one another, and this entire order has DQ > 0.  In this case, we include the entire order (for now, we may want to change this in the future though)."""
            all_wls += list(these_wls)
            all_fls += list(these_fls)
            all_dqs += list(these_dqs)
    """If every single order had all DQ flags, then we print out the warning."""
    if sum(all_dq_flags) == n_segments:
        return_title = "Warning: All fluxes have DQ > 0."
    all_wls = numpy.asarray(all_wls)
    all_fls = numpy.asarray(all_fls)
    all_dqs = numpy.asarray(all_dqs)
    sorted_indexes = numpy.argsort(all_wls)
    all_wls = all_wls[sorted_indexes]
    all_fls = all_fls[sorted_indexes]
    all_dqs = all_dqs[sorted_indexes]
    return all_wls, all_fls, all_dqs, return_title
