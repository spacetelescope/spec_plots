__version__ = '1.0'

"""
.. module:: specutils_stis
   :synopsis: Contains functions for reading and plotting HST STIS spectra.
.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

from astropy.io import fits
import math
import matplotlib.pyplot as pyplot
import numpy
import os

class STISAvoidRegion:
    """
    Defines a STIS avoid region, which is simply a section of wavelength space that should not be included when determining the optimal y-axis plot range.  The object consists of a starting wavelength, ending wavelength, and string description of what that region is.
    :raises: ValueError
    """
    def __init__(self, minwl=None, maxwl=None, description=""):
        if minwl is None:
            raise ValueError("Must specify a minimum wavelength for this STIS avoid region.")
        if maxwl is None:
            raise ValueError("Must specify a maximum wavelength for this STIS avoid region.")
        if minwl >= maxwl:
            raise ValueError("Minimum wavelength must be less than maximum wavelength for this STIS avoid region.  Given min. wavelength = "+str(minwl)+" and max. wavelength = "+str(maxwl)+".")
        self.minwl = minwl
        self.maxwl = maxwl
        self.description = description

class STIS1DSpectrum:
    """
    Defines a STIS 1D spectrum (either "x1d" extracted or "sx1" summed extracted), including wavelegnth, flux, and flux errors.  A STIS 1D spectrum object consists of N associations.  If the file is an association, then N > 1, otherwise N = 1.  Each of these N associations can contain M orders.  If the association is an Echelle spectrum, then 24 < M < 70, depending on instrument configuration, otherwise M = 1.  Each of these M orders contain typical spectral data (wavelengths, fluxes, etc.), stored as STISOrderSpectrum objects.  The final data structure is then <STIS1DSpectrum>.associations[n].order[m].wavelengths (or .fluxes, .fluxerrs, etc.).
    """
    def __init__(self, association_spectra=None, orig_file=None):
        """
        Create a STIS1DSpectrum object out of a list of STISExposureSpectrum objects, which themselves are lists of STISOrderSpectrum objects.
        :param association_spectra: A list whose length is equal to the number of associations (length = "N" associations).  Each element in this list is a STISExposureSpectrum object, which itself is a list (length = "M" orders) of STISOrderSpectrum objects.
        :type association_spectra: list
        :param orig_file: Original FITS file read to create the spectrum (includes full path).
        :type orig_file: str
        :raises: ValueError
        """
        if association_spectra is not None:
            self.orig_file = orig_file
            self.associations = association_spectra
        else:
            raise ValueError("Must provide a list of STISExposureSpectrum objects.")

class STISExposureSpectrum:
    """
    Defines a STIS exposure spectrum, which consists of "M" STISOrderSpectrum objects.
    """
    def __init__(self, order_spectra=None):
        """
        Create a STISExposureSpectrum object out of a list of STISOrderSpectrum objects.
        :param order_spectra: The STISOrderSpectrum objects to build the STISExposureSpectrum object out of.
        :type order_spectra: list
        :raises: ValueError
        """
        if order_spectra is None:
            raise ValueError("Must provide a list of STISOrderSpectrum objects.")
        elif len(order_spectra) > 0:
            self.orders = order_spectra
        else:
            raise ValueError("Must provide a list of at least one STISOrderSpectrum objects, input list is empty.")

class STISOrderSpectrum:
    """
    Defines a STIS order spectrum, including wavelength, flux, and flux errors, which are stored as numpy arrays.  A scalar int property provides the number of elements in this segment.
    """
    def __init__(self, nelem=None, wavelengths=None, fluxes=None, fluxerrs=None):
        """
        Create a STISOrderSpectrum object, default to empty values.  Allows user to preallocate space if they desire by setting "nelem" but not providing lists/arrays on input right away.
        :param nelem: Number of elements for this segment's spectrum.
        :type nelem: int
        :param wavelengths: List of wavelengths in this segment's spectrum.
        :type wavelengths: list
        :param fluxes: List of fluxes in this segment's spectrum.
        :type fluxes: list
        :param fluxerrs: List of flux uncertainties in this segment's spectrum.
        :type fluxerrs: list
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

def generate_stis_avoid_regions():
    """
    Creates a list of STISAvoidRegion objects for use in the plotting routine, specifically designed for STIS spectra.
    """
    lya1215_ar = STISAvoidRegion(1214.,1217., "Lyman alpha emission line.")
    return [lya1215_ar]

def plotspec(stis_spectrum, output_type, output_file):
    """
    Accepts a STIS1DSpectrum object from the READSPEC function and produces preview plots.
    :param stis_spectrum: STIS spectrum as returned by READSPEC.
    :type stis_spectrum: STIS1DSpectrum
    :param output_type: What kind of output to make?
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
    """Start plot figure."""
    this_figure = pyplot.figure(figsize=(800./dpi_val,600./dpi_val),dpi=dpi_val)
    this_figure.subplots_adjust(hspace=0.3,top=0.935)
    this_figure.suptitle(os.path.basename(stis_spectrum.orig_file))
    n_associations = len(stis_spectrum.associations)
    for i in xrange(n_associations):
        this_plotarea = pyplot.subplot(n_associations,1,i+1)
        for j in xrange(len(stis_spectrum.associations[i].orders)):
            wls = stis_spectrum.associations[i].orders[j].wavelengths
            fls = stis_spectrum.associations[i].orders[j].fluxes
            """Determine optimal x-axis."""
            x_axis_range = set_plot_xrange(wls, fls)
            if all(numpy.isfinite(x_axis_range)):
                """Create COS avoid regions."""
                avoid_regions = generate_stis_avoid_regions()
                """Determine optimal y-axis, but only provide it with fluxes from the part of the spectrum that will be plotted based on the x-axis trimming."""
                y_axis_range = set_plot_yrange(wls,fls,avoid_regions=avoid_regions,wl_range=x_axis_range)
                pyplot.plot(wls, fls, 'k')
                """Overplot the x-axis edges that are trimmed to define the y-axis plot range as a shaded area."""
                pyplot.axvspan(numpy.nanmin(wls), x_axis_range[0],facecolor="lightgrey",alpha=0.5)
                pyplot.axvspan(x_axis_range[1], numpy.nanmax(wls),facecolor="lightgrey",alpha=0.5)
                """Overplot the avoid regions in a light color as a shaded area."""
                for ar in avoid_regions:
                    pyplot.axvspan(ar.minwl,ar.maxwl,facecolor="lightgrey",alpha=0.5)
                """Update the x-axis and y-axis range, but only adjust the ranges if this isn't an all-zero flux case."""
                """Note that we change the x-axis range here to be the min. and max. wavelength of this segment, rather than using the truncated version, so that all the plots for a similar instrument setting will have the same starting and ending plot values.  But, we still calculate the trimmed starting and ending wavelengths above for other things, such as defining the y-plot range."""
                x_axis_range = [numpy.nanmin(wls),numpy.nanmax(wls)]
                this_plotarea.set_xlim(x_axis_range)
                this_plotarea.set_ylim(y_axis_range)
            else:
                this_plotarea.set_axis_bgcolor("lightgrey")
                this_plotarea.text(0.5,0.5,"Fluxes are all 0.",horizontalalignment="center",verticalalignment="center",transform=this_plotarea.transAxes,size="x-large")

    """Display or plot to the desired format."""
    if output_type != "screen":
        pyplot.savefig(output_file, format=output_type, dpi=dpi_val,bbox_inches='tight')
    elif output_type == "screen":
        pyplot.show()

def readspec(input_file):
    """
    Reads in a STIS spectrum FITS file (x1d, sx1) and returns the wavelengths, fluxes, and flux uncertainties for the two (FUV segments) or three (NUV stripes).
    :param input_file: Name of input FITS file.
    :type input_file: str
    :returns: STIS1DSpectrum -- The spectroscopic data (wavelength, flux, flux error, etc):
    """
    with fits.open(input_file) as hdulist:
        """ Read in the number of extensions from the primary header.  This will determine whether this is an association (N > 1) or not (N = 1)."""
        n_associations = hdulist[0].header["NEXTEND"]
        """Create an initially empty list that will contain each extension's (association's) spectrum object."""
        all_association_spectra = []
        for exten in xrange(n_associations):
            exten_data_table = hdulist[exten+1].data
            """How many orders (table rows) in this extension?"""
            n_orders = len(exten_data_table["sporder"])
            """Create a list of STISOrderSpectra for this extension."""
            all_order_spectra = [STISOrderSpectrum(nelem=exten_data_table["nelem"][order], wavelengths=exten_data_table["WAVELENGTH"][order], fluxes=exten_data_table["FLUX"][order], fluxerrs=exten_data_table["ERROR"][order]) for order in xrange(n_orders)]
            """Create a STISExposureSpectrum from the STISOrderSpectrum objects."""
            this_exposure_spectrum = STISExposureSpectrum(order_spectra=all_order_spectra)
            all_association_spectra.append(this_exposure_spectrum)
        return STIS1DSpectrum(association_spectra=all_association_spectra, orig_file=input_file)

def _set_plot_xrange_test(flux_values, median_flux):
    """
    Defines the test for an invalid part of the spectrum when trimming from the edges along the wavelength (x) axis.
    :param flux_values: A scalar float or list of fluxes to test.
    :type flux_values: float or list
    :param median_flux: A median flux value used in the test.
    :type median_flux: float
    :returns: float or list -- A scalar float or list of True/False values depening on whether the input flux values pass the test.  Return type matches the type of the input flux values.  Note that if a return value is True, then the flux value is considered PART OF THE SPECTRUM TO TRIM/SKIP OVER.
    """
    try:
        return_var = [x <= 0. or median_flux/x >= 5. for x in flux_values]
    except TypeError:
        return_var = flux_values <= 0. or median_flux/flux_values >= 5.
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
        print "***WARNING in SPECUTILS_STIS: Wavelength array contains NaN values.  Behavior has not been fully tested in this case."
    """Sort the wavelength and flux arrays.  Don't do it in-place since we aren't (yet) going to sort the fluxes, errors, etc. as well."""
    sorted_indexes = numpy.argsort(wavelengths)
    sorted_wavelengths = wavelengths[sorted_indexes]
    sorted_fluxes = fluxes[sorted_indexes]
    """Find the median flux value, ignoring any NaN values or fluxes that are 0.0."""
    where_finite_and_notzero = numpy.where( (numpy.isfinite(sorted_fluxes)) & (sorted_fluxes != 0.0) )
    median_flux = numpy.median(sorted_fluxes[where_finite_and_notzero])
    """Find the first element in the array that is NOT considered an "edge effect", and the last element in the array that is NOT considered an "edge effect".  If the input array is all zeroes, then it will find the last and first element, respectively.  Note that the trim does not just stop at the first index that satisfies this requirement, since there can be spikes at the edges that can fool the algorithm.  Instead, it requires that the next "n_consecutive" data points after each trial location also fail the test for "edge effect"."""
    n_consecutive = 20
    start_index = 0
    end_index = -1
    n_fluxes = len(sorted_fluxes)
    done_trimming = False
    while not done_trimming:
        if start_index > n_fluxes-n_consecutive-1:
            done_trimming = True
        elif not numpy.any(_set_plot_xrange_test(sorted_fluxes[start_index:start_index+n_consecutive+1], median_flux)):
            """Test if next "n_consecutive" points also *fail( the edge effect test, e.g., they are from the *good* part of the spectrum, and if so, then we have found a good location and can break out of the while loop."""
            done_trimming = True
        else:
            start_index += 1
    done_trimming = False
    while not done_trimming:
        if end_index < -1*(n_fluxes-n_consecutive):
            done_trimming = True
        elif end_index != -1 and not numpy.any(_set_plot_xrange_test(sorted_fluxes[end_index-n_consecutive:end_index+1], median_flux)):
            """Test if next "n_consecutive" points also *fail( the edge effect test, e.g., they are from the *good* part of the spectrum, and if so, then we have found a good location and can break out of the while loop."""
            done_trimming = True
        elif end_index == -1 and not numpy.any(_set_plot_xrange_test(sorted_fluxes[end_index-n_consecutive:], median_flux)):
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
    :type avoid_regions: list of STISAvoidRegion objects.
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
