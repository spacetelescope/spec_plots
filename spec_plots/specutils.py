__version__ = '1.30'

"""
.. module:: specutils

   :synopsis: Contains utility functions for reading and plotting HST spectra.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

import math
import numpy
import specutils_cos
import specutils_stis

#--------------------

class SpecUtilsError(Exception, object):
    """
    This class defines a generic Exception to use for errors raised in the SPECUTILS modules (specutils, specutils_cos, specutils_stis, etc.).  It simply prints the given value when raising the exception, e.g., 
    
    .. code-block:: python
    
         raise SpecUtilsError("Print this string") 
         SpecUtilsError: *** SPECUTILS ERROR: 'Print this string'
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
        return "*** SPECUTILS ERROR: "+repr(self.value)

#--------------------

class AvoidRegion(object):
    """
    Defines an avoid region, which is a section of wavelength space that should not be included when determining the optimal y-axis plot range.  The object consists of a starting wavelength, ending wavelength, and string description of what that region is.

    :raises: ValueError
    """
    def __init__(self, minwl=None, maxwl=None, description=""):

        if minwl is None:
            raise ValueError("Must specify a minimum wavelength for this avoid region.")

        if maxwl is None:
            raise ValueError("Must specify a maximum wavelength for this avoid region.")

        if minwl >= maxwl:
            raise ValueError("Minimum wavelength must be less than maximum wavelength for this avoid region.  Given min. wavelength = "+str(minwl)+" and max. wavelength = "+str(maxwl)+".")

        """ Assign the min. wl., max. wl., and description to the object. """
        self.minwl = minwl
        self.maxwl = maxwl
        self.description = description

#--------------------

def calc_plot_metrics(instrument, wls, fls, flerrs, dqs, n_consecutive, flux_scale_factor, fluxerr_scale_factor):
    """
    
    """

    """ Calculate statistics on the fluxes. """
    median_flux, median_fluxerr, fluxerr_95th = get_flux_stats(fls, flerrs)

    """ Determine optimal x-axis.  This is not the x-axis plot range used, but rather the area of the plot that is considered when scaling the y-axis. """
    optimal_xaxis_range = set_plot_xrange(instrument, wls, fls, flerrs, dqs, n_consecutive, flux_scale_factor, fluxerr_scale_factor, median_flux, median_fluxerr, fluxerr_95th)

    """ Create COS avoid regions. """
    avoid_regions = generate_avoid_regions(instrument)

    """ Determine the optimal y-axis. """
    if all(numpy.isfinite(optimal_xaxis_range)):
        y_axis_range = set_plot_yrange(wls, fls, avoid_regions=avoid_regions, wl_range=optimal_xaxis_range)
    else:
        y_axis_range = [numpy.nan, numpy.nan]

    """ Return the plot_metrics dict. """
    return {"median_flux":median_flux, "median_fluxerr":median_fluxerr, "fluxerr_95th":fluxerr_95th, "optimal_xaxis_range":optimal_xaxis_range, "avoid_regions":avoid_regions, "y_axis_range":y_axis_range}

#--------------------

def debug_oplot(this_plotarea, all_wls, all_fls, all_flerrs, all_dqs, median_flux, median_fluxerr, flux_scale_factor, fluxerr_scale_factor, fluxerr_95th, oplot_percentiles=False):
    """
    Creates plots of the spectra with color-coding and special annotation to identify which points were rejected by which tests.  Useful for debugging and understanding why a given plot had its plot axes defined the way it did.

    :param this_plotarea: The AxesSubplot object to plot on.

    :type this_plotarea: matplotlib.axes._subplots.AxesSubplot

    :param all_wls: Array of wavelengths.
    
    :type all_wls: numpy.ndarray

    :param all_fls: Array of fluxes.
    
    :type all_fls: numpy.ndarray

    :param all_flerrs: Array of flux uncertainties.
    
    :type all_flerrs: numpy.ndarray

    :param all_dqs: Array of data quality flags.
    
    :type all_dqs: numpy.ndarray

    :param median_flux: The median flux used in determining where the best part of the spectrum is.

    :type median_flux: float

    :param median_fluxerr: The median flux uncertainty used in determining where the best part of the spectrum is.

    :type median_fluxerr: float

    :param flux_scale_factor: Max. allowed ratio between the flux and a median flux value used in edge trimming.

    :type flux_scale_factor: float

    :param fluxerr_scale_factor: Max. allowed ratio between the flux uncertainty and a median flux uncertainty value used in edge trimming.

    :type fluxerr_scale_factor: float

    :param fluxerr_95th: The flux uncertainty corresponding to the 95th percentile.

    :type fluxerr_95th: float

    :param oplot_percentiles: Set this to True to overplot points where the flux uncertainties are greater than the 95th percentile.  Default = False.

    :type oplot_percentiles: bool
    """

    """ Plot all of the points in black. """
    this_plotarea.errorbar(all_wls, all_fls, yerr=all_flerrs, ecolor='c', color='k', label='Passed')

    """ Plot those fluxes that fail because they are much greater than the median in blue. """
    if numpy.isfinite(median_flux):
        where_fluxtoolarge = numpy.where( abs(all_fls/median_flux) > flux_scale_factor )
        if len(where_fluxtoolarge[0]) > 0:
            this_plotarea.plot(all_wls[where_fluxtoolarge], all_fls[where_fluxtoolarge], 'bo', label="Flux>>Median")

    """ Plot those fluxes that fail because they are exactly equal to 0 in green. """
    where_allzero = numpy.where(all_fls == 0.0)
    if len(where_allzero[0]) > 0:
        this_plotarea.plot(all_wls[where_allzero], all_fls[where_allzero], 'go', label="Flux=0")

    """ Plot those fluxes that fail because they have a data quality flag > 0 in red. """
    where_dqnotzero = numpy.where((all_dqs > 0) & (all_dqs != 16))
    if len(where_dqnotzero[0]) > 0:
        this_plotarea.plot(all_wls[where_dqnotzero], all_fls[where_dqnotzero], 'ro', label="DQ>0")

    """ Plot those fluxes that fail because their flux uncertainties are much greater than the median flux uncertainty in magenta. """
    if numpy.isfinite(median_fluxerr):
        where_bigerr = numpy.where(all_flerrs/median_fluxerr > fluxerr_scale_factor)
        if len(where_bigerr[0]) > 0:
            this_plotarea.plot(all_wls[where_bigerr], all_fls[where_bigerr], 'mo', label="FluxErr>>Median")

    """ Plot those fluxes that fail because their flux uncertainties are greater than the 95th percentile in yellow. """
    if oplot_percentiles:
        where_bigerrpercentile = numpy.where(all_flerrs > fluxerr_95th)
        if len(where_bigerrpercentile[0]) > 0:
            this_plotarea.plot(all_wls[where_bigerrpercentile], all_fls[where_bigerrpercentile], 'yo', label="FluxErr>95th %")

    this_plotarea.legend(loc="upper center", ncol=4)

#--------------------

def dq_has_flag(dq, flag_to_check):
    """
    Returns true/false if the given DQ value has a specific flag value set after unpacked into a 16-bit string.  For example:
    
    .. code-block:: python
    
         print dq_has_flag(48,16)
         True
         print dq_has_flag(40, 16)
         False

    :param dq: The DQ value.

    :type dq: int

    :param flag_to_check: The flag value to check if it's set to True.

    :type flag_to_check: int

    :returns: bool -- Returns True if `flag_to_check` is set to True inside `dq`.

    :raises: ValueError
    """

    """ Make sure `flag_to_check` is a power of 2. """
    if (flag_to_check & (flag_to_check-1)) == 0 and flag_to_check != 0:
        dq_16bit_str = "{0:b}".format(dq)
        flag_16bit_str = "{0:b}".format(flag_to_check)
        
        """ If the 16-bit string of the value to check is longer than 16-bit string version of the DQ value, then we know it can't be part of the DQ bitmask.  If not, then look for that bit to be set (by counting from the right). """
        if len(flag_16bit_str) <= len(dq_16bit_str) and dq_16bit_str[-1*len(flag_16bit_str)] == '1':
            return True
        else:
            return False
    else:
        raise ValueError("Flag to check must be a power of 2.  Asked to check whether flag " + str(flag_to_check) + " is set to True in bitmask value " + str(dq) + ", but " + str(flag_to_check) + " is not a power of 2.")

#--------------------

def edge_trim(instrument, fluxes, fluxerrs, dqs, n_consecutive, median_flux, flux_scale_factor, median_fluxerr, fluxerr_scale_factor, fluxerr_95th):
    """
    Returns start and end indexes (end indexes are negatively indexed) of the best part of the spectrum to use when defining the plot's wavelength ranges.  Returns two sets of start and end indexes: one set without taking into account DQ flags, and one set that does take into account DQ flags (since some spectra have all DQ flags set > 0).

    :param instrument: The instrument that is being tested.

    :type instrument: str

    :param fluxes: The fluxes to be plotted.

    :type fluxes: numpy.ndarray

    :param fluxerrs: The uncertainties of the fluxes to be plotted.

    :type fluxerrs: numpy.ndarray

    :param dqs: The DQ flags of the spectrum.

    :type dqs: numpy.ndarray

    :param n_consecutive: The consecutive number of fluxes that must belong to the "best" part of the spectrum for the start/end index to count.

    :type n_consecutive: int

    :param median_flux: The median flux, used in determining where the best part of the spectrum is.

    :type median_flux: float

    :param flux_scale_factor: Max. allowed ratio between the flux and a median flux value, used in edge trimming.  Default = 10.

    :type flux_scale_factor: float

    :param median_fluxerr: The median flux uncertainty, used in determining where the best part of the spectrum is.

    :type median_fluxerr: float

    :param fluxerr_scale_factor: Max. allowed ratio between the flux uncertainty and a median flux uncertainty value, used in edge trimming.  Default = 5.

    :type fluxerr_scale_factor: float

    :param fluxerr_95th: The flux uncertainty corresponding to the 95th percentile.

    :type fluxerr_95th: float

    :returns: int tuple -- Indexes that define the best part of the spectrum, in the order of (start_index_nodq, end_index_nodq, start_index_withdq, end_index_withdq).
    """

    n_fluxes = len(fluxes)

    """ Determine the indices where each flux point satisfies our test criteria, ignoring DQ values. """
    where_good_fluxes_nodq = numpy.where(numpy.asarray(_set_plot_xrange_test(instrument, fluxes, fluxerrs, median_flux, flux_scale_factor, median_fluxerr, fluxerr_scale_factor, fluxerr_95th, dqs, checkFluxes=True, checkFluxRatios=False, checkFluxErrRatios=True, checkFluxErrPercentile=False, checkDQs=False)))[0]
    n_good_fluxes_nodq = len(where_good_fluxes_nodq)

    """ Handle the trivial case where the number of valid points is less than the `n_consecutive` desired.  In this case, the start and end indices are set to the last and first indices, simulating what would happen if the entire array was traversed. """
    if n_good_fluxes_nodq < n_consecutive:
        start_index_nodq = 0
        end_index_nodq = -1
    else:
        """ The start index is then the min. index for which the `n_consecutive`'th index is exactly equal to `n_consecutive`-1. """
        start_indexes_nodq_with_good_diffs = numpy.asarray([x for i,x in enumerate(where_good_fluxes_nodq[0:n_good_fluxes_nodq-n_consecutive+1]) \
                                                                if where_good_fluxes_nodq[i+n_consecutive-1] - \
                                                                where_good_fluxes_nodq[i] == n_consecutive-1])

        if len(start_indexes_nodq_with_good_diffs) > 0:
            start_index_nodq = numpy.min(start_indexes_nodq_with_good_diffs)
            end_index_nodq = numpy.max(start_indexes_nodq_with_good_diffs) - n_fluxes + n_consecutive - 1
        else:
            start_index_nodq = 0
            end_index_nodq = -1

    """ Determine the indices where each flux point satisfies our test criteria, taking into account DQ values. """
    where_good_fluxes_withdq = numpy.where(numpy.asarray(_set_plot_xrange_test(instrument, fluxes, fluxerrs, median_flux, flux_scale_factor, median_fluxerr, fluxerr_scale_factor, fluxerr_95th, dqs, checkFluxes=True, checkFluxRatios=False, checkFluxErrRatios=True, checkFluxErrPercentile=False, checkDQs=True)))[0]
    n_good_fluxes_withdq = len(where_good_fluxes_withdq)

    """ Handle the trivial case where the number of valid points is less than the `n_consecutive` desired.  In this case, the start and end indices are set to the last and first indices, simulating what would happen if the entire array was traversed. """
    if n_good_fluxes_withdq < n_consecutive:
        start_index_withdq = 0
        end_index_withdq = -1
    else:
        """ The start index is then the min. index for which the `n_consecutive`'th index is exactly equal to `n_consecutive`-1. """
        start_indexes_withdq_with_good_diffs = numpy.asarray([x for i,x in enumerate(where_good_fluxes_withdq[0:n_good_fluxes_withdq-n_consecutive+1]) \
                                                                if where_good_fluxes_withdq[i+n_consecutive-1] - \
                                                                where_good_fluxes_withdq[i] == n_consecutive-1])
        if len(start_indexes_withdq_with_good_diffs) > 0:
            start_index_withdq = numpy.min(start_indexes_withdq_with_good_diffs)
            end_index_withdq = numpy.max(start_indexes_withdq_with_good_diffs) - n_fluxes + n_consecutive - 1
        else:
            start_index_withdq = 0
            end_index_withdq = -1

    return start_index_nodq, end_index_nodq, start_index_withdq, end_index_withdq

#--------------------

def generate_avoid_regions(instrument):
    """
    Creates a list of AvoidRegion objects for use in the plotting routine, specifically designed for the given instrument.

    :param instrument: The instrument to generate the Avoid Region for.

    :type instrument: str

    """

    if instrument == "cos":
        lya1215_ar = AvoidRegion(1214.,1217., "Lyman alpha emission line.")
        return [lya1215_ar]

    elif instrument == "stis":
        lya1215_ar = AvoidRegion(1214.,1217., "Lyman alpha emission line.")
        return [lya1215_ar]

    else:
        return []

#--------------------

def get_flux_stats(fluxes, fluxerrs):
    """
    Calculates median values of the fluxes and flux uncertainties, and the 95th percentile flux uncertainty.

    :param fluxes: The fluxes to be plotted.

    :type fluxes: numpy.ndarray

    :param fluxerrs: The uncertainties of the fluxes to be plotted.

    :type fluxerrs: numpy.ndarray
    """

    """ Find the median flux value, ignoring any NaN values or fluxes that are 0.0. """
    where_finite_and_notzero = numpy.where( (numpy.isfinite(fluxes)) & (fluxes != 0.0) )

    if len(where_finite_and_notzero[0]) > 0:
        median_flux = numpy.median(fluxes[where_finite_and_notzero])
        median_fluxerr = numpy.median(fluxerrs[where_finite_and_notzero])
    else:
        median_flux = numpy.nan
        median_fluxerr = numpy.nan

    """ Get the 95th percentile flux uncertainty value. """
    fluxerr_95th = numpy.percentile(fluxerrs, 95.)

    return median_flux, median_fluxerr, fluxerr_95th

#--------------------

def rms(values, offset=0.):
    """
    Calculates the RMS about some offset (default offset is 0.)

    :param values: Array of values to compute the rms of.

    :type values: numpy.ndarray

    :param offset: Optional offset to compute the rms about.  Defaults to 0.0.

    :type offset: float

    :returns: float -- A scalar float containing the rms about the offset.
    """
    return math.sqrt(numpy.nanmean([(x-offset)**2 for x in values]))

#--------------------

def _set_plot_xrange_test(instrument, flux_values, flux_err_values, median_flux, flux_scale_factor, median_fluxerr, fluxerr_scale_factor, fluxerr_95th, dqs, checkFluxes=False, checkFluxRatios=False, checkFluxErrRatios=False, checkFluxErrPercentile=False, checkDQs=False):
    """
    Defines the test for an invalid part of the spectrum when trimming from the edges along the wavelength (x) axis.

    :param instrument: The instrument that is being tested.

    :type instrument: str

    :param flux_values: An array of fluxes to test.

    :type flux_values: numpy.ndarray

    :param flux_err_values: An array of associated flux uncertainties.

    :type flux_err_values: numpy.ndarray

    :param median_flux: A median flux value, used in the test.

    :type median_flux: float

    :param flux_scale_factor: Max. allowed ratio between the flux and a median flux value, used in edge trimming.  Default = 10.

    :type flux_scale_factor: float

    :param median_fluxerr: A median flux uncertainty value, used in the test.

    :type median_fluxerr: float

    :param fluxerr_scale_factor: Max. allowed ratio between the flux uncertainty and a median flux uncertainty value, used in edge trimming.  Default = 5.

    :type fluxerr_scale_factor: float

    :param fluxerr_95th: The flux uncertainty corresponding to the 95th percentile.

    :type fluxerr_95th: float

    :param dqs: The array of DQ flags to use in the test.

    :type dqs: numpy.ndarray

    :param checkFluxes: Should the value of the fluxes be used to test edge trimming?

    :type checkFluxes: bool

    :param checkFluxRatios: Should the ratio of the fluxes to the median flux value be used to test edge trimming?

    :type checkFluxRatios: bool

    :param checkFluxErrRatios: Should the ratio of the flux uncertainties to the median flux uncertainty be used to test edge trimming?

    :type checkFluxErrRatios: bool

    :param checkFluxErrPercentile: Should the highest percentile flux uncertainties be used to test edge trimming?

    :type checkFluxErrPercentile: bool

    :param checkDQs: Should the highest percentile flux uncertainties be used to test edge trimming?

    :type checkDQs: bool

    :returns: list -- A list of True/False values depening on whether the input flux values pass the test.  Note that if a return value is True, then the flux value is considered PART OF THE SPECTRUM TO KEEP.  If median_flux is input as NaN, then this function returns False for all flux_values (i.e., skip all of them since median_flux is not defined).

    :raises: ValueError
    """
    
    """ Make sure input parameters are all numpy arrays. """
    if not isinstance(flux_values, numpy.ndarray) or not isinstance(flux_err_values, numpy.ndarray) or not isinstance(dqs, numpy.ndarray):
        raise ValueError("The flux, flux uncertainty, and DQ values must be passed as a numpy.ndarray object.")

    """ Return array of boolean values for the edge_trim test. """
    if numpy.isfinite(median_flux):
        bool_results = [ ( (instrument == "cos" and x != 0. and checkFluxes) or (instrument == "stis" and x != 0. and checkFluxes) or (not checkFluxes) ) \
                             and ( (abs(x/median_flux) < flux_scale_factor and checkFluxRatios) or (not checkFluxRatios) ) \
                             and ( (y/median_fluxerr < fluxerr_scale_factor and checkFluxErrRatios) or (not checkFluxErrRatios) ) \
                             and ( (y <= fluxerr_95th and checkFluxErrPercentile) or (not checkFluxErrPercentile) ) \
                             and ( (instrument == "stis" and z == 0 or z == 16 and checkDQs) or (instrument == "cos" and z == 0 and checkDQs) or (not checkDQs) ) \
                             for x,y,z in zip(flux_values, flux_err_values, dqs) ]
    else:
        bool_results = [False] * len(flux_values)

    return bool_results
    
#--------------------

def set_plot_xrange(instrument, wavelengths, fluxes, fluxerrs, dqs, n_consecutive, flux_scale_factor, fluxerr_scale_factor, median_flux, median_fluxerr, fluxerr_95th):
    """
    Given an array of wavelengths and fluxes, returns a list of [xmin,xmax] to define an optimal x-axis plot range.

    :param wavelengths: The wavelengths to be plotted.

    :param instrument: The instrument that is being tested.

    :type instrument: str

    :type wavelengths: numpy.ndarray

    :param fluxes: The fluxes to be plotted.

    :type fluxes: numpy.ndarray

    :param fluxerrs: The uncertainties of the fluxes to be plotted.

    :type fluxerrs: numpy.ndarray

    :param dqs: The DQ flags of the spectrum to be plotted.

    :type dqs: numpy.ndarray

    :param n_consecutive: How many consecutive points must pass the test for the index to count as the valid start/end of the spectrum?  Default = 20.

    :type n_consecutive: int

    :param flux_scale_factor: Max. allowed ratio between the flux and a median flux value, used in edge trimming.  Default = 10.

    :type flux_scale_factor: float

    :param fluxerr_scale_factor: Max. allowed ratio between the flux uncertainty and a median flux uncertainty value, used in edge trimming.  Default = 5.

    :type fluxerr_scale_factor: float

    :param median_flux: A median flux value, used in the test.

    :type median_flux: float

    :param median_fluxerr: A median flux uncertainty value, used in the test.

    :type median_fluxerr: float

    :param fluxerr_95th: The flux uncertainty corresponding to the 95th percentile.

    :type fluxerr_95th: float

    :returns: list -- Two-element list containing the optimal [xmin,xmax] values to define the x-axis plot range.
    """

    """ Test if there are any NaN's in the wavelength array.  If so, issue a warning for now... """
    if numpy.isnan(numpy.sum(wavelengths)):
        print "***WARNING in SPECUTILS_STIS: Wavelength array contains NaN values.  Behavior has not been fully tested in this case."

    """ Find the first element in the array that is NOT considered an "edge effect", and the last element in the array that is NOT considered an "edge effect".  If the input array is all zeroes, then it will find the last and first element, respectively.  Note that the trim does not just stop at the first index that satisfies this requirement, since there can be spikes at the edges that can fool the algorithm.  Instead, it requires that the next "n_consecutive" data points after each trial location also fail the test for "edge effect".  It returns start and stop indexes taking into account DQ values, as well as start and stop indexes without taking into account DQ values. """
    start_index_nodq, end_index_nodq, start_index_withdq, end_index_withdq = edge_trim(instrument, fluxes, fluxerrs, dqs, n_consecutive, median_flux, flux_scale_factor, median_fluxerr, fluxerr_scale_factor, fluxerr_95th)

    """ Return the optimal start and end wavelength values for defining the x-axis plot range.  Note that if the fluxes are all zeroes, then start index will be past end index, so we return NaN values to indicate a special plot should be made in that case.  The odd conditional below checks to make sure the end index (working from the back of the list via negative indexes) stops before reaching the start index (which works from the front using zero-based, positive indexes), otherwise return NaN values because the array is all zeroes. """
    if len(fluxes) + end_index_withdq > start_index_withdq:
        return [wavelengths[start_index_withdq], wavelengths[end_index_withdq]]
    elif len(fluxes) + end_index_nodq > start_index_nodq:
        return [wavelengths[start_index_nodq], wavelengths[end_index_nodq]]
    else:
        return [numpy.nan, numpy.nan]

#--------------------

def set_plot_yrange(wavelengths, fluxes, avoid_regions=None, wl_range=None):
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

    """ Default the wavelength range to be the entire spectrum if not specified. """
    if wl_range is None:
        wl_range = [numpy.nanmin(wavelengths), numpy.nanmax(wavelengths)]

    """ This list will keep track of which fluxes to retain when defining the y-axis plot range.  Setting the value to 1 means keep this flux for consideration. """
    keep_indices = numpy.asarray([1] * len(wavelengths))

    """ Set any fluxes to 0 (don't keep) if they fall within an Avoid Region. """
    if avoid_regions is not None:
        for i,ar in enumerate(avoid_regions):
            if i == 0:
                """ If this is the first Avoid Region, then we need to check if the wavelengths are within the specified bounds supplied through the `wl_range` parameter, in addition to checking if they are within the Avoid Region itself. """
                reject_indices = [i for i in range(len(wavelengths)) if wavelengths[i] >= ar.minwl and wavelengths[i] <= ar.maxwl or wavelengths[i] < wl_range[0] or wavelengths[i] > wl_range[1]]
            else:
                """ Don't need to worry about checking wavelengths within bounds after the first Avoid Region is examined. """
                reject_indices = [i for i in range(len(wavelengths)) if wavelengths[i] >= ar.minwl and wavelengths[i] <= ar.maxwl]
            """ Set indices that we don't want to keep to 0.  Note that if reject_indices is an empty list then nothing will change. """
            keep_indices[reject_indices] = 0

    """ After all indices have been set to keep or reject, pull out just the fluxes that should be kept. """
    keep_fluxes = numpy.asarray([f for ii,f in enumerate(fluxes) if keep_indices[ii] == 1 and numpy.isfinite(fluxes[ii])])

    """ Don't just take the pure min and max, since large outliers can affect the calculation.  Instead, take the 1th and 99th percentile fluxes within the region to calculate the `min` and `max` fluxes. """
    min_flux = numpy.percentile(keep_fluxes,1.)
    max_flux = numpy.percentile(keep_fluxes,99.)

    """ Determine a y-buffer based on the difference between the max. and min. fluxes. """
    ybuffer = 0.3 * (max_flux-min_flux)

    """ Make sure the min. and max. fluxes aren't identical (both 0., or both the same exact value.  If so, just return the min. and max.  value nudged by 1.0. """
    if min_flux != max_flux:
        return [min_flux-ybuffer, max_flux+ybuffer]
    else:
        return [min_flux-1., max_flux+1.]

#--------------------

def stitch_components(input_exposure, n_consecutive, flux_scale_factor, fluxerr_scale_factor, segment_names=None):
    """
    Given a COSSpectrum or STISExposureSpectrum object, stitches each segment/order, respectively, into a contiguous array.

    :param input_exposure: The COS segment or STIS exposure spectrum to stitch.

    :type input_exposure: COSSpectrum or STISExposureSpectrum

    :param n_consecutive: How many consecutive points must pass the test for the index to count as the valid start/end of the spectrum?  Default = 20.

    :type n_consecutive: int

    :param flux_scale_factor: Max. allowed ratio between the flux and a median flux value, used in edge trimming.  Default = 10.

    :type flux_scale_factor: float

    :param fluxerr_scale_factor: Max. allowed ratio between the flux uncertainty and a median flux uncertainty value, used in edge trimming.  Default = 5.

    :type fluxerr_scale_factor: float

    :param segment_names: List of segment names if input_exposure is a COSSpectrum object.

    :type segment_names: list

    :returns: dict{numpy array, numpy array, numpy array, str} -- The stitched wavelengths, fluxes, flux errors, and an informational plot title in the event that all the fluxes had the DQ flag set.  These are packaged in a dict.

    :raises: ValueError
    """

    """ These lits will contain the stitched spectrum. """
    all_wls = [] ; all_fls = [] ; all_flerrs = [] ; all_dqs = []
    
    """ Determine how many pieces there are to stitch, and how to loop through them (different depending on instrument type. """
    if isinstance(input_exposure, specutils_cos.COSSpectrum):
        if segment_names is not None:
            n_components = len(segment_names)
            loop_iterable = segment_names
            inst_type = "cos"
        else:
            raise ValueError("Must provide a list of segment names for COS spectra.")

    elif isinstance(input_exposure, specutils_stis.STISExposureSpectrum):
        n_components = len(input_exposure.orders)
        loop_iterable = xrange(n_components)
        inst_type = "stis"

    else:
        raise ValueError("Input must be either a COSSpectrum or STISExposureSpectrum object.")
    
    """ Create an array that will keep track whether each component to stitch is filled with bad DQ flags.  The return title will be set to a non-empty string if every component to stitch is filled with bad DQ flags. """
    all_dq_flags = numpy.zeros(n_components)
    return_title = ""

    """ Get the wavelengths, fluxes, flux uncertainties, and DQ flags for this component. """
    for jj,j in enumerate(loop_iterable):
        if inst_type == "stis":
            these_wls = input_exposure.orders[j].wavelengths
            these_fls = input_exposure.orders[j].fluxes
            these_flerrs = input_exposure.orders[j].fluxerrs
            these_dqs = input_exposure.orders[j].dqs
            where_bad_dq = numpy.where( (these_dqs > 0) & (these_dqs != 16))[0]
        elif inst_type == "cos":
            these_wls = input_exposure.segments[j].wavelengths
            these_fls = input_exposure.segments[j].fluxes
            these_flerrs = input_exposure.segments[j].fluxerrs
            these_dqs = input_exposure.segments[j].dqs
            where_bad_dq = numpy.where(these_dqs > 0)[0]
        n_elems = len(these_wls)

        """ Calculate some statistics for this component. """
        median_flux, median_fluxerr, fluxerr_95th = get_flux_stats(these_fls, these_flerrs)

        """ Find the start and end indices using the edge trimming function. """
        start_index_nodq, end_index_nodq, start_index_withdq, end_index_withdq = edge_trim(inst_type, these_fls, these_flerrs, these_dqs, n_consecutive, median_flux, flux_scale_factor, median_fluxerr, fluxerr_scale_factor, fluxerr_95th)

        """ Check if this component has all bad DQ flags. """
        if len(where_bad_dq) == n_elems:
            all_dq_flags[jj] = 1

        """ Only append the parts of this order's spectrum that pass the edge trimming, but do NOT trim too much. """
        if (start_index_withdq+1) / float(n_elems) <= 0.1:
            start_index_to_use = start_index_withdq
        elif (start_index_nodq+1) / float(n_elems) <= 0.1:
            start_index_to_use = start_index_nodq
        else:
            start_index_to_use = 0

        if abs(end_index_withdq) / float(n_elems) <= 0.1:
            end_index_to_use = end_index_withdq
        elif abs(end_index_nodq) / float(n_elems) <= 0.1:
            end_index_to_use = end_index_nodq
        else:
            end_index_to_use = -1

        """ Append the portion of this component to the spectrum we are building up. """
        if end_index_to_use == -1:
            all_wls += list(these_wls[start_index_to_use:])
            all_fls += list(these_fls[start_index_to_use:])
            all_flerrs += list(these_flerrs[start_index_to_use:])
            all_dqs += list(these_dqs[start_index_to_use:])
        else:
            all_wls += list(these_wls[start_index_to_use:end_index_to_use+1])
            all_fls += list(these_fls[start_index_to_use:end_index_to_use+1])
            all_flerrs += list(these_flerrs[start_index_to_use:end_index_to_use+1])
            all_dqs += list(these_dqs[start_index_to_use:end_index_to_use+1])

    """ If every single order had all DQ flags, then we want to print a warning on the plot. """
    if sum(all_dq_flags) == n_components:
        if inst_type == 'stis':
            return_title = "Warning: All fluxes have DQ > 0 and != 16."
        elif inst_type == 'cos':
            return_title = "Warning: All fluxes have DQ > 0."

    """ Convert these to numpy arrays. """
    """ <DEVEL> Should these just be numpy arrays to begin with? </DEVEL> """
    all_wls = numpy.asarray(all_wls)
    all_fls = numpy.asarray(all_fls)
    all_flerrs = numpy.asarray(all_flerrs)
    all_dqs = numpy.asarray(all_dqs)

    """ Make sure the spectrum is sorted in wavelength. """
    """ <DEVEL> Not sure how to handle components that might overlap in wavelength space, but are not monotonically increasing/decreasing as each component is stitched. </DEVEL> """
    sorted_indexes = numpy.argsort(all_wls)
    all_wls = all_wls[sorted_indexes]
    all_fls = all_fls[sorted_indexes]
    all_flerrs = all_flerrs[sorted_indexes]
    all_dqs = all_dqs[sorted_indexes]

    return {"wls":all_wls, "fls":all_fls, "flerrs":all_flerrs, "dqs":all_dqs, "title":return_title}

#--------------------
