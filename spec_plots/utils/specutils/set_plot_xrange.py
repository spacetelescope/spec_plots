__version__ = '1.31'

"""
.. module:: set_plot_xrange

   :synopsis: Determines the optimal x-axis plot range, in wavelength space.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

import numpy

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

    :param flux_scale_factor: Max. allowed ratio between the flux and a median flux value, used in edge trimming.

    :type flux_scale_factor: float

    :param median_fluxerr: A median flux uncertainty value, used in the test.

    :type median_fluxerr: float

    :param fluxerr_scale_factor: Max. allowed ratio between the flux uncertainty and a median flux uncertainty value, used in edge trimming.

    :type fluxerr_scale_factor: float

    :param fluxerr_95th: The flux uncertainty corresponding to the 95th percentile.

    :type fluxerr_95th: float

    :param dqs: The array of DQ flags to use in the test.  For COS, these are the DQ_WGT bits from the header.

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
                             and ( (instrument == "stis" and not is_bad_dq(instrument, z) and checkDQs) or (instrument == "cos" and not is_bad_dq(instrument, z) and checkDQs) or (not checkDQs) ) \
                             for x,y,z in zip(flux_values, flux_err_values, dqs) ]
    else:
        bool_results = [False] * len(flux_values)

    return bool_results
#--------------------

#--------------------
def set_plot_xrange(instrument, wavelengths, fluxes, fluxerrs, dqs, n_consecutive, flux_scale_factor, fluxerr_scale_factor, median_flux, median_fluxerr, fluxerr_95th):
    """
    Given an array of wavelengths and fluxes, returns a list of [xmin,xmax] to define an optimal x-axis plot range.

    :param instrument: The instrument that is being tested.

    :type instrument: str

    :param wavelengths: The wavelengths to be plotted.

    :type wavelengths: numpy.ndarray

    :param fluxes: The fluxes to be plotted.

    :type fluxes: numpy.ndarray

    :param fluxerrs: The uncertainties of the fluxes to be plotted.

    :type fluxerrs: numpy.ndarray

    :param dqs: The DQ flags of the spectrum to be plotted.  For COS, these are the DQ_WGT bits from the header.

    :type dqs: numpy.ndarray

    :param n_consecutive: How many consecutive points must pass the test for the index to count as the valid start/end of the spectrum?

    :type n_consecutive: int

    :param flux_scale_factor: Max. allowed ratio between the flux and a median flux value, used in edge trimming.

    :type flux_scale_factor: float

    :param fluxerr_scale_factor: Max. allowed ratio between the flux uncertainty and a median flux uncertainty value, used in edge trimming.

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
