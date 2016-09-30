"""
.. module:: edge_trim
   :synopsis: Returns the start and end indices of a spectrum that represent the
       `best` wavelengths to use when defining the optimal plot axis ranges.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

#--------------------
# Built-In Imports
#--------------------
from __future__ import absolute_import
from __future__ import division
from builtins import zip
#--------------------
# External Imports
#--------------------
import numpy
#--------------------
# Package Imports
#--------------------
from spec_plots.utils.specutils.is_bad_dq import is_bad_dq
from spec_plots import __version__

#--------------------

def _set_plot_xrange_test(instrument, flux_values, flux_err_values, median_flux,
                          flux_scale_factor, median_fluxerr,
                          fluxerr_scale_factor, fluxerr_95th, dqs,
                          check_fluxes=False, check_flux_ratios=False,
                          check_flux_err_ratios=False,
                          check_flux_err_percentile=False, check_dqs=False):
    """
    Defines the test for an invalid part of the spectrum when trimming from the
    edges along the wavelength (x) axis.

    :param instrument: The instrument that is being tested.

    :type instrument: str

    :param flux_values: An array of fluxes to test.

    :type flux_values: numpy.ndarray

    :param flux_err_values: An array of associated flux uncertainties.

    :type flux_err_values: numpy.ndarray

    :param median_flux: A median flux value, used in the test.

    :type median_flux: float

    :param flux_scale_factor: Max. allowed ratio between the flux and a median
        flux value, used in edge trimming.

    :type flux_scale_factor: float

    :param median_fluxerr: A median flux uncertainty value, used in the test.

    :type median_fluxerr: float

    :param fluxerr_scale_factor: Max. allowed ratio between the flux uncertainty
        and a median flux uncertainty value, used in edge trimming.

    :type fluxerr_scale_factor: float

    :param fluxerr_95th: The flux uncertainty corresponding to the 95th
        percentile.

    :type fluxerr_95th: float

    :param dqs: The array of DQ flags to use in the test.  For COS, these are
        the DQ_WGT bits from the header.

    :type dqs: numpy.ndarray

    :param check_fluxes: Should the value of the fluxes be used to test edge
        trimming?

    :type check_fluxes: bool

    :param check_flux_ratios: Should the ratio of the fluxes to the median flux
        value be used to test edge trimming?

    :type check_flux_ratios: bool

    :param check_flux_err_ratios: Should the ratio of the flux uncertainties to
        the median flux uncertainty be used to test edge trimming?

    :type check_flux_err_ratios: bool

    :param check_flux_err_percentile: Should the highest percentile flux
        uncertainties be used to test edge trimming?

    :type check_flux_err_percentile: bool

    :param check_dqs: Should the highest percentile flux uncertainties be used
        to test edge trimming?

    :type check_dqs: bool

    :returns: list -- A list of True/False values depening on whether the input
    flux values pass the test.  Note that if a return value is True, then the
    flux value is considered PART OF THE SPECTRUM TO KEEP.  If median_flux is
    input as NaN, then this function returns False for all flux_values (i.e.,
    skip all of them since median_flux is not defined).

    :raises: ValueError
    """

    # Make sure input parameters are all numpy arrays.
    if (not isinstance(flux_values, numpy.ndarray) or
            not isinstance(flux_err_values, numpy.ndarray) or
            not isinstance(dqs, numpy.ndarray)):
        raise ValueError("The flux, flux uncertainty, and DQ values must be"
                         " passed as a numpy.ndarray object.")

    # Return array of boolean values for the edge_trim test.
    if numpy.isfinite(median_flux):
        bool_results = [((instrument == "cos" and x != 0. and check_fluxes) or
                         (instrument == "stis" and x != 0. and check_fluxes) or
                         (instrument == "miri" and x != 0. and check_fluxes) or
                         (not check_fluxes))
                        and ((abs(x/median_flux) < flux_scale_factor and
                              check_flux_ratios) or (not check_flux_ratios))
                        and ((y/median_fluxerr < fluxerr_scale_factor and
                              check_flux_err_ratios) or
                             (not check_flux_err_ratios))
                        and ((y <= fluxerr_95th and check_flux_err_percentile)
                             or (not check_flux_err_percentile))
                        and ((instrument == "stis" and not
                              is_bad_dq(instrument, z) and check_dqs) or
                             (instrument == "cos" and not
                              is_bad_dq(instrument, z) and check_dqs) or
                             (instrument == "miri" and not
                              is_bad_dq(instrument, z) and check_dqs) or
                             (not check_dqs))
                        for x, y, z in zip(flux_values, flux_err_values, dqs)]
    else:
        bool_results = [False] * len(flux_values)

    return bool_results
#--------------------

#--------------------
def edge_trim(instrument, fluxes, fluxerrs, dqs, n_consecutive, median_flux,
              flux_scale_factor, median_fluxerr, fluxerr_scale_factor,
              fluxerr_95th):
    """
    Returns start and end indexes (end indexes are negatively indexed) of the
    best part of the spectrum to use when defining the plot's wavelength
    ranges.  Returns two sets of start and end indexes: one set without
    taking into account DQ flags, and one set that does take into account DQ
    flags (since some spectra have all DQ flags set > 0).

    :param instrument: The instrument that is being tested.

    :type instrument: str

    :param fluxes: The fluxes to be plotted.

    :type fluxes: numpy.ndarray

    :param fluxerrs: The uncertainties of the fluxes to be plotted.

    :type fluxerrs: numpy.ndarray

    :param dqs: The DQ flags of the spectrum.  For COS, these are the DQ_WGT
        bits from the header.

    :type dqs: numpy.ndarray

    :param n_consecutive: The consecutive number of fluxes that must belong to
        the "best" part of the spectrum for the start/end index to count.

    :type n_consecutive: int

    :param median_flux: The median flux, used in determining where the best part
        of the spectrum is.

    :type median_flux: float

    :param flux_scale_factor: Max. allowed ratio between the flux and a median
        flux value, used in edge trimming.

    :type flux_scale_factor: float

    :param median_fluxerr: The median flux uncertainty, used in determining
        where the best part of the spectrum is.

    :type median_fluxerr: float

    :param fluxerr_scale_factor: Max. allowed ratio between the flux uncertainty
        and a median flux uncertainty value, used in edge trimming.

    :type fluxerr_scale_factor: float

    :param fluxerr_95th: The flux uncertainty corresponding to the 95th
        percentile.

    :type fluxerr_95th: float

    :returns: int tuple -- Indexes that define the best part of the spectrum,
        in the order of (start_index_nodq, end_index_nodq, start_index_withdq,
        end_index_withdq).
    """

    # How many fluxes are there in total?
    n_fluxes = len(fluxes)

    # This defines the size of the spectrum at the edges to examine as a
    # "cheat", without having to look at the entire spectrum to find the good
    # boundaries.  It will look at the whole spectrum as a fallback option.
    edge_size = 100

    # SECTION FOR INDICES IGNORING DQ VALUES

    # First run the xrange test on the first and last 100 elements (for
    # speed).  If no good indices are found within those ranges, then and only
    # then will we run it on the entire array. """
    where_good_fluxes_nodq = numpy.where(numpy.asarray(_set_plot_xrange_test(
        instrument, numpy.concatenate((fluxes[0:edge_size],
                                       fluxes[-1*edge_size:])),
        numpy.concatenate((fluxerrs[0:edge_size],
                           fluxerrs[-1*edge_size:])), median_flux,
        flux_scale_factor, median_fluxerr, fluxerr_scale_factor,
        fluxerr_95th, numpy.concatenate((dqs[0:edge_size], dqs[-1*edge_size:])),
        check_fluxes=True, check_flux_ratios=False,
        check_flux_err_ratios=True, check_flux_err_percentile=False,
        check_dqs=False)))[0]
    if len(where_good_fluxes_nodq) > 0:
        where_good_fluxes_nodq[where_good_fluxes_nodq > edge_size-1] += (
            len(fluxes) - 2*edge_size)
        # Get the start and end indices, ignoring DQ values.
        start_index_nodq, end_index_nodq = find_good_indices(
            where_good_fluxes_nodq, n_consecutive, n_fluxes, first_pass=True)
    else:
        # If there were no results returned from the where call above, make
        # sure the start and end indices will fail the test.
        start_index_nodq = edge_size
        end_index_nodq = -2*edge_size-1

    # If the start and stop indices are not located within the range at the
    # end, run the where on the entire array.
    if (start_index_nodq > edge_size-1 or end_index_nodq < -2*edge_size+
            edge_size+1):
        where_good_fluxes_nodq = numpy.where(numpy.asarray(
            _set_plot_xrange_test(instrument, fluxes, fluxerrs, median_flux,
                                  flux_scale_factor, median_fluxerr,
                                  fluxerr_scale_factor, fluxerr_95th, dqs,
                                  check_fluxes=True,
                                  check_flux_ratios=False,
                                  check_flux_err_ratios=True,
                                  check_flux_err_percentile=False,
                                  check_dqs=False)))[0]
        start_index_nodq, end_index_nodq = find_good_indices(
            where_good_fluxes_nodq, n_consecutive, n_fluxes)

    # SECTION FOR INDICES INCLUDING DQ VALUES

    # First run the xrange test on the first and last 100 elements (for
    # speed).  If no good indices are found within those ranges, then and only
    # then will we run it on the entire array.
    where_good_fluxes_withdq = numpy.where(numpy.asarray(_set_plot_xrange_test(
        instrument, numpy.concatenate((fluxes[0:edge_size],
                                       fluxes[-1*edge_size:])),
        numpy.concatenate((fluxerrs[0:edge_size],
                           fluxerrs[-1*edge_size:])), median_flux,
        flux_scale_factor, median_fluxerr, fluxerr_scale_factor,
        fluxerr_95th, numpy.concatenate((dqs[0:edge_size],
                                         dqs[-1*edge_size:])),
        check_fluxes=True, check_flux_ratios=False,
        check_flux_err_ratios=True, check_flux_err_percentile=False,
        check_dqs=True)))[0]
    if len(where_good_fluxes_withdq) > 0:
        where_good_fluxes_withdq[where_good_fluxes_withdq > edge_size-1] += (
            len(fluxes) - 2*edge_size)
        # Get the start and end indices, including DQ values.
        start_index_withdq, end_index_withdq = find_good_indices(
            where_good_fluxes_withdq, n_consecutive, n_fluxes, first_pass=True)
    else:
        # If there were no results returned from the where call above, make
        # sure the start and end indices will fail the test.
        start_index_withdq = edge_size
        end_index_withdq = -2*edge_size-1

    # If the start and stop indices are not located within the range at the
    # end, run the where on the entire array.
    if (start_index_withdq > edge_size-1 or end_index_withdq < -2*edge_size+
            edge_size+1):
        where_good_fluxes_withdq = numpy.where(numpy.asarray(
            _set_plot_xrange_test(instrument, fluxes, fluxerrs, median_flux,
                                  flux_scale_factor, median_fluxerr,
                                  fluxerr_scale_factor, fluxerr_95th, dqs,
                                  check_fluxes=True,
                                  check_flux_ratios=False,
                                  check_flux_err_ratios=True,
                                  check_flux_err_percentile=False,
                                  check_dqs=True)))[0]
        start_index_withdq, end_index_withdq = find_good_indices(
            where_good_fluxes_withdq, n_consecutive, n_fluxes)

    return (start_index_nodq, end_index_nodq, start_index_withdq,
            end_index_withdq)
#--------------------

#--------------------
def find_good_indices(indices_arr, n_consecutive, n_fluxes, first_pass=False):
    """
    Given an array of indices, determines the start and stop indices of the
    `best` part of the spectrum based on the required number of consecutive
    `good` fluxes.

    :param indices_arr: Something

    :type indices_arr: list

    :param n_consecutive: The consecutive number of fluxes that must belong to
        the "best" part of the spectrum for the start/end index to count.

    :type n_consecutive: int

    :param n_fluxes: The total number of data points in the spectrum.

    :type n_fluxes: int

    :param first_pass: Used to determine if this is the first attempt to find
        good start and stop indices or not.  The first pass normally uses only a
        part of the spectrum (pieces at the beginning and end), while the next
        pass uses the entire spectrum.  If it's not the first pass and no good
        start and stop indices are found, then it must return the last and first
        index.

    :type first_pass: bool

    :returns: int tuple -- Indices that define the best part of the spectrum
        (start and stop).
    """

    # How many total good indices are there to work with?
    n_good_fluxes = len(indices_arr)

    # Handle the trivial case where the number of valid points is less than
    # the `n_consecutive` desired.  In this case, the start and end indices are
    # set to the last and first indices, simulating what would happen if the
    # entire array was traversed.
    if n_good_fluxes < n_consecutive:
        if not first_pass:
            start_index = 0
            end_index = -1
        else:
            start_index = n_fluxes
            end_index = -1 * n_fluxes
    else:
        # The start index is then the min. index for which the
        # `n_consecutive`'th index is exactly equal to `n_consecutive`-1.
        start_indexes_with_good_diffs = numpy.asarray(
            [x for i, x in enumerate(indices_arr[0:n_good_fluxes-n_consecutive+
                                                 1]) if
             indices_arr[i+n_consecutive-1]-indices_arr[i] == n_consecutive-1])

        if len(start_indexes_with_good_diffs) > 0:
            start_index = numpy.min(start_indexes_with_good_diffs)
            end_index = (numpy.max(start_indexes_with_good_diffs) - n_fluxes +
                         n_consecutive - 1)
        else:
            if not first_pass:
                start_index = 0
                end_index = -1
            else:
                start_index = n_fluxes
                end_index = -1 * n_fluxes

    return start_index, end_index
#--------------------
