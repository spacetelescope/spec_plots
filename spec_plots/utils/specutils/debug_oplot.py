"""
.. module:: debug_oplot
   :synopsis: Plots spectral data points color-coded based on a suite of
       rejection criteria, primarily for debugging purposes.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

#--------------------
# Built-In Imports
#--------------------
from __future__ import absolute_import
from __future__ import division
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

def debug_oplot(this_plotarea, instrument, all_wls, all_fls, all_flerrs,
                all_dqs, median_flux, median_fluxerr, flux_scale_factor,
                fluxerr_scale_factor, fluxerr_95th, oplot_percentiles=False):
    """
    Creates plots of the spectra with color-coding and special annotation to
    identify which points were rejected by which tests.  Useful for
    debugging and understanding why a given plot had its plot axes defined
    the way it did.

    :param this_plotarea: The AxesSubplot object to plot on.

    :type this_plotarea: matplotlib.axes._subplots.AxesSubplot

    :param instrument: The instrument that is being tested.

    :type instrument: str

    :param all_wls: Array of wavelengths.

    :type all_wls: numpy.ndarray

    :param all_fls: Array of fluxes.

    :type all_fls: numpy.ndarray

    :param all_flerrs: Array of flux uncertainties.

    :type all_flerrs: numpy.ndarray

    :param all_dqs: Array of data quality flags.  For COS, these are the DQ_WGT
        bits from the header.

    :type all_dqs: numpy.ndarray

    :param median_flux: The median flux used in determining where the best part
        of the spectrum is.

    :type median_flux: float

    :param median_fluxerr: The median flux uncertainty used in determining where
        the best part of the spectrum is.

    :type median_fluxerr: float

    :param flux_scale_factor: Max. allowed ratio between the flux and a median
        flux value used in edge trimming.

    :type flux_scale_factor: float

    :param fluxerr_scale_factor: Max. allowed ratio between the flux uncertainty
        and a median flux uncertainty value used in edge trimming.

    :type fluxerr_scale_factor: float

    :param fluxerr_95th: The flux uncertainty corresponding to the 95th
        percentile.

    :type fluxerr_95th: float

    :param oplot_percentiles: Set this to True to overplot points where the flux
        uncertainties are greater than the 95th percentile.  Default = False.

    :type oplot_percentiles: bool
    """

    # Plot all of the points in black.
    this_plotarea.errorbar(all_wls, all_fls, yerr=all_flerrs, ecolor='c',
                           color='k', label='Passed')

    # Plot those fluxes that fail because they are much greater than the median
    # in blue.
    if numpy.isfinite(median_flux):
        where_fluxtoolarge = numpy.where(abs(all_fls/median_flux) >
                                         flux_scale_factor)
        if len(where_fluxtoolarge[0]) > 0:
            this_plotarea.plot(all_wls[where_fluxtoolarge],
                               all_fls[where_fluxtoolarge], 'bo',
                               label="Flux>>Median")

    # Plot those fluxes that fail because they are exactly equal to 0 in green.
    where_allzero = numpy.where(all_fls == 0.0)
    if len(where_allzero[0]) > 0:
        this_plotarea.plot(all_wls[where_allzero], all_fls[where_allzero],
                           'go', label="Flux=0")

    # Plot those fluxes that fail because they have a bad DQ value in red.
    where_bad_dq = numpy.where(is_bad_dq(instrument, all_dqs))[0]
    if len(where_bad_dq) > 0:
        this_plotarea.plot(all_wls[where_bad_dq], all_fls[where_bad_dq], 'ro',
                           label="BAD DQ_WGT")

    # Plot those fluxes that fail because their flux uncertainties are much
    # greater than the median flux uncertainty in magenta.
    if numpy.isfinite(median_fluxerr):
        where_bigerr = numpy.where(all_flerrs/median_fluxerr >
                                   fluxerr_scale_factor)
        if len(where_bigerr[0]) > 0:
            this_plotarea.plot(all_wls[where_bigerr], all_fls[where_bigerr],
                               'mo', label="FluxErr>>Median")

    # Plot those fluxes that fail because their flux uncertainties are greater
    # than the 95th percentile in yellow.
    if oplot_percentiles:
        where_bigerrpercentile = numpy.where(all_flerrs > fluxerr_95th)
        if len(where_bigerrpercentile[0]) > 0:
            this_plotarea.plot(all_wls[where_bigerrpercentile],
                               all_fls[where_bigerrpercentile], 'yo',
                               label="FluxErr>95th %")

    this_plotarea.legend(loc="upper center", ncol=4)
#--------------------
