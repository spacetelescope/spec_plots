"""
.. module:: calc_plot_metrics
   :synopsis: Calculates metrics related to a plot, including flux statistics,
       avoid regions, and optimal x- and y-axis ranges for plots.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

#--------------------
# Built-In Imports
#--------------------
from __future__ import absolute_import
#--------------------
# External Imports
#--------------------
import matplotlib.collections
import numpy
#--------------------
# Package Imports
#--------------------
from spec_plots.utils.specutils.get_flux_stats import get_flux_stats
from spec_plots.utils.specutils.avoidregion import generate_avoid_regions
from spec_plots.utils.specutils.set_plot_xrange import set_plot_xrange
from spec_plots.utils.specutils.set_plot_yrange import set_plot_yrange
from spec_plots import __version__

#--------------------

def calc_plot_metrics(instrument, wls, fls, flerrs, dqs, n_consecutive,
                      flux_scale_factor, fluxerr_scale_factor):
    """
    Calculates a variety of plot metrics, including flux statistics, avoid
    regions, and optimal x- and y-axis ranges.

    :param instrument: The instrument that is being tested.

    :type instrument: str

    :param wls: The wavelengths to be plotted.

    :type wls: numpy.ndarray

    :param fls: The fluxes to be plotted.

    :type fls: numpy.ndarray

    :param flerrs: The uncertainties of the fluxes to be plotted.

    :type flerrs: numpy.ndarray

    :param dqs: The DQ flags of the spectrum to be plotted.  For COS, these are
        the DQ_WGT bits from the header.

    :type dqs: numpy.ndarray

    :param n_consecutive: How many consecutive points must pass the test for
        the index to count as the valid start/end of the spectrum?

    :type n_consecutive: int

    :param flux_scale_factor: Max. allowed ratio between the flux and a median
        flux value, used in edge trimming.

    :type flux_scale_factor: float

    :param fluxerr_scale_factor: Max. allowed ratio between the flux uncertainty
        and a median flux uncertainty value, used in edge trimming.

    :type fluxerr_scale_factor: float

    :returns: dict -- A container for the plot metrics.
    """

    # Calculate statistics on the fluxes.
    median_flux, median_fluxerr, fluxerr_95th = get_flux_stats(fls, flerrs)

    # Determine optimal x-axis.  This is not the x-axis plot range used, but
    # rather the area of the plot that is considered when scaling the y-axis.
    optimal_xaxis_range = set_plot_xrange(instrument, wls, fls, flerrs, dqs,
                                          n_consecutive, flux_scale_factor,
                                          fluxerr_scale_factor, median_flux,
                                          median_fluxerr, fluxerr_95th)

    # Create COS avoid regions.
    avoid_regions = generate_avoid_regions(instrument)

    # Determine the optimal y-axis.
    if all(numpy.isfinite(optimal_xaxis_range)):
        y_axis_range = set_plot_yrange(wls, fls, avoid_regions=avoid_regions,
                                       wl_range=optimal_xaxis_range)
    else:
        y_axis_range = [numpy.nan, numpy.nan]

    # Construct the LineCollection of segments for this curve.   Only do so if
    # the plot transparency will be < 1.0.
    points = numpy.array([wls, fls]).T.reshape(-1, 1, 2)
    segments = numpy.concatenate([points[:-1], points[1:]], axis=1)
    linecoll = matplotlib.collections.LineCollection(segments)

    # Return the plot_metrics dict.
    return {"median_flux":median_flux, "median_fluxerr":median_fluxerr,
            "fluxerr_95th":fluxerr_95th, "optimal_xaxis_range":
                optimal_xaxis_range, "avoid_regions":avoid_regions,
            "y_axis_range":y_axis_range, "line_collection":linecoll}
#--------------------

