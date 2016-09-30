"""
.. module:: get_flux_stats
   :synopsis: Calculates statistics on an arrays of fluxes and corresponding
       flux uncertainties, such as the median values and the 95th percentile
       values.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

#--------------------
# Built-In Imports
#--------------------
from __future__ import absolute_import
#--------------------
# External Imports
#--------------------
import numpy
#--------------------
# Package Imports
#--------------------
from spec_plots import __version__

#--------------------

def get_flux_stats(fluxes, fluxerrs):
    """
    Calculates median values of the fluxes and flux uncertainties, and the 95th
    percentile flux uncertainty.

    :param fluxes: The fluxes to be plotted.

    :type fluxes: numpy.ndarray

    :param fluxerrs: The uncertainties of the fluxes to be plotted.

    :type fluxerrs: numpy.ndarray

    :returns: tuple -- The median flux, median flux uncertainty, and the 95th
        percentile flux uncertainty.
    """

    # Find the median flux value, ignoring any NaN values or fluxes that are
    # 0.0.
    where_finite_and_notzero = numpy.where((numpy.isfinite(fluxes)) &
                                           (fluxes != 0.0))

    if len(where_finite_and_notzero[0]) > 0:
        median_flux = numpy.median(fluxes[where_finite_and_notzero])
        median_fluxerr = numpy.median(fluxerrs[where_finite_and_notzero])
    else:
        median_flux = numpy.nan
        median_fluxerr = numpy.nan

    # Get the 95th percentile flux uncertainty value.
    fluxerr_95th = numpy.percentile(fluxerrs, 95.)

    return median_flux, median_fluxerr, fluxerr_95th
#--------------------
