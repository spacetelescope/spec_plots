"""
.. module:: set_plot_yrange
   :synopsis: Determines the optimal y-axis plot range.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

#--------------------
# Built-In Imports
#--------------------
from __future__ import absolute_import
from builtins import range
#--------------------
# External Imports
#--------------------
import numpy
#--------------------
# Package Imports
#--------------------
from spec_plots import __version__

#--------------------

def set_plot_yrange(wavelengths, fluxes, avoid_regions=None, wl_range=None):
    """
    Given an array of wavelengths, fluxes, and avoid regions, returns a list of
    [ymin,ymax] to define an optimal y-axis plot range.

    :param wavelengths: The wavelengths to be plotted.

    :type wavelengths: numpy.ndarray

    :param fluxes: The fluxes to be plotted.

    :type fluxes: numpy.ndarray

    :param avoid_regions: A list of wavelength ranges to avoid when calculating
        optimal y-axis plot range.

    :type avoid_regions: list of STISAvoidRegion objects.

    :param wl_range: The min. and max. wavelength that defines the x-axis plot
        range.  The default is None, in which case the min. and max. if the
        input wavelength array will be used.

    :type wl_range: list

    :returns: list -- Two-element list containing the optimal [ymin,ymax] values
        to define the y-axis plot range.

    .. note::

         This function makes use of an internal look-up table of wavelength
         regions where known contaminating emission lines or other strong UV
         artifacts can affect the zoom level of the plot.
    """

    # Default the wavelength range to be the entire spectrum if not specified.
    if wl_range is None:
        wl_range = [numpy.nanmin(wavelengths), numpy.nanmax(wavelengths)]

    # This list will keep track of which fluxes to retain when defining the
    # y-axis plot range.  Setting the value to 1 means keep this flux for
    # consideration.
    keep_indices = numpy.asarray([1] * len(wavelengths))

    # Set any fluxes to 0 (don't keep) if they fall within an Avoid Region.
    if avoid_regions is not None:
        for i, region in enumerate(avoid_regions):
            if i == 0:
                # If this is the first Avoid Region, then we need to check if
                # the wavelengths are within the specified bounds supplied
                # through the `wl_range` parameter, in addition to checking if
                # they are within the Avoid Region itself.
                reject_indices = [i for i in range(len(wavelengths)) if (
                    wavelengths[i] >= region.minwl and
                    wavelengths[i] <= region.maxwl or
                    wavelengths[i] < wl_range[0] or
                    wavelengths[i] > wl_range[1])]
            else:
                # Don't need to worry about checking wavelengths within bounds
                # after the first Avoid Region is examined.
                reject_indices = [i for i in range(len(wavelengths)) if (
                    wavelengths[i] >= region.minwl and
                    wavelengths[i] <= region.maxwl)]
            # Set indices that we don't want to keep to 0.  Note that if
            # reject_indices is an empty list then nothing will change.
            keep_indices[reject_indices] = 0

    # After all indices have been set to keep or reject, pull out just the
    # fluxes that should be kept.
    keep_fluxes = numpy.asarray([f for ii, f in enumerate(fluxes) if (
        keep_indices[ii] == 1 and numpy.isfinite(fluxes[ii]))])

    # Don't just take the pure min and max, since large outliers can affect
    # the calculation.  Instead, take the 1th and 99th percentile fluxes within
    # the region to calculate the `min` and `max` fluxes.
    min_flux = numpy.percentile(keep_fluxes, 1.)
    max_flux = numpy.percentile(keep_fluxes, 99.)

    # Determine a y-buffer based on the difference between the max. and min.
    # fluxes.
    ybuffer = 0.3 * (max_flux-min_flux)

    # Make sure the min. and max. fluxes aren't identical (both 0., or both
    # the same exact value.  If so, just return the min. and max. value nudged
    # by 1.0.
    if min_flux != max_flux:
        if min_flux - ybuffer > 0.:
            return [min_flux-ybuffer, max_flux+ybuffer]
        else:
            # We don't want the y-axis range to be TOO far negative, so limit
            # it to be close to the lowest data point. """
            return [1.1*min_flux, max_flux+ybuffer]
    else:
        return [min_flux-1., max_flux+1.]
#--------------------
