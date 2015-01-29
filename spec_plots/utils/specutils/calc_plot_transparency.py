__version__ = '1.32.0'

"""
.. module:: calc_plot_transparency

   :synopsis: Calculates the transparency ("alpha") value to be used in the plot.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

import numpy
import warnings

#--------------------

def calc_flux_crossings(fls):
    """
    Counts the number of flux crossings through the median flux value.

    :param fls: The fluxes to count the number of crossing on.

    :type fls: list

    :returns: int -- The number of median flux crossings.
    """

    """ Compute the median value for this section of spectrum, ignoring NaN values. """
    med_fl = numpy.nanmedian(fls)

    """ Count the number of crossings.  A crossing is counted once the sign of the flux - median value changes (going from over/under the median to under/over the median).  If the flux is exactly 0.0 then it is not counted as a crossing. """
    bin_crosses_total = 0
    cur_sign = numpy.sign(fls[0] - med_fl)
    for j in fls[1:]:
        this_sign = numpy.sign(j - med_fl)
        if cur_sign != 0:
            if this_sign != cur_sign:
                bin_crosses_total += 1
                cur_sign = this_sign
        else:
            cur_sign = this_sign

    return bin_crosses_total

#--------------------

def calc_plot_transparency(wls, fls, y_axis_range):
    """
    Calculates the plot transparency (alpha value) to use.

    :param wls: The wavelengths to be plotted.

    :type wls: numpy.ndarray

    :param fls: The fluxes to be plotted.

    :type fls: numpy.ndarray

    :param y_axis_range: The optimal y-axis plot range for this spectrum.

    :type path_length_sq: list

    :returns: float -- The transparency (alpha) value to use in the plot.
    """

    """ This is the minimum transparency value allowed. """
    min_alpha = 0.10

    """ Compute the y-range. """
    yrange = y_axis_range[1]-y_axis_range[0]
    midflux = yrange/2. + y_axis_range[0]
    upper_flux_bound = 0.80 * yrange + y_axis_range[0]

    """ Divide the fluxes into wavelength bins. """
    digitized_indices, unique_digitized_indices, med_npoints_per_bin = divide_into_bins(wls)
    n_uniq_bins = len(unique_digitized_indices)

    bin_crosses = [0] * n_uniq_bins
    """ How many times do the fluxes cross the local median in the bin? """
    for a,i in enumerate(unique_digitized_indices):
        these_wls = wls[numpy.where(digitized_indices == i)[0]]
        these_fls = fls[numpy.where(digitized_indices == i)[0]]
        bin_crosses_total = calc_flux_crossings(these_fls)        
        bin_crosses[a] = float(bin_crosses_total) / float(len(these_fls)) * 100.
    bin_crosses = numpy.asarray(bin_crosses)

    """ Count flux statistics """
    n_bins_gt_per = len(numpy.where(bin_crosses >= 25.)[0])
    n_bins_gt_per_frac = float(n_bins_gt_per) / float(n_uniq_bins) * 100.
    n_flux_gt_per = len(numpy.where(abs(fls - midflux) >= upper_flux_bound - midflux)[0])
    n_flux_gt_per_frac = float(n_flux_gt_per) / float(len(fls)) * 100.

    """ Calculate the alpha value to return. """
    if n_bins_gt_per_frac >= 75. and n_flux_gt_per_frac >= 15.:
        return_alpha = min_alpha
    else:
        return_alpha = 1.0

    return return_alpha, n_bins_gt_per_frac, n_flux_gt_per_frac

#--------------------

def divide_into_bins(wls):
    """
    Divides the input wavelengths into bins appropriate for calculating the number of zero crossings through the median flux.

    :param wls: The wavelengths to be plotted.

    :type wls: numpy.ndarray

    :returns: numpy.ndarray, numpy.ndarray, float -- Array of digitized indices for each wavelength (i.e., the bin it belongs to), an array of unique digitized indices, and the median number of data points in each bin.

    """

    """ Initialize some counter variables. """
    npoints_per_bin = [0] ; while_counter = 0.

    """ Use a binsize of 5 Angstroms to start with.  This works well for most COS spectra, but not always for STIS.  If there aren't enough points in each bin, increase the bin size by an order of magnitude until there are enough. """
    while numpy.median(numpy.asarray(npoints_per_bin)) <= 30:
        """ Define bin size. """
        binsize = 5. * 10.**while_counter

        bins = numpy.arange(numpy.nanmin(wls), numpy.nanmax(wls)+binsize, binsize)
        digitized_indices = numpy.digitize(wls, bins)
        unique_digitized_indices = numpy.unique(digitized_indices)
        
        """ Count the number of points in each bin. """
        npoints_per_bin = [0] * len(unique_digitized_indices)
        for a,i in enumerate(unique_digitized_indices):
            npoints_per_bin[a] = len(numpy.where(digitized_indices == i)[0])

        """ Increment counter in case we need to use a larger binsize. """
        while_counter += 1.

    return digitized_indices, unique_digitized_indices, numpy.median(numpy.asarray(npoints_per_bin))

#--------------------
