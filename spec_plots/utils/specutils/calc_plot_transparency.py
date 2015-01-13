__version__ = '1.31'

"""
.. module:: calc_plot_transparency

   :synopsis: Calculates the transparency ("alpha") value to be used in the plot, based on the number of "peaks" in the data.  A peak a counted if the data point before and after it are less than the value being considered.  A minimum transparency is used if the number of peaks is half the total number of points (largest number of peaks possible), and a value of 1.0 (no transparency) is used if the number of peaks is zero (it will almost never be exactly 0).  Values between these two boundaries are scaled linearly between the minimum transparency value and 1.0.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

import numpy

#--------------------

def calc_plot_transparency(wls, fls, y_axis_range):
    """
    Calculates the plot transparency (alpha value) to use, based on the y-axis range compared to the average scatter in the spectrum.

    :param wls: The wavelengths to be plotted.

    :type wls: numpy.ndarray

    :param fls: The fluxes to be plotted.

    :type fls: numpy.ndarray

    :param y_axis_range: The optimal y-axis plot range for this spectrum.

    :type path_length_sq: list

    :returns: float -- The transparency (alpha) value to use in the plot.
    """

    """ This is the minimum transparency value allowed. """
    min_alpha = 0.1

    """ Compute the y-range. """
    yrange = y_axis_range[1]-y_axis_range[0]

    """ What is the median difference between adjacent fluxes? """
    flux_deltas = abs(numpy.subtract(fls[1:], fls[0:-1]))
    median_flux_delta = numpy.median(flux_deltas)

    """ Define the y-range ratio as the ratio between the median_flux_delta and the yrange. """
    y_range_ratio = median_flux_delta / yrange
    y_range_ratio_alpha = numpy.percentile(flux_deltas, 45.) / yrange
    y_range_ratio_beta =  numpy.percentile(flux_deltas, 55.) / yrange

    """ Calculate the alpha value to return. """
    if y_range_ratio_alpha > 0.08 and y_range_ratio_beta > 0.08 and y_range_ratio - y_range_ratio_alpha > 0.01 and y_range_ratio_beta - y_range_ratio > 0.01:
        return_alpha = min_alpha
    else:
        return_alpha = 1.0

    return return_alpha, y_range_ratio

#--------------------
