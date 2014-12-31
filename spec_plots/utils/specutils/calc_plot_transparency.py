__version__ = '1.31'

"""
.. module:: calc_plot_transparency

   :synopsis: Calculates the transparency ("alpha") value to be used in the plot, based on the number of "peaks" in the data.  A peak a counted if the data point before and after it are less than the value being considered.  A minimum transparency is used if the number of peaks is half the total number of points (largest number of peaks possible), and a value of 1.0 (no transparency) is used if the number of peaks is zero (it will almost never be exactly 0).  Values between these two boundaries are scaled linearly between the minimum transparency value and 1.0.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

#--------------------

def calc_plot_transparency(path_length_sq):
    """
    Calculates the plot transparency (alpha value) to use, based on the total path length of the line segments that make up the spectrum.

    :param path_length_sq: The square of the path length for all line segments in this spectrum.

    :type path_length_sq: float

    :returns: float -- The transparency (alpha) value to use in the plot.
    """

    """ This is the minimum transparency value allowed. """
    min_alpha = 0.1

    """ Calculate the alpha value to return. """
    if path_length_sq > 50.:
        return_alpha = 0.25
    elif path_length_sq <= 50. and path_length_sq > 1.5:
        return_alpha = 0.25
    else:
        return_alpha = 0.25

    return return_alpha

#--------------------
