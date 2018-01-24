"""
.. module:: calc_covering_fraction
   :synopsis: Calculates the ratio of the number of pixels containing plot data
       vs. empty (background) plot pixels.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

#--------------------
# Built-In Imports
#--------------------
from __future__ import absolute_import
from __future__ import division
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

def count_blue_red(buf):
    """
    Counts the number of "blue" vs. "red" pixels in the plot area.

    :param buf: The plot pixels read into a buffer of RGB values (
        one-dimensional, [r,g,b,r,g,b,...]).

    :type buf: numpy.ndarray

    :returns: tuple -- The number of blue and red pixels, respectively.  Blue
        represent pixels that have plot data in them, red represent background
        pixels that do not have any plot lines going through them.
    """

    blue_count = 0
    red_count = 0
    for i in range(0, len(buf), 3):
        # First make sure this isn't a grey pixel, where R=G=B.
        if not (buf[i] == buf[i+1] and buf[i] == buf[i+2]):
            if buf[i] > 0 and buf[i+2] == 0:
                red_count += 1
            elif buf[i+2] > 0:
                blue_count += 1

    return blue_count, red_count
#--------------------

#--------------------
def calc_covering_fraction(fig, subplots, subplot_num, optimize=True):
    """
    Calculates the ratio of pixels that contain plot data vs. background (empty)
    pixels by temporarily filling the background with RGB value = red and
    then identifying blue (plot data) pixels.  This ratio is used to
    determine whether to make the plot lines transparent or not for readability.

    :param fig: The plot figure.

    :type fig: Figure

    :param subplots: The list of subplots.

    :type subplots: list

    :param subplot_num: The index indicating which subplot is currently being
        examined.

    :type subplot_num: int

    :param optimize: If set to True, will use a slightly optimized version of
        determining the plot covering fraction.  This has been tested to
        *almost* always give the same result as the non-optimized version.
        There has been (rare) instances where the covering fractions differed
        slightly, so if the preference is to have more accurate/robust covering
        fractions, set optimize=False!  The speed benefit is on the order of
        10-20%, but is largest for plots with multiple subplots.

    :type optimize: bool

    :returns: float -- The ratio of pixels containing plot data ("blue") vs.
        empty ("red") pixels, as a percentage (0. <= p <= 100.).
    """

    # Turn off all other subplots since the RGB array is for the entire
    # figure, not just this subplot part of it.
    for i, subplot in enumerate(subplots):
        if i != subplot_num:
            subplot.set_visible(False)
        else:
            subplot.set_visible(True)

    # Temporarily make the background color red so that it's easy to identify
    # empty pixels.
    subplots[subplot_num].set_facecolor('red')

    # Draw the plot, convert into numpy array of RGB values.
    fig.canvas.draw()
    buf = numpy.fromstring(fig.canvas.tostring_rgb(), dtype=numpy.uint8)

    # Count the number of "blue" and "red" pixels.  The optimized version
    # attempts to split the RGB array into equal slices so that the counting of
    # red and blue pixels can be done quicker.  Early tests have shown this
    # usually produces the same results.
    if optimize:
        n_per_subplot = int(numpy.ceil(len(buf) / len(subplots)))
        while n_per_subplot % 3 != 0:
            n_per_subplot += 1
        startindex = subplot_num * n_per_subplot
        endindex = startindex + n_per_subplot
        if subplot_num < len(subplots)-1:
            blue_count, red_count = count_blue_red(buf[startindex:endindex])
        else:
            blue_count, red_count = count_blue_red(buf[startindex:])
    else:
        blue_count, red_count = count_blue_red(buf)

    # If this is the last subplot, then make sure we make all the other
    # subplots visible again.
    if subplot_num == subplots.shape[0]-1:
        for subplot in subplots:
            subplot.set_visible(True)

    # Make the background color white again for this subplot.
    subplots[subplot_num].set_facecolor('white')

    return float(blue_count) / float(red_count+blue_count) * 100.
#--------------------
