__version__ = '1.32.0'

"""
.. module:: calc_covering_fraction

   :synopsis: Calculates the ratio of the number of pixels containing plot data vs. empty (background) plot pixels.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

import numpy

#--------------------

def testisGrey(rgb):
    return rgb[0] == rgb[1] and rgb[0] == rgb[2]

#--------------------

def count_blue_red(buf):
    blue_count = 0 ; red_count = 0
    for ii in xrange(0,len(buf),3):
        if not testisGrey(buf[ii:ii+3]):
            if buf[ii] > 0 and buf[ii+2] == 0:
                red_count += 1
            elif buf[ii+2] > 0:
                blue_count += 1
    return blue_count, red_count

#--------------------

def calc_covering_fraction(fig, subplots, subplot_num):
    """
    Calculates the ratio of pixels that contain plot data vs. background (empty) pixels by temporarily filling the background with RGB value = red and then identifying blue (plot data) pixels.  This ratio is used to determine whether to make the plot lines transparent or not for readability.

    :param fig: The plot figure.

    :type fig: Figure

    :param subplots: The list of subplots.

    :type fls: 

    :param subplot_num: The index indicating which subplot is currently being examined.

    :type subplot_num: int

    :returns: float -- The ratio of pixels containing plot data ("blue") vs. empty ("red") pixels, as a percentage (0. <= p <= 100.).
    """

    """ Turn off all other subplots since the RGB array is for the entire figure, not just this subplot part of it. """
    for i,x in enumerate(subplots):
        if i != subplot_num:
            x.set_visible(False)
        else:
            x.set_visible(True)

    """ Temporarily make the background color red so that it's easy to identify empty pixels. """
    subplots[subplot_num].set_axis_bgcolor('red')

    """ Draw the plot, convert into numpy array of RGB values. """
    fig.canvas.draw()
    buf = numpy.fromstring( fig.canvas.tostring_rgb(), dtype=numpy.uint8 )

    """ Count the number of "blue" and "red" pixels. """
    blue_count, red_count = count_blue_red(buf)

    """ If this is the last subplot, then make sure we make all the other subplots visible again. """
    if subplot_num == subplots.shape[0]-1:
        for x in subplots:
            x.set_visible(True)

    """ Make the background color white again for this subplot. """
    subplots[subplot_num].set_axis_bgcolor('white')
    
    print blue_count, red_count

    return blue_count / (red_count+blue_count) * 100.

#--------------------
