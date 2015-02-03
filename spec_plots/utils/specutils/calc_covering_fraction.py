__version__ = '1.32.0'

"""
.. module:: calc_covering_fraction

   :synopsis: Calculates the ratio of the number of pixels containing plot data vs. empty (background) plot pixels.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

import numpy

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
    w,h = fig.canvas.get_width_height()
    buf = numpy.fromstring( fig.canvas.tostring_rgb(), dtype=numpy.uint8 )
    buf.shape = (w, h, 3)
    buf = numpy.roll( buf, 2, axis=1)
    
    """ Now count the number of "blue" vs. "red" pixels in the entire figure (which should just be this subplot since we hid all the other ones. """
    blue_count = 0.
    red_count = 0.
    for dim1 in xrange(buf.shape[0]):
        for dim2 in xrange(buf.shape[1]):
            """ Any grey pixels have equal R, G, and B values, and almost always represnt the axis titles, boundaries, etc., and are ignored. """
            isGrey = buf[dim1,dim2,0] == buf[dim1,dim2,1] and buf[dim1,dim2,0] == buf[dim1,dim2,2]
            
            """ We define ANY data point that has ANY blue component and is not grey as containing data. """
            if buf[dim1,dim2,2] > 0 and not isGrey:
                blue_count += 1.
            elif buf[dim1,dim2,0] > 0 and buf[dim1,dim2,2] == 0 and not isGrey:
                """ We define any empty plot pixel as containing ANY red component but absolutely NO blue component as background/empty plot pixels.  Note that the grid lines (even if present) would count as red pixels in this definition. """
                red_count += 1.

    """ If this is the last subplot, then make sure we make all the other subplots visible again. """
    if subplot_num == subplots.shape[0]-1:
        for x in subplots:
            x.set_visible(True)

    """ Make the background color white again for this subplot. """
    subplots[subplot_num].set_axis_bgcolor('white')

    return blue_count / (red_count+blue_count) * 100.

#--------------------
