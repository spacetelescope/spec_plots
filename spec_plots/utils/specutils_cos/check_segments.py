__version__ = '1.32.0'

"""
.. module:: check_segments

   :synopsis: Checks that the array of segment names in the COS FITS extension table are some combination of allowed values.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

import numpy

""" <DEVEL> Note that this hack to make it so that the user can import `check_segments` directly as a module or run it from the command line as __main__ has the side effect of importing this module twice, despite my best efforts to work around it.  I don't think it will be a major issue, but worth thinking about in the future. </DEVEL> """
if __package__ is None:
    import sys, os
    specutils_cos_dir = os.path.dirname(os.path.abspath(__file__))
    utils_dir = os.path.dirname(specutils_cos_dir)
    parent_dir = os.path.dirname(utils_dir)
    sys.path.insert(1, parent_dir)
    import utils
    __package__ = str("utils.specutils")
    __name__ = str(__package__+"."+__name__)
    del sys, os

from ..specutils import SpecUtilsError

#--------------------

def check_segments(segments_list, input_file):
    """
    Checks that the array of "segments" in the COS spectrum header are expected values.  It returns a scalar string representing the band (either "FUV" or "NUV").  If the array is not what's expected, an Exception is raised.

    :param segments_list: List of segment labels.

    :type segments_list: list

    :param input_file: Name of input FITS file.

    :type input_file: str

    :returns: str -- A string representation of the band, either "FUV" or "NUV".

    :raises: utils.specutils.SpecUtilsError

    .. note::

         This function does not attempt to access the input file, it only requires the name of the file for error reporting purposes.
    """

    """ Get number of segments. """
    n_segments = len(segments_list)

    """ Segment array must contain either one, two, or three elements. """
    if n_segments == 1:
        """ This can happen with FUV data for some reason, so it can be either "FUVA" or "FUVB". """
        if numpy.array_equal(segments_list, ["FUVA"]) or numpy.array_equal(segments_list, ["FUVB"]):
            this_band = "FUV"
        else:
            raise SpecUtilsError("The array of SEGMENT strings contains 1 value, but is not equal to [\"FUVA\"] or [\"FUVB\"] in file " + input_file)

    elif n_segments == 2:
        """ Must be ["FUVA", "FUVB"]. """
        if numpy.array_equal(segments_list, ["FUVA", "FUVB"]):
            this_band = "FUV"
        else:
            raise SpecUtilsError("The array of SEGMENT strings contains 2 values, but is not equal to [\"FUVA\", \"FUVB\"] in file " + input_file)

    elif n_segments == 3:
        """ Must be ["NUVA", "NUVB", "NUVC"]. """
        if numpy.array_equal(segments_list, ["NUVA", "NUVB", "NUVC"]):
            this_band = "NUV"
        else:
            raise SpecUtilsError("The array of SEGMENT strings contains 3 values, but is not equal to [\"NUVA\", \"NUVB\", \"NUVC\"] in file " + input_file)

    else:
        raise SpecUtilsError("The array of SEGMENT strings should contain 1, 2, or 3 values, found " + str(n_segments) + " in file " + input_file)

    return this_band

#--------------------
