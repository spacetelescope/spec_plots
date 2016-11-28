"""
.. module:: check_segments
   :synopsis: Checks that the array of segment names in the COS FITS extension
       table are some combination of allowed values.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

#--------------------
# Built-In Imports
#--------------------
from __future__ import absolute_import
import os
import sys
from builtins import str
#--------------------
# External Imports
#--------------------
import numpy
#--------------------
# Package Imports
#--------------------
from spec_plots.utils.specutils.specutilserror import SpecUtilsError
from spec_plots import __version__
#--------------------

#--------------------
def check_segments(segments_list, input_file):
    """
    Checks that the array of "segments" in the COS spectrum header are expected
    values.  It returns a scalar string representing the band (either "FUV"
    or "NUV").  If the array is not what's expected, an Exception is raised.

    :param segments_list: List of segment labels.

    :type segments_list: list

    :param input_file: Name of input FITS file.

    :type input_file: str

    :returns: str -- A string representation of the band, either "FUV" or "NUV".

    :raises: utils.specutils.SpecUtilsError

    .. note::

         This function does not attempt to access the input file, it only
         requires the name of the file for error reporting purposes.
    """

    # Get number of segments.
    n_segments = len(segments_list)

    # Segment array must contain either one, two, or three elements.
    if n_segments == 1:
        # This can happen with FUV data for some reason, so it can be either
        # "FUVA" or "FUVB".
        if (numpy.array_equal(segments_list, ["FUVA"]) or
                numpy.array_equal(segments_list, ["FUVB"])):
            this_band = "FUV"
        else:
            raise SpecUtilsError("The array of SEGMENT strings contains one"
                                 " value, but is not equal to [\"FUVA\"] or"
                                 " [\"FUVB\"] in file " + input_file)

    elif n_segments == 2:
        # Must be ["FUVA", "FUVB"].
        if numpy.array_equal(segments_list, ["FUVA", "FUVB"]):
            this_band = "FUV"
        else:
            raise SpecUtilsError("The array of SEGMENT strings contains two"
                                 " values, but is not equal to [\"FUVA\", "
                                 "\"FUVB\"] in file " + input_file)

    elif n_segments == 3:
        # Must be ["NUVA", "NUVB", "NUVC"].
        if numpy.array_equal(segments_list, ["NUVA", "NUVB", "NUVC"]):
            this_band = "NUV"
        else:
            raise SpecUtilsError("The array of SEGMENT strings contains three"
                                 " values, but is not equal to [\"NUVA\","
                                 " \"NUVB\", \"NUVC\"] in file " + input_file)

    else:
        raise SpecUtilsError("The array of SEGMENT strings should contain 1,"
                             " 2, or 3 values, found " + str(n_segments) +
                             " in file " + input_file)

    return this_band
#--------------------
