"""
.. module:: rms
   :synopsis: Calculates the root-mean-square of an array of input values.
       Optionally, a constant offset value can be subtracted before calculating
       the rms.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

#--------------------
# Built-In Imports
#--------------------
from __future__ import absolute_import
import math
#--------------------
# External Imports
#--------------------
import numpy
#--------------------
# Package Imports
#--------------------
from spec_plots import __version__

#--------------------

def rms(values, offset=0.):
    """
    Calculates the RMS about some offset (default offset is 0.)

    :param values: Array of values to compute the rms of.

    :type values: numpy.ndarray

    :param offset: Optional offset to compute the rms about.  Defaults to 0.0.

    :type offset: float

    :returns: float -- A scalar float containing the rms about the offset.
    """
    return math.sqrt(numpy.nanmean([(x-offset)**2 for x in values]))
#--------------------
