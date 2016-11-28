"""
.. module:: is_bad_dq
   :synopsis: Returns True if a given DQ (or array of DQ values) fail(s) a test
       defined in this module, otherwise returns False.  If an array of DQ
       values was given, then an array of True/False values is returned for each
       element in the array, otherwise a single True/False value is returned.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

#--------------------
# Built-In Imports
#--------------------
from __future__ import absolute_import
#--------------------
# External Imports
#--------------------
import numpy
#--------------------
# Package Imports
#--------------------
from spec_plots import __version__
#--------------------

def is_bad_dq(instrument, dqs):
    """
    Returns True/False whether the DQ values are from a good part of the
    spectrum.  If the input DQs is a scalar value, the return result is also
    a scalar value.  If the input DQs are a numpy array, the return result
    is also a numpy array.

    :param instrument: The instrument the DQs come from.

    :type instrument: str

    :param dqs: Array of DQ (STIS) or DQ_WGT (COS) values.

    :type dqs: numpy.ndarray

    :returns: bool or numpy.ndarray -- A scalar boolean or array of booleans.
    """

    if instrument == "cos":
        if isinstance(dqs, numpy.ndarray):
            return numpy.asarray([x < 1 for x in dqs])
        else:
            return dqs < 1

    elif instrument == "stis":
        if isinstance(dqs, numpy.ndarray):
            return numpy.asarray([x != 0 and x != 16 for x in dqs])
        else:
            return dqs != 0 and dqs != 16

    elif instrument in ["miri", "nirspec", "niriss"]:
        if isinstance(dqs, numpy.ndarray):
            return numpy.asarray([x < 1 for x in dqs])
        else:
            return dqs < 1

    else:
        return numpy.asarray([])
#--------------------
