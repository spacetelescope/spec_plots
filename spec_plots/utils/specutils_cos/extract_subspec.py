"""
.. module:: extract_subspec
   :synopsis: Extracts a sub-spectrum from the provided COS spectrum object.
       Modifies the COS spectrum object to replace that segment's spectrum with
       the specified sub-spectrum.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

#--------------------
# Built-In Imports
#--------------------
from __future__ import absolute_import
from __future__ import print_function
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
def extract_subspec(cos_spectrum, segment, min_wl=None, max_wl=None):
    """
    Extracts a sub-spectrum of a COS segment's spectrum.  If both min. and max.
    wavelengths are specified, then it keeps the part of the spectrum within
    those bounds.  If only a min. or max. is specified, then those are treated
    as lower or upper bounds, and the part of the spectrum redder/bluer than the
    bound is retained.

    :param cos_spectrum: COS spectrum as returned by READSPEC.

    :type cos_spectrum: COSSpectrum

    :param segment: The segment to do the extraction on.

    :type segment: str

    :param min_wl: The minimum wavlength value to keep.  If None, then there
        will be no lower bound.

    :type min_wl: float

    :param max_wl: The maximum wavelength value to keep.  If None, then there
        will be no upper bound.

    :type max_wl: float
    """

    if segment in cos_spectrum.segments:
        if min_wl is not None and max_wl is not None:
            where_within = numpy.where(
                (cos_spectrum.segments[segment].wavelengths >= min_wl) &
                (cos_spectrum.segments[segment].wavelengths <= max_wl))[0]
        elif min_wl is not None and max_wl is None:
            where_within = numpy.where(
                (cos_spectrum.segments[segment].wavelengths >= min_wl))[0]
        elif min_wl is None and max_wl is not None:
            where_within = numpy.where(
                (cos_spectrum.segments[segment].wavelengths <= max_wl))[0]

        n_within = len(where_within)
        if n_within > 0:
            # Then we extract the subspectrum for this object by modifying the
            # cos_spectrum object.
            cos_spectrum.segments[segment].nelem = n_within
            cos_spectrum.segments[segment].wavelengths = (
                cos_spectrum.segments[segment].wavelengths[where_within])
            cos_spectrum.segments[segment].fluxes = (
                cos_spectrum.segments[segment].fluxes[where_within])
            cos_spectrum.segments[segment].fluxerrs = (
                cos_spectrum.segments[segment].fluxerrs[where_within])
            cos_spectrum.segments[segment].dqs = (
                cos_spectrum.segments[segment].dqs[where_within])
        else:
            print("*** WARNING in SPECUTILS_COS: Requested subspectrum does"
                  " not overlap with this segment's spectrum.  No extraction"
                  " will be done.")
    else:
        raise SpecUtilsError("The segment where you want to perform the"
                             " subspectrum extraction is not present."
                             "  Specified \"" + segment +
                             "\", available segments for this spectrum are: " +
                             ', '.join(cos_spectrum.segments)+".")
#--------------------
