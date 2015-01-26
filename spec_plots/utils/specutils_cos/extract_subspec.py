__version__ = '1.32.0'

"""
.. module:: extract_subspec

   :synopsis: Extracts a sub-spectrum from the provided COS spectrum object.  Modifies the COS spectrum object to replace that segment's spectrum with the specified sub-spectrum.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

import numpy

""" <DEVEL> Note that this hack to make it so that the user can import `extract_subspec` directly as a module or run it from the command line as __main__ has the side effect of importing this module twice, despite my best efforts to work around it.  I don't think it will be a major issue, but worth thinking about in the future. </DEVEL> """
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

def extract_subspec(cos_spectrum, segment, min_wl=None, max_wl=None):
    """
    Extracts a sub-spectrum of a COS segment's spectrum.  If both min. and max. wavelengths are specified, then it keeps the part of the spectrum within those bounds.  If only a min. or max. is specified, then those are treated as lower or upper bounds, and the part of the spectrum redder/bluer than the bound is retained.

    :param cos_spectrum: COS spectrum as returned by READSPEC.

    :type cos_spectrum: COSSpectrum

    :param segment: The segment to do the extraction on.

    :type segment: str

    :param min_wl: The minimum wavlength value to keep.  If None, then there will be no lower bound.

    :type min_wl: float

    :param max_wl: The maximum wavelength value to keep.  If None, then there will be no upper bound.

    :type max_wl: float
    """
    if segment in cos_spectrum.segments:
        if min_wl is not None and max_wl is not None:
            where_within = numpy.where( (cos_spectrum.segments[segment].wavelengths >= min_wl) & (cos_spectrum.segments[segment].wavelengths <= max_wl) )[0]
        elif min_wl is not None and max_wl is None:
            where_within = numpy.where( (cos_spectrum.segments[segment].wavelengths >= min_wl) )[0]
        elif min_wl is None and max_wl is not None:
            where_within = numpy.where( (cos_spectrum.segments[segment].wavelengths <= max_wl) )[0]
        
        n_within = len(where_within)
        if n_within > 0:
            """ Then we extract the subspectrum for this object by modifying the cos_spectrum object. """
            cos_spectrum.segments[segment].nelem = n_within
            cos_spectrum.segments[segment].wavelengths = cos_spectrum.segments[segment].wavelengths[where_within]
            cos_spectrum.segments[segment].fluxes = cos_spectrum.segments[segment].fluxes[where_within]
            cos_spectrum.segments[segment].fluxerrs = cos_spectrum.segments[segment].fluxerrs[where_within]
            cos_spectrum.segments[segment].dqs = cos_spectrum.segments[segment].dqs[where_within]
        else:
            print "*** WARNING in SPECUTILS_COS: Requested subspectrum does not overlap with this segment's spectrum.  No extraction will be done."
    else:
        raise SpecUtilsError("The segment where you want to perform the subspectrum extraction is not present.  Specified \"" + segment + "\", available segments for this spectrum are: " + ', '.join(cos_spectrum.segments)+".")

#--------------------
