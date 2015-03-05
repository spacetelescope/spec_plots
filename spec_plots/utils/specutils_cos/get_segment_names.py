__version__ = '1.32.0'

"""
.. module:: get_segment_names

   :synopsis: Returns a list of segment names such that the bluest segment comes first.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

#--------------------

def get_segment_names(cos_spectrum):
    """
    Returns a list of segment names sorted such that the bluest segment come first.

    :param cos_spectrum: COS spectrum as returned by READSPEC.

    :type cos_spectrum: COSSpectrum
    """

    """ Get an initial list of segment names. """
    segment_names = cos_spectrum.segments.keys()

    """ Reverse the list of segment names *FOR FUV DATA*, becaue the bluest segment is the latter in the alphabet, but only for the FUV spectra.  If NUV, then just make sure the segments are sorted alphabetically. """
    if cos_spectrum.band == 'FUV':
        segment_names.sort(reverse=True)
    else:
        segment_names.sort(reverse=False)

    return segment_names

#--------------------
