"""
.. module:: cosspectrum
   :synopsis: Class definitions for a COS Spectrum object.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

#--------------------
# Built-In Imports
#--------------------
from __future__ import absolute_import
from builtins import object
from builtins import str
#--------------------
# External Imports
#--------------------
import numpy
#--------------------
# Package Imports
#--------------------
from spec_plots import __version__

#--------------------

class COSSpectrum(object):
    """
    Defines a COS spectrum, including wavelegnth, flux, flux errors, and DQ_WGT
    values.  A COS spectrum consists of N segments (N = {2,3}) stored as a dict
    object.  Each of these dicts contain a COSSegment object that contains the
    wavelengths, fluxes, flux errors, etc.

    :raises: ValueError
    """
    def __init__(self, optical_element, band=None, cos_segments=None,
                 orig_file=None):
        """
        Create a COSSpectrum object given a band choice (must be "FUV" or
        "NUV").

        :param optical_element: The string representation of the optical
            element used for this observation, e.g., "G140L".

        :type optical_element: str

        :param band: Which band is this spectrum for ("FUV" or "NUV")?

        :type band: str

        :param cos_segments: [Optional] COSSegment objects to populate the
            COSSpectrum with.

        :type cos_segments: dict

        :param orig_file: Original FITS file read to create the spectrum
            (includes full path).

        :type orig_file: str

        :raises: ValueError
        """

        # Record the optical element.
        self.optical_element = optical_element

        # Record the original file name.
        self.orig_file = orig_file

        if band.strip().upper() == "FUV":
            # Record the band name.
            self.band = band

            if cos_segments is not None:
                if len(cos_segments) == 2:
                    self.segments = {'FUVA':cos_segments['FUVA'],
                                     'FUVB':cos_segments['FUVB']}
                elif len(cos_segments) == 1 and 'FUVA' in cos_segments:
                    self.segments = {'FUVA':cos_segments['FUVA']}
                elif len(cos_segments) == 1 and 'FUVB' in cos_segments:
                    self.segments = {'FUVB':cos_segments['FUVB']}
                else:
                    raise ValueError("Band is specified as " +
                                     band.strip().upper() +
                                     ", expected 1 or 2 COSSegment objects as"
                                     " a list but received " +
                                     str(len(cos_segments))+".")
            else:
                # Otherwise `cos_segments` was not supplied, so create with
                # an empty `segments` property.
                self.segments = {'FUVA':COSSegment(), 'FUVB':COSSegment()}

        elif band.strip().upper() == "NUV":
            # Record the band name.
            self.band = band

            if cos_segments is not None:
                if len(cos_segments) == 3:
                    self.segments = {'NUVA':cos_segments['NUVA'],
                                     'NUVB':cos_segments['NUVB'],
                                     'NUVC':cos_segments['NUVC']}
                elif len(cos_segments) == 2:
                    if 'NUVA' in cos_segments and 'NUVB' in cos_segments:
                        self.segments = {'NUVA':cos_segments['NUVA'],
                                         'NUVB':cos_segments['NUVB']}
                    if 'NUVA' in cos_segments and 'NUVC' in cos_segments:
                        self.segments = {'NUVA':cos_segments['NUVA'],
                                         'NUVC':cos_segments['NUVC']}
                    if 'NUVB' in cos_segments and 'NUVC' in cos_segments:
                        self.segments = {'NUVB':cos_segments['NUVB'],
                                         'NUVC':cos_segments['NUVC']}
                elif len(cos_segments) == 1:
                    if 'NUVA' in cos_segments:
                        self.segments = {'NUVA':cos_segments['NUVA']}
                    if 'NUVB' in cos_segments:
                        self.segments = {'NUVB':cos_segments['NUVB']}
                    if 'NUVC' in cos_segments:
                        self.segments = {'NUVC':cos_segments['NUVC']}
                else:
                    raise ValueError("Band is specified as " +
                                     band.strip().upper() +
                                     ", expected 1, 2, or 3 COSSegment objects"
                                     " as a list but received " +
                                     str(len(cos_segments))+".")
            else:
                # Otherwise `cos_segments` was not supplied, so create with
                # an empty `segments` property.
                self.segments = {'NUVA':COSSegment(), 'NUVB':COSSegment(),
                                 'NUVC':COSSegment()}

        else:
            raise ValueError("Must specify band=\"FUV\" or band=\"NUV\".")
#--------------------

#--------------------
class COSSegment(object):
    """
    Defines a spectrum from a COS segment.  The data (wavelength, flux, flux
    errors) are stored as numpy ndarrays.  A scalar int property provides the
    number of elements in this segment.
    """

    def __init__(self, nelem=None, wavelengths=None, fluxes=None, fluxerrs=None,
                 dqs=None):
        """
        Create a COSSegment object, default to empty values.  Allows user to
        preallocate space if they desire by setting "nelem" but not providing
        lists/arrays on input right away.

        :param nelem: Number of elements for this segment's spectrum.

        :type nelem: int

        :param wavelengths: [Optional] List of wavelengths in this segment's
        spectrum.

        :type wavelengths: list

        :param fluxes: [Optional] List of fluxes in this segment's spectrum.

        :type fluxes: list

        :param fluxerrs: [Optional] List of flux uncertainties in this segment's
        spectrum.

        :type fluxerrs: list

        :param dqs: [Optional] List of Data Quality (DQ) flags in this segment's
        spectrum.

        :type dqs: list
        """

        # <DEVEL> Should it be required to have `nelem` > 0 *OR* specify
        # arrays on input?  Otherwise they are pre-allocated to empty lists.
        # </DEVEL>
        if nelem is not None:
            self.nelem = nelem
        else:
            self.nelem = 0

        if wavelengths is not None:
            self.wavelengths = numpy.asarray(wavelengths)
        else:
            self.wavelengths = numpy.zeros(self.nelem)

        if fluxes is not None:
            self.fluxes = numpy.asarray(fluxes)
        else:
            self.fluxes = numpy.zeros(self.nelem)

        if fluxerrs is not None:
            self.fluxerrs = numpy.asarray(fluxerrs)
        else:
            self.fluxerrs = numpy.zeros(self.nelem)

        if dqs is not None:
            self.dqs = numpy.asarray(dqs)
        else:
            self.dqs = numpy.zeros(self.nelem)
#--------------------
