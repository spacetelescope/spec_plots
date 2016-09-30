"""
.. module:: stis1dspectrum
   :synopsis: Defines the class for STIS 1D spectra.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

#--------------------
# Built-In Imports
#--------------------
from __future__ import absolute_import
from builtins import object
#--------------------
# External Imports
#--------------------
import numpy
#--------------------
# Package Imports
#--------------------
from spec_plots import __version__

#--------------------

class STIS1DSpectrum(object):
    """
    Defines a STIS 1D spectrum (either "x1d" extracted or "sx1" summed
    extracted), including wavelegnth, flux, and flux errors.  A STIS 1D
    spectrum object consists of N associations.  If the file is an association,
    then N > 1, otherwise N = 1.  Each of these N associations can contain M
    orders.  If the association is an Echelle spectrum, then 24 < M < 70,
    depending on instrument configuration, otherwise M = 1.  Each of these M
    orders contain typical spectral data (wavelengths, fluxes, etc.), stored
    as STISOrderSpectrum objects.  The final data structure is then
    <STIS1DSpectrum>.associations[n].order[m].wavelengths (or .fluxes,
    .fluxerrs, etc.).

    :raises: ValueError
    """

    def __init__(self, association_spectra, orig_file=None):
        """
        Create a STIS1DSpectrum object out of a list of STISExposureSpectrum
        objects, which themselves are lists of STISOrderSpectrum objects.

        :param association_spectra: A list whose length is equal to the number
            of associations (length = "N" associations).  Each element in this
            list is a STISExposureSpectrum object, which itself is a list
            (length = "M" orders) of STISOrderSpectrum objects.

        :type association_spectra: list

        :param orig_file: Original FITS file read to create the spectrum
            (includes full path).

        :type orig_file: str
        """

        # Record the original file name along with the list of associations.
        self.orig_file = orig_file
        self.associations = association_spectra
#--------------------

#--------------------
class STISExposureSpectrum(object):
    """
    Defines a STIS exposure spectrum, which consists of "M" STISOrderSpectrum
    objects.
    """
    def __init__(self, order_spectra):
        """
        Create a STISExposureSpectrum object out of a list of STISOrderSpectrum
        objects.

        :param order_spectra: The STISOrderSpectrum objects to build the
            STISExposureSpectrum object out of.

        :type order_spectra: list

        :raises: ValueError
        """

        if len(order_spectra) > 0:
            self.orders = order_spectra
        else:
            raise ValueError("Must provide a list of at least one"
                             " STISOrderSpectrum object, input list is empty.")
#--------------------

#--------------------
class STISOrderSpectrum(object):
    """
    Defines a STIS order spectrum, including wavelength, flux, flux errors,
    and data quality flags, which are stored as numpy arrays.  A scalar int
    property provides the number of elements in this segment.
    """

    def __init__(self, nelem=None, wavelengths=None, fluxes=None, fluxerrs=None,
                 dqs=None):

        """
        Create a STISOrderSpectrum object, default to empty values.  Allows user
        to preallocate space if they desire by setting "nelem" but not providing
        lists/arrays on input right away.

        :param nelem: Number of elements for this segment's spectrum.

        :type nelem: int

        :param wavelengths: List of wavelengths in this segment's spectrum.

        :type wavelengths: list

        :param fluxes: List of fluxes in this segment's spectrum.

        :type fluxes: list

        :param fluxerrs: List of flux uncertainties in this segment's spectrum.

        :type fluxerrs: list

        :param dqs: List of data quality flags.

        :type dqs: list
        """

        # <DEVEL> Should it be required to have `nelem` > 0 *OR* specify
        # arrays on input?  Otherwise they are pre-allocated to empty
        # lists. </DEVEL>
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
