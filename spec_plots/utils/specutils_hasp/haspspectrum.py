"""
.. module:: haspspectrum
   :synopsis: Defines the class for HASP spectra.

.. moduleauthor:: Rob Swaters <rswaters@stsci.edu>
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

class HASPSpectrum():
    """
    Defines a HASP 1D spectrum, including wavelegnth, flux, and flux errors.

    :raises: ValueError
    """

    def __init__(self, nelem=None, wavelengths=None, fluxes=None, fluxerrs=None,
                 dqs=None, plot_info=None, orig_file=None):
        """
        Create a HASPSpectrum object, default to empty values.  Allows user
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

        :param plot_info: dictionary with information for the preview plot.

        :type plot_info: dict

        :param orig_file: Original FITS file read to create the spectrum
            (includes full path).

        :type orig_file: str
        """

        # Record the original file name along with the list of associations.
        self.orig_file = orig_file

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

        if plot_info is not None:
            self.plot_info = plot_info
        else:
            self.plot_info = {}

#--------------------
