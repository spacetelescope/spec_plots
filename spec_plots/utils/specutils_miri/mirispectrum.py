"""
.. module:: mirispectrum
   :synopsis: Class definitions for a MIRI Spectrum object.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

#--------------------
# Built-In Imports
#--------------------
from __future__ import absolute_import
from builtins import object
#--------------------
# Package Imports
#--------------------
from spec_plots import __version__

#--------------------

class MIRISpectrum(object):
    """
    Defines a MIRI spectrum, including wavelegnth, flux, flux errors, and DQ
    values.

    :raises: ValueError
    """
    def __init__(self, wl_arr, fl_arr, fle_arr, dq_arr, orig_file=None):
        """
        Create a MIRISpectrum object.

        :param wl_arr: The array of wavelength values.

        :type wl_arr: numpy.ndarray

        :param fl_arr: The array of flux values.

        :type fl_arr: numpy.ndarray

        :param fle_arr: The array of flux uncertainty values.

        :type fle_arr: numpy.ndarray

        :param dq_arr: The array of DQ values.

        :type dq_arr: numpy.ndarray

        :param orig_file: Original FITS file read to create the spectrum
            (includes full path).

        :type orig_file: str
        """

        # Record the original file name.
        self.orig_file = orig_file

        # Record the wavelength, flux, flux uncertainty, and DQ values.
        self.wavelengths = wl_arr
        self.fluxes = fl_arr
        self.fluxerrs = fle_arr
        self.dqs = dq_arr
#--------------------
