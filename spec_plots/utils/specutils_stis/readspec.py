"""
.. module:: readspec
   :synopsis: Reads in a STIS spectrum from a FITS file.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

#--------------------
# Built-In Imports
#--------------------
from __future__ import absolute_import
from builtins import range
#--------------------
# External Imports
#--------------------
from astropy.io import fits
#--------------------
# Package Imports
#--------------------
from spec_plots.utils.specutils_stis.stis1dspectrum import (
    STIS1DSpectrum, STISExposureSpectrum, STISOrderSpectrum)
from spec_plots import __version__

#--------------------

#--------------------
def readspec(input_file):
    """
    Reads in a STIS spectrum FITS file (x1d, sx1) and returns the wavelengths,
    fluxes, and flux uncertainties for the two (FUV segments) or three (NUV
    stripes).

    :param input_file: Name of input FITS file.

    :type input_file: str

    :returns: STIS1DSpectrum -- The spectroscopic data (wavelength, flux, flux
        error, etc):
    """

    with fits.open(input_file) as hdulist:
        # Create an initially empty list that will contain each extension's
        # (association's) spectrum object.
        all_association_spectra = []

        # Loop over each association and create the COS spectrum objects.
        for exten in hdulist[1:]:
            exten_data_table = exten.data

            # How many orders (table rows) in this extension?
            n_orders = len(exten_data_table["sporder"])

            # Create a list of STISOrderSpectra for this extension.
            all_order_spectra = [STISOrderSpectrum(
                nelem=exten_data_table["nelem"][order],
                wavelengths=exten_data_table["WAVELENGTH"][order],
                fluxes=exten_data_table["FLUX"][order],
                fluxerrs=exten_data_table["ERROR"][order],
                dqs=exten_data_table["DQ"][order]) for order in range(n_orders)]

            # Create a STISExposureSpectrum from the STISOrderSpectrum
            # objects.  Append to the running list of them.
            this_exposure_spectrum = STISExposureSpectrum(all_order_spectra)
            all_association_spectra.append(this_exposure_spectrum)

        return STIS1DSpectrum(all_association_spectra, orig_file=input_file)
#--------------------
