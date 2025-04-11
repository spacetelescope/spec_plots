"""
.. module:: readspec
   :synopsis: Reads in a JWST spectrum from a FITS file.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

#--------------------
# Built-In Imports
#--------------------
from __future__ import absolute_import
from __future__ import print_function
import sys
#--------------------
# External Imports
#--------------------
from astropy.io import fits
from astropy.table import Table
#--------------------
# Package Imports
#--------------------
from spec_plots import __version__
from spec_plots.utils.specutils_jwst.jwstspectrum import JWSTSpectrum

#--------------------

#--------------------
def readspec(input_file):
    """
    Reads in a JWST spectrum FITS file (x1d, x1dints) and
    returns the wavelengths, fluxes, flux uncertainties, and DQ values.

    :param input_file: Name of input FITS file.

    :type input_file: str

    :returns: JWSTSpectrum -- The spectroscopic data (wavelength, flux,
        flux error, etc):

    :raises: KeyError
    """

    with fits.open(input_file) as hdulist:
        # Read the data from the "EXTRACT1D" FITS extension.
        try:
            jwst_table = Table.read(hdulist["EXTRACT1D"],
                                        unit_parse_strict='silent')
        except KeyError:
            print("*** MAKE_JWST_SPEC_PREVIEWS ERROR: EXTRACT1D extension not"
                  " found in FITS file.")
            sys.exit()

        # Extract wavelength, fluxes, flux uncertainties, and DQ flags for
        # each segment.
        try:
            wavelength_col = jwst_table["WAVELENGTH"]
        except KeyError:
            print("*** MAKE_JWST_SPEC_PREVIEWS ERROR: WAVELENGTH column not"
                  " found in first extension's binary table.")
            sys.exit()
        wavelength_unit = wavelength_col.unit

        try:
            flux_col = jwst_table["FLUX"]
        except KeyError:
            print("*** MAKE_JWST_SPEC_PREVIEWS ERROR: FLUX column not found in"
                  " first extension's binary table.")
            sys.exit()
        flux_unit = flux_col.unit

        try:
            fluxerr_col = jwst_table["ERROR"]
        except KeyError:
            try:
                fluxerr_col = jwst_table["FLUX_ERROR"]
            except KeyError:
                print("*** MAKE_JWST_SPEC_PREVIEWS ERROR: neither ERROR "
                          "nor FLUX_ERROR column found in first "
                          "extension's binary table.")
                sys.exit()

        try:
            dq_col = jwst_table["DQ"]
        except KeyError:
            print("*** MAKE_JWST_SPEC_PREVIEWS ERROR: DQ column not found"
                  " in first extension's binary table.")
            sys.exit()

        # Create JWSTSpectrum object.
        return_spec = JWSTSpectrum(wavelength_col.data, flux_col.data,
                                       fluxerr_col.data, dq_col.data,
                                       wavelength_unit, flux_unit,
                                       orig_file=input_file)

        return return_spec
#--------------------
