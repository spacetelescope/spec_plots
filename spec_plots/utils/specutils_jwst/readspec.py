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
            sys.exit("*** MAKE_JWST_SPEC_PREVIEWS ERROR: EXTRACT1D extension"
                         " not found in FITS file.")

        # For time series observations (TSO), many exposures are included in
        # a single EXTRACT1D. In this case, a preview plot is only made for
        # a single exposure (the first one). Remake a simplified Table that
        # contains the columns of interest as columns rather than lists in a
        # row entry.
        extract_one_row = False
        if len(jwst_table["WAVELENGTH"].shape) == 2:
            extract_one_row = True

        # Extract wavelength, fluxes, flux uncertainties, and DQ flags for
        # each segment.
        try:
            if not extract_one_row:
                wavelength_col = jwst_table["WAVELENGTH"].data
            else:
                wavelength_col = jwst_table[0]["WAVELENGTH"]
            wavelength_unit = jwst_table["WAVELENGTH"].unit
        except KeyError:
            sys.exit("*** MAKE_JWST_SPEC_PREVIEWS ERROR: WAVELENGTH column not"
                  " found in first extension's binary table.")

        try:
            if not extract_one_row:
                flux_col = jwst_table["FLUX"].data
            else:
                flux_col = jwst_table[0]["FLUX"]
            flux_unit = jwst_table["FLUX"].unit
        except KeyError:
            sys.exit("*** MAKE_JWST_SPEC_PREVIEWS ERROR: FLUX column not"
                         " found in first extension's binary table.")

        try:
            if not extract_one_row:
                fluxerr_col = jwst_table["ERROR"].data
            else:
                fluxerr_col = jwst_table[0]["ERROR"]
        except KeyError:
            try:
                if not extract_one_row:
                    fluxerr_col = jwst_table["FLUX_ERROR"].data
                else:
                    fluxerr_col = jwst_table[0]["FLUX_ERROR"]
            except KeyError:
                sys.exit("*** MAKE_JWST_SPEC_PREVIEWS ERROR: neither ERROR "
                          "nor FLUX_ERROR column found in first "
                          "extension's binary table.")

        try:
            if not extract_one_row:
                dq_col = jwst_table["DQ"].data
            else:
                dq_col = jwst_table[0]["DQ"]
        except KeyError:
            sys.exit("*** MAKE_JWST_SPEC_PREVIEWS ERROR: DQ column not found"
                  " in first extension's binary table.")

        # Create JWSTSpectrum object.
        return_spec = JWSTSpectrum(wavelength_col, flux_col,
                                       fluxerr_col, dq_col,
                                       wavelength_unit, flux_unit,
                                       orig_file=input_file)

        return return_spec
#--------------------
