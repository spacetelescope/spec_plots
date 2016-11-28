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
#--------------------
# External Imports
#--------------------
from astropy.io import fits
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

        # Read the data from the first extension.  For JWST, the spectra are
        # always stored as tables in the first FITS extension.
        jwst_tabledata = hdulist[1].data

        # Extract wavelength, fluxes, flux uncertainties, and DQ flags for
        # each segment.
        try:
            wavelength_table = jwst_tabledata.field("WAVELENGTH")
        except KeyError:
            print("*** MAKE_JWST_SPEC_PREVIEWS ERROR: WAVELENGTH column not"
                  " found in first extension's binary table.")
            exit(1)

        try:
            flux_table = jwst_tabledata.field("FLUX")
        except KeyError:
            print("*** MAKE_JWST_SPEC_PREVIEWS ERROR: FLUX column not found in"
                  " first extension's binary table.")
            exit(1)

        try:
            fluxerr_table = jwst_tabledata.field("ERROR")
        except KeyError:
            print("*** MAKE_JWST_SPEC_PREVIEWS ERROR: ERROR column not found"
                  " in first extension's binary table.")
            exit(1)

        try:
            dq_table = jwst_tabledata.field("DQ")
        except KeyError:
            print("*** MAKE_JWST_SPEC_PREVIEWS ERROR: DQ column not found"
                  " in first extension's binary table.")
            exit(1)

        # Create JWSTSpectrum object.
        return_spec = JWSTSpectrum(wavelength_table, flux_table, fluxerr_table,
                                   dq_table, orig_file=input_file)

        return return_spec
#--------------------
