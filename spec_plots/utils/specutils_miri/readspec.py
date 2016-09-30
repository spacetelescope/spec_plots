"""
.. module:: readspec
   :synopsis: Reads in a MIRI spectrum from a FITS file.

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
from spec_plots.utils.specutils_miri.mirispectrum import MIRISpectrum

#--------------------

#--------------------
def readspec(input_file):
    """
    Reads in a MIRI spectrum FITS file (x1d, x1dints) and
    returns the wavelengths, fluxes, flux uncertainties, and DQ values.

    :param input_file: Name of input FITS file.

    :type input_file: str

    :returns: MIRISpectrum -- The spectroscopic data (wavelength, flux,
        flux error, etc):

    :raises: KeyError
    """

    with fits.open(input_file) as hdulist:

        # Read the data from the first extension.  For MIRI, the spectra are
        # always stored as tables in the first FITS extension.
        miri_tabledata = hdulist[1].data

        # Extract wavelength, fluxes, flux uncertainties, and DQ flags for
        # each segment.
        try:
            wavelength_table = miri_tabledata.field("WAVELENGTH")
        except KeyError:
            print("*** MAKE_JWST_SPEC_PREVIEWS ERROR: WAVELENGTH column not"
                  " found in first extension's binary table.")
            exit(1)

        try:
            flux_table = miri_tabledata.field("FLUX")
        except KeyError:
            print("*** MAKE_JWST_SPEC_PREVIEWS ERROR: FLUX column not found in"
                  " first extension's binary table.")
            exit(1)

        try:
            fluxerr_table = miri_tabledata.field("ERR")
        except KeyError:
            print("*** MAKE_JWST_SPEC_PREVIEWS ERROR: ERR column not found"
                  " in first extension's binary table.")
            exit(1)

        try:
            dq_table = miri_tabledata.field("DQ")
        except KeyError:
            print("*** MAKE_JWST_SPEC_PREVIEWS ERROR: DQ column not found"
                  " in first extension's binary table.")
            exit(1)

        # Create MIRISpectrum object.
        return_spec = MIRISpectrum(wavelength_table, flux_table, fluxerr_table,
                                   dq_table, orig_file=input_file)

        return return_spec
#--------------------
