"""
.. module:: readspec
   :synopsis: Reads in a HASP spectrum from a FITS file.

.. moduleauthor:: Rob Swaters <rswaters@stsci.edu>
"""

#--------------------
# Built-In Imports
#--------------------
from __future__ import absolute_import
#--------------------
# External Imports
#--------------------
from astropy.io import fits
#--------------------
# Package Imports
#--------------------
from spec_plots.utils.specutils_hasp.haspspectrum import HASPSpectrum
from spec_plots import __version__

#--------------------

#--------------------
def readspec(input_file):
    """
    Reads in a HASP spectrum FITS file (*_cspec.fits) and returns the
    wavelengths, fluxes, and flux uncertainties.

    :param input_file: Name of input FITS file.

    :type input_file: str

    :returns: HASPSpectrum -- The spectroscopic data (wavelength, flux, flux
        error, etc):
    """

    with fits.open(input_file) as hdulist:
        # Create an initially empty list that will contain each extension's
        # (association's) spectrum object.

        try:
            title = hdulist[0].header['TARGNAME']
        except KeyError:
            title = ""

        plot_info = {"title": title, "orig_file": input_file}

        # Data are in the first extension.
        data_table = hdulist[1].data

        return HASPSpectrum(
            nelem=len(data_table["WAVELENGTH"][0]),
            wavelengths=data_table["WAVELENGTH"][0],
            fluxes=data_table["FLUX"][0],
            fluxerrs=data_table["ERROR"][0],
            dqs=[0 if x > 0 else 1 for x in data_table["EFF_EXPTIME"][0]],
            plot_info=plot_info
            )

#--------------------
