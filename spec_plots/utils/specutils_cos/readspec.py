"""
.. module:: readspec
   :synopsis: Reads in a COS spectrum from a FITS file.

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
import numpy
#--------------------
# Package Imports
#--------------------
from spec_plots import __version__
from spec_plots.utils.specutils_cos.check_segments import check_segments
from spec_plots.utils.specutils_cos.cosspectrum import COSSpectrum, COSSegment

#--------------------

#--------------------
def readspec(input_file):
    """
    Reads in a COS spectrum FITS file (x1d, x1dsum, or x1dsum{1,2,3,4}) and
    returns the wavelengths, fluxes, flux uncertainties, and DQ_WGT values for
    the two (FUV segments) or three (NUV stripes).

    :param input_file: Name of input FITS file.

    :type input_file: str

    :returns: COSSpectrum -- The spectroscopic data (wavelength, flux,
        flux error, etc):

    :raises: KeyError
    """

    with fits.open(input_file) as hdulist:

        # Read the data from the first extension.  For COS, the spectra are
        # always stored as tables in the first FITS extension.
        cos_tabledata = hdulist[1].data

        # Extract the SEGMENTS.  This is either an (up-to) 2-element array
        # of ["FUVA", "FUVB"], or an (up-to) 3-element array of ["NUVA", "NUVB",
        # "NUVC"].
        try:
            segment_arr = cos_tabledata.field("SEGMENT")
        except KeyError:
            print("*** MAKE_HST_SPEC_PREVIEWS ERROR: SEGMENT column not found"
                  " in first extension's binary table.")
            exit(1)

        # Determine which band this is (NUV, FUV).
        band = check_segments(segment_arr, input_file)

        # Extract the optical element from the primary header.
        try:
            optical_element = hdulist[0].header["OPT_ELEM"]
        except KeyError:
            print("*** MAKE_HST_SPEC_PREVIEWS ERROR: OPT_ELEM keyword not"
                  " found in the primary header.")
            exit(1)

        # Extract the number of elements (n_wavelengths, n_fluxes, etc.) for
        # each segment.  The dimension will match the array of segment names.
        try:
            nelems_arr = cos_tabledata.field("NELEM")
        except KeyError:
            print("*** MAKE_HST_SPEC_PREVIEWS ERROR: NELEM column not found"
                  " in first extension's binary table.")
            exit(1)

        # Extract wavelength, fluxes, flux uncertainties, and DQ flags for
        # each segment.  These will be mxn tables, where m is the number of
        # segment names, and n is the number of elements from the nelems array.
        try:
            wavelength_table = cos_tabledata.field("WAVELENGTH")
        except KeyError:
            print("*** MAKE_HST_SPEC_PREVIEWS ERROR: WAVELENGTH column not"
                  " found in first extension's binary table.")
            exit(1)

        try:
            flux_table = cos_tabledata.field("FLUX")
        except KeyError:
            print("*** MAKE_HST_SPEC_PREVIEWS ERROR: FLUX column not found in"
                  " first extension's binary table.")
            exit(1)

        try:
            fluxerr_table = cos_tabledata.field("ERROR")
        except KeyError:
            print("*** MAKE_HST_SPEC_PREVIEWS ERROR: ERROR column not found"
                  " in first extension's binary table.")
            exit(1)

        try:
            dq_table = cos_tabledata.field("DQ_WGT")
        except KeyError:
            print("*** MAKE_HST_SPEC_PREVIEWS ERROR: DQ_WGT column not found"
                  " in first extension's binary table.")
            exit(1)

        # Create COSSegment objects to populate the COSSpectrum object with.
        if band == 'FUV':
            # Try creating FUVA COSSegment object.
            try:
                fuva_index = numpy.where(segment_arr == "FUVA")[0][0]
            except IndexError:
                # Then there is no FUVA segment, normally this happends if
                # there is only one segment present (in which case it's the
                # other one of the two).
                fuva_index = None
            if fuva_index is not None:
                fuva_cossegment = COSSegment(nelem=nelems_arr[fuva_index],
                                             wavelengths=
                                             wavelength_table[fuva_index, :],
                                             fluxes=flux_table[fuva_index, :],
                                             fluxerrs=fluxerr_table[
                                                 fuva_index, :],
                                             dqs=dq_table[fuva_index, :])

            # Try creating FUVB COSSegment object.
            try:
                fuvb_index = numpy.where(segment_arr == "FUVB")[0][0]
            except IndexError:
                # Then there is no FUVB segment, normally this happens if
                # there is only one segment present (in which case it's the
                # other one of the two).
                fuvb_index = None
            if fuvb_index is not None:
                fuvb_cossegment = COSSegment(nelem=nelems_arr[fuvb_index],
                                             wavelengths=wavelength_table[
                                                 fuvb_index, :],
                                             fluxes=flux_table[fuvb_index, :],
                                             fluxerrs=fluxerr_table[
                                                 fuvb_index, :],
                                             dqs=dq_table[fuvb_index, :])
        elif band == 'NUV':
            # Try creating NUVA COSSegment object.
            try:
                nuva_index = numpy.where(segment_arr == "NUVA")[0][0]
            except IndexError:
                # Then there is no NUVA segment, not sure this is possible
                # but handle the case anyways.
                nuva_index = None
            if nuva_index is not None:
                nuva_cossegment = COSSegment(nelem=nelems_arr[nuva_index],
                                             wavelengths=wavelength_table[
                                                 nuva_index, :],
                                             fluxes=flux_table[nuva_index, :],
                                             fluxerrs=fluxerr_table[
                                                 nuva_index, :],
                                             dqs=dq_table[nuva_index, :])

            # Try creating NUVB COSSegment object.
            try:
                nuvb_index = numpy.where(segment_arr == "NUVB")[0][0]
            except IndexError:
                # Then there is no NUVB segment, not sure this is possible
                # but handle the case anyways.
                nuvb_index = None
            if nuvb_index is not None:
                nuvb_cossegment = COSSegment(nelem=nelems_arr[nuvb_index],
                                             wavelengths=wavelength_table[
                                                 nuvb_index, :],
                                             fluxes=flux_table[nuvb_index, :],
                                             fluxerrs=fluxerr_table[
                                                 nuvb_index, :],
                                             dqs=dq_table[nuvb_index, :])

            # Try creating NUVC COSSegment object.
            try:
                nuvc_index = numpy.where(segment_arr == "NUVC")[0][0]
            except IndexError:
                # Then there is no NUVC segment, not sure this is possible
                # but handle the case anyways.
                nuvc_index = None
            if nuvc_index is not None:
                nuvc_cossegment = COSSegment(nelem=nelems_arr[nuvc_index],
                                             wavelengths=wavelength_table[
                                                 nuvc_index, :],
                                             fluxes=flux_table[nuvc_index, :],
                                             fluxerrs=fluxerr_table[
                                                 nuvc_index, :],
                                             dqs=dq_table[nuvc_index, :])

        # Create COSSpectrum object.
        if band == 'FUV':
            # Handle case where both are supplied.
            if fuva_index is not None and fuvb_index is not None:
                return_spec = COSSpectrum(optical_element, band=band,
                                          cos_segments={'FUVA':fuva_cossegment,
                                                        'FUVB':fuvb_cossegment},
                                          orig_file=input_file)
            elif fuva_index is not None:
                # Handle cases where only one is supplied.
                return_spec = COSSpectrum(optical_element, band=band,
                                          cos_segments={'FUVA':fuva_cossegment},
                                          orig_file=input_file)
            elif fuvb_index is not None:
                return_spec = COSSpectrum(optical_element, band=band,
                                          cos_segments={'FUVB':fuvb_cossegment},
                                          orig_file=input_file)
            else:
                raise ValueError("Neither FUVA or FUVB segments were found,"
                                 " unable to create COS spectrum object.")

        elif band == 'NUV':
            # Handle case where all three are supplied.
            if (nuva_index is not None and nuvb_index is not None and
                    nuvc_index is not None):
                return_spec = COSSpectrum(optical_element, band=band,
                                          cos_segments={'NUVA':nuva_cossegment,
                                                        'NUVB':nuvb_cossegment,
                                                        'NUVC':nuvc_cossegment},
                                          orig_file=input_file)
            elif nuva_index is not None and nuvb_index is not None:
                # Handle cases where only two are supplied.
                return_spec = COSSpectrum(optical_element, band=band,
                                          cos_segments={'NUVA':nuva_cossegment,
                                                        'NUVB':nuvb_cossegment},
                                          orig_file=input_file)
            elif nuva_index is not None and nuvc_index is not None:
                return_spec = COSSpectrum(optical_element, band=band,
                                          cos_segments={'NUVA':nuva_cossegment,
                                                        'NUVC':nuvc_cossegment},
                                          orig_file=input_file)
            elif nuvb_index is not None and nuvc_index is not None:
                return_spec = COSSpectrum(optical_element, band=band,
                                          cos_segments={'NUVB':nuvb_cossegment,
                                                        'NUVC':nuvc_cossegment},
                                          orig_file=input_file)
            elif nuva_index is not None:
                # Handle cases where only one is supplied.
                return_spec = COSSpectrum(optical_element, band=band,
                                          cos_segments={'NUVA':nuva_cossegment},
                                          orig_file=input_file)
            elif nuvb_index is not None:
                return_spec = COSSpectrum(optical_element, band=band,
                                          cos_segments={'NUVB':nuvb_cossegment},
                                          orig_file=input_file)
            elif nuvc_index is not None:
                return_spec = COSSpectrum(optical_element, band=band,
                                          cos_segments={'NUVC':nuvc_cossegment},
                                          orig_file=input_file)
            else:
                raise ValueError("None of the NUVA, NUVB, or NUVC segments"
                                 " were found, unable to create COS spectrum"
                                 " object.")

        return return_spec
#--------------------
