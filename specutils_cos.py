############################################################################################
## SPECUTILS_COS
##
## Given an HST COS spectrum file name, this script will read in data and create preview plots of various sizes and formats.
##
############################################################################################



############################################################################################
## Place import commands and logging options.
############################################################################################
import os
import pyfits
import numpy
import matplotlib.pyplot as pyplot
from make_hst_spec_previews import HSTSpecPrevError
############################################################################################



############################################################################################
## This function checks that the array of "segments" in the COS spectrum header are expected values.  It returns a scalar string representing the band (either "FUV" or "NUV").  If the array is not what's expected for COS, an Exception is raised.
############################################################################################
def check_segments(segments_arr, input_file):
    arr_len = len(segments_arr)
    ## Segment array must contain either two or three elements.
    if arr_len == 2:
        ## Must be ["FUVA", "FUVB"]
        if numpy.array_equal(segments_arr, ["FUVA", "FUVB"]):
            this_band = "FUV"
        else:
            raise HSTSpecPrevError("The array of SEGMENT strings contains 2 values, but is not equal to [\"FUVA\", \"FUVB\"] in file " + input_file)
    elif arr_len == 3:
        ## Must be ["NUVA", "NUVB", "NUVC"]
        if numpy.array_equal(segments_arr, ["NUVA", "NUVB", "NUVC"]):
            this_band = "NUV"
        else:
            raise HSTSpecPrevError("The array of SEGMENT strings contains 3 values, but is not equal to [\"NUVA\", \"NUVB\", \"NUVC\"] in file " + input_file)
    else:
        raise HSTSpecPrevError("The array of SEGMENT strings should contain 2 or 3 values, found " + str(arr_len) + " in file " + input_file)
    return this_band
############################################################################################



############################################################################################
## This function reads in COS spectrum FITS file (x1d, x1dsum, or x1dsum{1,2,34} FITS files) and returns the wavelengths, fluxes, and flux uncertainties for the two (FUV segments) or three (NUV stripes).
##
## The return value is a tuple containing {band, wavelength table, flux table, flux uncertainty table}, where:
##      band = Scalar string, either "FUV" or "NUV".
##      wavelength table = Table of wavelengths for each segment/slice, where a given segment/slice is a row in the table.  FUV has 2 segments, NUV has 3 slices.
##      flux table = Table of fluxes for each segment/slice, where a given segment/slice is a row in the table.  FUV has 2 segments, NUV has 3 slices.
##      flux uncertainty table = Table of flux uncertainties for each segment/slice, where a given segment/slice is a row in the table.  FUV has 2 segments, NUV has 3 slices.
############################################################################################
def readspec(input_file, verbosity):
    hdulist = pyfits.open(input_file)

    ## Read the data from the first extension.  For COS, the spectra are always stored as tables in the first FITS extension.
    cos_tabledata = hdulist[1].data

    ## Extract the SEGMENTS.  This is either a 2-element array of ["FUVA", "FUVB"], or a 3-element array of ["NUVA", "NUVB", "NUVC"].
    segment_arr = cos_tabledata.field("SEGMENT")
    band = check_segments(segment_arr, input_file)

    ## Extract the number of elements (n_wavelengths, n_fluxes, etc.) for each segment.  This will also be either a 2-element array (FUV) or 3-element array (NUV).
    nelems_arr = cos_tabledata.field("NELEM")

    ## Extract wavelength, fluxes, and flux uncertainties for each segment.  These will be either 2xn (FUV) or 3xn (NUV) tables.
    wavelength_table = cos_tabledata.field("WAVELENGTH")
    flux_table = cos_tabledata.field("FLUX")
    fluxerr_table = cos_tabledata.field("ERROR")

    return (band, wavelength_table, flux_table, fluxerr_table)
############################################################################################



############################################################################################
## This function accepts a COS spectrum tuple from the READSPEC function and produces preview plots.
############################################################################################
def plotspec(this_spec_tuple, output_type, output_file):
    ## Make sure the output path exists, if not, create it.
    if not os.path.exists(os.path.dirname(output_file)):
        os.mkdir(os.path.dirname(output_file))
    ## Extract the band from the tuple.
    this_band = this_spec_tuple[0]
    ## Start plot figure.
    pyplot.figure(1)
    if this_band == "FUV":
        ## Then there are two segments, so make a 2-panel plot.
        ## First plot area (top) is for *FUVB*.
        segment_1_plotarea = pyplot.subplot(2, 1, 1)
        pyplot.plot(this_spec_tuple[1][0], this_spec_tuple[2][0], 'k')
        ## Second plot area (bottom) is for *FUVA*.
        segment_2_plotarea = pyplot.subplot(2, 1, 2)
        pyplot.plot(this_spec_tuple[1][1], this_spec_tuple[2][1], 'k')
        if output_type == "png":
            pyplot.savefig(output_file, format=output_type)
        elif output_type == "screen":
            pyplot.show()
    elif this_band == "NUV":
        pass
    else:
        raise HSTSpecPrevError("The band is not understood, expect either \"FUV\" or \"NUV\", found " + this_band + " in file " + input_file)
############################################################################################
