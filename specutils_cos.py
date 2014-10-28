__version__ = '1.2'

"""

.. module:: specutils_cos

   :synopsis: Contains functions for reading and plotting HST COS spectra.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

from astropy.io import fits
import math
from matplotlib import rc
import matplotlib.pyplot as pyplot
from matplotlib.ticker import FormatStrFormatter
import numpy
import os
import specutils

class COSSpectrum:
    """
    Defines a COS spectrum, including wavelegnth, flux, and flux errors.  A COS spectrum consists of N segments (N = {2,3}) stored as a dict object.  Each of these dicts contain a COSSegment object that contains the wavelengths, fluxes, flux errors, etc.

    :raises: ValueError
    """
    def __init__(self, band=None, cos_segments=None, orig_file=None):
        """
        Create a COSSpectrum object given a band choice (must be "FUV" or "NUV").

        :param band: Which band is this spectrum for ("FUV" or "NUV")?

        :type band: str

        :param cos_segments: [Optional] COSSegment objects to populate the COSSpectrum with.

        :type cos_segments: dict

        :param orig_file: Original FITS file read to create the spectrum (includes full path).

        :type orig_file: str
        """
        self.orig_file = orig_file
        if band.strip().upper() == "FUV":
            self.band = band
            if cos_segments is not None:
                if len(cos_segments) == 2:
                    self.segments = {'FUVA':cos_segments['FUVA'], 'FUVB':cos_segments['FUVB']}
                elif len(cos_segments) == 1 and 'FUVA' in cos_segments:
                    self.segments = {'FUVA':cos_segments['FUVA']}
                elif len(cos_segments) == 1 and 'FUVB' in cos_segments:
                    self.segments = {'FUVB':cos_segments['FUVB']}
                else:
                    raise ValueError("Band is specified as "+band.strip().upper()+", expected 1 or 2 COSSegment objects as a list but received "+str(len(cos_segments))+".")
            else:
                self.segments = {'FUVA':COSSegment(), 'FUVB':COSSegment()}
        elif band.strip().upper() == "NUV":
            self.band = band
            if cos_segments is not None:
                if len(cos_segments) == 3:
                    self.segments = {'NUVA':cos_segments['NUVA'], 'NUVB':cos_segments['NUVB'], 'NUVC':cos_segments['NUVC']}
                else:
                    raise ValueError("Band is specified as "+band.strip().upper()+", expected 3 COSSegment objects as a list but received "+str(len(cos_segments))+".")
            else:
                self.segments = {'NUVA':COSSegment(), 'NUVB':COSSegment(), 'NUVC':COSSegment()}
        else:
            raise ValueError("Must specify band=\"FUV\" or band=\"NUV\".")


class COSSegment:
    """
    Defines a spectrum from a COS segment.  The data (wavelength, flux, flux errors) are stored as numpy ndarrays.  A scalar int property provides the number of elements in this segment.
    """
    def __init__(self, nelem=None, wavelengths=None, fluxes=None, fluxerrs=None, dqs=None):
        """
        Create a COSSegment object, default to empty values.  Allows user to preallocate space if they desire by setting "nelem" but not providing lists/arrays on input right away.

        :param nelem: Number of elements for this segment's spectrum.

        :type nelem: int

        :param wavelengths: List of wavelengths in this segment's spectrum.

        :type wavelengths: list

        :param fluxes: List of fluxes in this segment's spectrum.

        :type fluxes: list

        :param fluxerrs: List of flux uncertainties in this segment's spectrum.

        :type fluxerrs: list

        :param dqs: List of Data Quality (DQ) flags in this segment's spectrum.

        :type dqs: list
        """
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


def check_segments(segments_list, input_file):
    """
    Checks that the array of "segments" in the COS spectrum header are expected values.  It returns a scalar string representing the band (either "FUV" or "NUV").  If the array is not what's expected for COS, an Exception is raised.

    :param segments_list: List of segment labels.

    :type segments_list: list

    :param input_file: Name of input FITS file.

    :type input_file: str

    :returns: str -- A string representation of the band, either "FUV" or "NUV".

    :raises: specutils.SpecUtilsError

    .. note::

         This function does not attempt to access the input file, it only requires the name of the file for error reporting purposes.
    """

    list_len = len(segments_list)
    """Segment array must contain either one, two, or three elements."""
    if list_len == 1:
        """Can happen with FUV data for some reason, so it can be either "FUVA" or "FUVB"."""
        if numpy.array_equal(segments_list, ["FUVA"]) or numpy.array_equal(segments_list, ["FUVB"]):
            this_band = "FUV"
        else:
            try:
                raise specutils.SpecUtilsError("The array of SEGMENT strings contains 1 value, but is not equal to [\"FUVA\"] or [\"FUVB\"] in file " + input_file)
            except specutils.SpecUtilsError as error_string:
                print error_string
                exit(1)
    elif list_len == 2:
        """Must be ["FUVA", "FUVB"]."""
        if numpy.array_equal(segments_list, ["FUVA", "FUVB"]):
            this_band = "FUV"
        else:
            try:
                raise specutils.SpecUtilsError("The array of SEGMENT strings contains 2 values, but is not equal to [\"FUVA\", \"FUVB\"] in file " + input_file)
            except specutils.SpecUtilsError as error_string:
                print error_string
                exit(1)
    elif list_len == 3:
        """Must be ["NUVA", "NUVB", "NUVC"]."""
        if numpy.array_equal(segments_list, ["NUVA", "NUVB", "NUVC"]):
            this_band = "NUV"
        else:
            try:
                raise specutils.SpecUtilsError("The array of SEGMENT strings contains 3 values, but is not equal to [\"NUVA\", \"NUVB\", \"NUVC\"] in file " + input_file)
            except specutils.SpecUtilsError as error_string:
                print error_string
                exit(1)
    else:
        try:
            raise specutils.SpecUtilsError("The array of SEGMENT strings should contain 1, 2, or 3 values, found " + str(list_len) + " in file " + input_file)
        except specutils.SpecUtilsError as error_string:
                print error_string
                exit(1)
    return this_band

def generate_cos_avoid_regions():
    """
    Creates a list of AvoidRegion objects for use in the plotting routine, specifically designed for COS spectra.
    """
    lya1215_ar = specutils.AvoidRegion(1214.,1217., "Lyman alpha emission line.")
    return [lya1215_ar]

def plotspec(cos_spectrum, output_type, output_file, n_consecutive, flux_scale_factor, fluxerr_scale_factor, output_size=None, debug=False, full_ylabels=False):
    """
    Accepts a COSSpectrum object from the READSPEC function and produces preview plots.

    :param cos_spectrum: COS spectrum as returned by READSPEC.

    :type cos_spectrum: COSSpectrum

    :param output_type: What kind of output to make?

    :type output_type: str

    :param output_file: Name of output file (including full path).

    :type output_file: str

    :param output_size: Size of plot in pixels (plots are square in dimensions).  Defaults to 1024.

    :param output_size: int

    :param full_ylabels: Should the y-labels contain the full values (including the power of 10 in scientific notation)?  Default = False.

    :type full_ylabels: bool

    :raises: OSError

    .. note::

         This function assumes a screen resolution of 96 DPI in order to generate plots of the desired sizes.  This is because matplotlib works in units of inches and DPI rather than pixels.
    """
    dpi_val = 96.
    if output_size is not None:
        if not isinstance(output_size, int):
            output_size = int(round(output_size))
    else:
        output_size = 1024
    """Make sure the output path exists, if not, create it."""
    if output_type != 'screen':
        if not os.path.isdir(os.path.dirname(output_file)):
            try:
                os.mkdir(os.path.dirname(output_file))
            except OSError as this_error:
                if this_error.errno == 13: 
                    print "*** MAKE_HST_SPEC_PREVIEWS ERROR: Output directory could not be created, "+repr(this_error.strerror)
                    exit(1)
                else:
                    raise
    """Extract the band type (FUV, NUV)."""
    this_band = cos_spectrum.band
    segment_names = cos_spectrum.segments.keys()
    """Reverse the list of segment names FOR FUV DATA, becaue the bluest segment is latter in the alphabet. but only for the FUV spectra.  If NUV, then just make sure the segments are sorted alphabetically."""
    if this_band == 'FUV':
        segment_names.sort(reverse=True)
    else:
        segment_names.sort(reverse=False)
    n_segments = len(segment_names)

    if output_size > 128:
        is_bigplot = True
        n_subplots = n_segments
        subplot_segment_names = segment_names
    else:
        is_bigplot = False
        """Change segment_names to be a one-element array, because we need to join all the segments into a single array if making a thumbnail plot.  Make sure it's still a list..."""
        n_subplots = 1
        subplot_segment_names = ["-".join(segment_names)]

    """Start plot figure."""
    this_figure, these_plotareas = pyplot.subplots(nrows=n_subplots, ncols=1, figsize=(output_size/dpi_val, output_size/dpi_val), dpi=dpi_val)

    if not isinstance(these_plotareas, numpy.ndarray):
        these_plotareas = numpy.asarray([these_plotareas])
    if is_bigplot:
        this_figure.subplots_adjust(hspace=0.3,top=0.915)
        this_figure.suptitle(os.path.basename(cos_spectrum.orig_file))
    else:
        this_figure.subplots_adjust(top=0.85,bottom=0.3,left=0.25,right=0.9)

    for i,s in enumerate(subplot_segment_names):
        this_plotarea = these_plotareas[i]
        if is_bigplot:
            all_wls = cos_spectrum.segments[s].wavelengths
            all_fls = cos_spectrum.segments[s].fluxes
            all_flerrs = cos_spectrum.segments[s].fluxerrs
            all_dqs = cos_spectrum.segments[s].dqs
            title_addendum = ""
        else:
            all_wls, all_fls, all_flerrs, all_dqs, title_addendum = specutils.stitch_components(cos_spectrum, n_consecutive, flux_scale_factor, fluxerr_scale_factor, segment_names=segment_names)

        median_flux, median_fluxerr, fluxerr_95th = specutils.get_flux_stats(all_fls, all_flerrs)

        if is_bigplot:
            this_plotarea.set_title(title_addendum, loc="right", size="small", color="red")

        """Determine optimal x-axis."""
        x_axis_range = specutils.set_plot_xrange("cos",all_wls, all_fls, all_flerrs, all_dqs, n_consecutive, flux_scale_factor, fluxerr_scale_factor, median_flux, median_fluxerr, fluxerr_95th)
        """Plot the spectrum, but only if valid wavelength ranges for x-axis are returned, otherwise plot a special "Fluxes Are All Zero" plot."""
        if all(numpy.isfinite(x_axis_range)):
            """Create COS avoid regions."""
            avoid_regions = generate_cos_avoid_regions()
            """Determine optimal y-axis, but only provide it with fluxes from the part of the spectrum that will be plotted based on the x-axis trimming."""
            y_axis_range = specutils.set_plot_yrange(all_wls, all_fls, avoid_regions=avoid_regions, wl_range=x_axis_range)
            this_plotarea.plot(all_wls, all_fls, 'b')
            this_plotarea.grid(True)
            """Get the values of the x- and y-axes plot limits that *would* have been used by pyplot automatically."""
            pyplot_xrange = this_plotarea.get_xlim()
            pyplot_yrange = this_plotarea.get_ylim()
            if debug:
                """Overplot points color-coded based on rejection criteria."""
                specutils.debug_oplot(this_plotarea, all_wls, all_fls, all_flerrs, all_dqs, median_flux, median_fluxerr, flux_scale_factor, fluxerr_scale_factor, fluxerr_95th)
                """Overplot the x-axis edges that are trimmed to define the y-axis plot range as a shaded area."""
                this_plotarea.axvspan(numpy.nanmin(all_wls), x_axis_range[0],facecolor="lightgrey",alpha=0.5)
                this_plotarea.axvspan(x_axis_range[1], numpy.nanmax(all_wls),facecolor="lightgrey",alpha=0.5)
                """Overplot the avoid regions in a light color as a shaded area."""
                for ar in avoid_regions:
                    this_plotarea.axvspan(ar.minwl,ar.maxwl,facecolor="lightgrey",alpha=0.5)
            """Note that we change the x-axis range here to be the min. and max. wavelength of this segment, rather than using the truncated version, so that all the plots for a similar instrument setting will have the same starting and ending plot values.  But, we still calculate the trimmed starting and ending wavelengths above for other things, such as defining the y-plot range."""
            this_plotarea.set_xlim(pyplot_xrange)
            """Change x-axis units to microns if a small plot, because there isn't enough space."""
            if not is_bigplot:
                rc('font', size=10)
                this_plotarea.set_xticklabels(this_plotarea.get_xticks()/10000.,rotation=45.)
                this_plotarea.locator_params(axis="both", nbins=4, steps=[1,2,4,6,8,10])
                this_figure.suptitle(r'$\lambda (\mu$m)', position=(0.83,0.99))
            else:
                """Make sure the font properties go back to normal."""
                pyplot.rcdefaults()
                this_plotarea.set_xlabel(r"Wavelength $(\AA)$")
                this_plotarea.set_ylabel(r"Flux $\mathrm{(erg/s/cm^2\!/\AA)}$")
                if full_ylabels:
                    this_plotarea.yaxis.set_major_formatter(FormatStrFormatter('%3.2E'))
            """Update y-axis range, but only adjust the ranges if this isn't an all-zero flux case (and not in debug mode, in which case I want to see the entire y-axis range)."""
            if not debug:
                this_plotarea.set_ylim(y_axis_range)
        else:
            x_axis_range = [numpy.nanmin(all_wls),numpy.nanmax(all_wls)]
            this_plotarea.set_xlim(x_axis_range)
            this_plotarea.set_axis_bgcolor("lightgrey")
            this_plotarea.set_yticklabels([])
            if not is_bigplot:
                rc('font', size=10)
                this_plotarea.set_xticklabels(this_plotarea.get_xticks()/10000.,rotation=45.)
                this_plotarea.locator_params(axis="both", nbins=4, steps=[1,2,4,6,8,10])
                this_figure.suptitle(r'$\lambda (\mu$m)', position=(0.83,0.99))
                textsize = "small"
                plottext = "Fluxes are \n all 0."
            else:
                """Make sure the font properties go back to normal."""
                pyplot.rcdefaults()
                this_plotarea.set_xlabel(r"Wavelength $(\AA)$")
                this_plotarea.set_ylabel(r"Flux $\mathrm{(erg/s/cm^2\!/\AA)}$")
                if full_ylabels:
                    this_plotarea.yaxis.set_major_formatter(FormatStrFormatter('%3.2E'))
                textsize = "x-large"
                plottext = "Fluxes are all 0."
            this_plotarea.text(0.5,0.5,plottext,horizontalalignment="center",verticalalignment="center",transform=this_plotarea.transAxes,size=textsize)

    """Display or plot to the desired format."""
    if output_type != "screen":
        """Deconstruct output file to include plot size information."""
        output_splits = os.path.split(output_file)
        file_splits = os.path.splitext(output_splits[1])
        revised_output_file = output_splits[0]+os.path.sep+file_splits[0]+'_{0:04d}'.format(output_size)+file_splits[1]
        this_figure.savefig(revised_output_file, format=output_type, dpi=dpi_val)
    elif output_type == "screen":
        pyplot.show()

def readspec(input_file):
    """
    Reads in a COS spectrum FITS file (x1d, x1dsum, or x1dsum{1,2,3,4}) and returns the wavelengths, fluxes, and flux uncertainties for the two (FUV segments) or three (NUV stripes).

    :param input_file: Name of input FITS file.

    :type input_file: str

    :returns: COSSpectrum -- The spectroscopic data (wavelength, flux, flux error, etc):

    :raises: KeyError
    """
    with fits.open(input_file) as hdulist:
        """Read the data from the first extension.  For COS, the spectra are always stored as tables in the first FITS extension."""
        cos_tabledata = hdulist[1].data
        """Extract the SEGMENTS.  This is either a 2-element array of ["FUVA", "FUVB"], or a 3-element array of ["NUVA", "NUVB", "NUVC"]."""
        try:
            segment_arr = cos_tabledata.field("SEGMENT")
        except KeyError:
            print "*** MAKE_HST_SPEC_PREVIEWS ERROR: SEGMENT column not found in first extension's binary table."
            exit(1)
        band = check_segments(segment_arr, input_file)
        """Extract the number of elements (n_wavelengths, n_fluxes, etc.) for each segment.  This will also be either a 2-element array (FUV) or 3-element array (NUV)."""
        try:
            nelems_arr = cos_tabledata.field("NELEM")
        except KeyError:
            print "*** MAKE_HST_SPEC_PREVIEWS ERROR: NELEM column not found in first extension's binary table."
            exit(1)
        """Extract wavelength, fluxes, and flux uncertainties for each segment.  These will be either 2xn (FUV) or 3xn (NUV) tables."""
        try:
            wavelength_table = cos_tabledata.field("WAVELENGTH")
        except KeyError:
            print "*** MAKE_HST_SPEC_PREVIEWS ERROR: WAVELENGTH column not found in first extension's binary table."
            exit(1)
        try:
            flux_table = cos_tabledata.field("FLUX")
        except KeyError:
            print "*** MAKE_HST_SPEC_PREVIEWS ERROR: FLUX column not found in first extension's binary table."
            exit(1)
        try:
            fluxerr_table = cos_tabledata.field("ERROR")
        except KeyError:
            print "*** MAKE_HST_SPEC_PREVIEWS ERROR: ERROR column not found in first extension's binary table."
            exit(1)
        try:
            dq_table = cos_tabledata.field("DQ")
        except KeyError:
            print "*** MAKE_HST_SPEC_PREVIEWS ERROR: DQ column not found in first extension's binary table."
            exit(1)
        """Create COSSegment objects to populate the COSSpectrum object with"""
        if band == 'FUV':
            try:
                fuva_index = numpy.where(segment_arr == "FUVA")[0][0]
            except IndexError:
                """Then there is no FUVA segment, normally this happends if there is only one segment present (in which case it's the other one of the two)."""
                fuva_index = None
            if fuva_index is not None:
                fuva_cossegment = COSSegment(nelem=nelems_arr[fuva_index], wavelengths=wavelength_table[fuva_index,:], fluxes=flux_table[fuva_index,:], fluxerrs=fluxerr_table[fuva_index,:], dqs=dq_table[fuva_index,:])
            try:
                fuvb_index = numpy.where(segment_arr == "FUVB")[0][0]
            except IndexError:
                """Then there is no FUVB segment, normally this happends if there is only one segment present (in which case it's the other one of the two)."""
                fuvb_index = None
            if fuvb_index is not None:
                fuvb_cossegment = COSSegment(nelem=nelems_arr[fuvb_index], wavelengths=wavelength_table[fuvb_index,:], fluxes=flux_table[fuvb_index,:], fluxerrs=fluxerr_table[fuvb_index,:], dqs=dq_table[fuvb_index,:])
        elif band == 'NUV':
            nuva_index = numpy.where(segment_arr == "NUVA")[0][0]
            nuva_cossegment = COSSegment(nelem=nelems_arr[nuva_index], wavelengths=wavelength_table[nuva_index,:], fluxes=flux_table[nuva_index,:], fluxerrs=fluxerr_table[nuva_index,:], dqs=dq_table[nuva_index,:])
            nuvb_index = numpy.where(segment_arr == "NUVB")[0][0]
            nuvb_cossegment = COSSegment(nelem=nelems_arr[nuvb_index], wavelengths=wavelength_table[nuvb_index,:], fluxes=flux_table[nuvb_index,:], fluxerrs=fluxerr_table[nuvb_index,:], dqs=dq_table[nuvb_index,:])
            nuvc_index = numpy.where(segment_arr == "NUVC")[0][0]
            nuvc_cossegment = COSSegment(nelem=nelems_arr[nuvc_index], wavelengths=wavelength_table[nuvc_index,:], fluxes=flux_table[nuvc_index,:], fluxerrs=fluxerr_table[nuvc_index,:], dqs=dq_table[nuvc_index,:])
        """Create COSSpectrum object."""
        if band == 'FUV':
            if fuva_index is not None and fuvb_index is not None:
                return_spec = COSSpectrum(band=band, cos_segments={'FUVA':fuva_cossegment,'FUVB':fuvb_cossegment}, orig_file=input_file)
            elif fuva_index is not None:
                return_spec = COSSpectrum(band=band, cos_segments={'FUVA':fuva_cossegment}, orig_file=input_file)
            else:
                return_spec = COSSpectrum(band=band, cos_segments={'FUVB':fuvb_cossegment}, orig_file=input_file)
        elif band == 'NUV':
            return_spec = COSSpectrum(band=band, cos_segments={'NUVA':nuva_cossegment,'NUVB':nuvb_cossegment,'NUVC':nuvc_cossegment}, orig_file=input_file)
        return return_spec
