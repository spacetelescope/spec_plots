__version__ = '1.2'

"""
.. module:: specutils_stis

   :synopsis: Contains functions for reading and plotting HST STIS spectra.

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

#--------------------

class STIS1DSpectrum:
    """
    Defines a STIS 1D spectrum (either "x1d" extracted or "sx1" summed extracted), including wavelegnth, flux, and flux errors.  A STIS 1D spectrum object consists of N associations.  If the file is an association, then N > 1, otherwise N = 1.  Each of these N associations can contain M orders.  If the association is an Echelle spectrum, then 24 < M < 70, depending on instrument configuration, otherwise M = 1.  Each of these M orders contain typical spectral data (wavelengths, fluxes, etc.), stored as STISOrderSpectrum objects.  The final data structure is then <STIS1DSpectrum>.associations[n].order[m].wavelengths (or .fluxes, .fluxerrs, etc.).
    """
    def __init__(self, association_spectra=None, orig_file=None):
        """
        Create a STIS1DSpectrum object out of a list of STISExposureSpectrum objects, which themselves are lists of STISOrderSpectrum objects.

        :param association_spectra: A list whose length is equal to the number of associations (length = "N" associations).  Each element in this list is a STISExposureSpectrum object, which itself is a list (length = "M" orders) of STISOrderSpectrum objects.

        :type association_spectra: list

        :param orig_file: Original FITS file read to create the spectrum (includes full path).

        :type orig_file: str

        :raises: ValueError
        """
        if association_spectra is not None:
            self.orig_file = orig_file
            self.associations = association_spectra
        else:
            raise ValueError("Must provide a list of STISExposureSpectrum objects.")

#--------------------

class STISExposureSpectrum:
    """
    Defines a STIS exposure spectrum, which consists of "M" STISOrderSpectrum objects.
    """
    def __init__(self, order_spectra=None):
        """
        Create a STISExposureSpectrum object out of a list of STISOrderSpectrum objects.

        :param order_spectra: The STISOrderSpectrum objects to build the STISExposureSpectrum object out of.

        :type order_spectra: list

        :raises: ValueError
        """
        if order_spectra is None:
            raise ValueError("Must provide a list of STISOrderSpectrum objects.")
        elif len(order_spectra) > 0:
            self.orders = order_spectra
        else:
            raise ValueError("Must provide a list of at least one STISOrderSpectrum objects, input list is empty.")

#--------------------

class STISOrderSpectrum:
    """
    Defines a STIS order spectrum, including wavelength, flux, and flux errors, which are stored as numpy arrays.  A scalar int property provides the number of elements in this segment.
    """
    def __init__(self, nelem=None, wavelengths=None, fluxes=None, fluxerrs=None, dqs=None):
        """
        Create a STISOrderSpectrum object, default to empty values.  Allows user to preallocate space if they desire by setting "nelem" but not providing lists/arrays on input right away.

        :param nelem: Number of elements for this segment's spectrum.

        :type nelem: int

        :param wavelengths: List of wavelengths in this segment's spectrum.

        :type wavelengths: list

        :param fluxes: List of fluxes in this segment's spectrum.

        :type fluxes: list

        :param fluxerrs: List of flux uncertainties in this segment's spectrum.

        :type fluxerrs: list

        :param dqs: List of data quality flags.

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

#--------------------

def generate_stis_avoid_regions():
    """
    Creates a list of AvoidRegion objects for use in the plotting routine, specifically designed for STIS spectra.
    """
    lya1215_ar = specutils.AvoidRegion(1214.,1217., "Lyman alpha emission line.")
    return [lya1215_ar]

#--------------------

def plotspec(stis_spectrum, output_type, output_file, n_consecutive, flux_scale_factor, fluxerr_scale_factor, output_size=None, debug=False, full_ylabels=False):
    """
    Accepts a STIS1DSpectrum object from the READSPEC function and produces preview plots.

    :param stis_spectrum: STIS spectrum as returned by READSPEC.

    :type stis_spectrum: STIS1DSpectrum

    :param output_type: What kind of output to make?

    :type output_type: str

    :param output_file: Name of output file (including full path).

    :type output_file: str

    :param n_consecutive: How many consecutive points must pass the test for the index to count as the valid start/end of the spectrum?  Default = 20.

    :type n_consecutive: int

    :param flux_scale_factor: Max. allowed ratio between the flux and a median flux value, used in edge trimming.  Default = 10.

    :type flux_scale_factor: float

    :param fluxerr_scale_factor: Max. allowed ratio between the flux uncertainty and a median flux uncertainty value, used in edge trimming.  Default = 5.

    :type fluxerr_scale_factor: float

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

    """If the figure size is large, then plot up to three associations, otherwise force only one association on the plot."""
    n_associations = len(stis_spectrum.associations)
    if output_size > 128:
        is_bigplot = True
        if n_associations <= 3:
            n_subplots = n_associations
            subplot_indexes = range(n_associations)
        else:
            n_subplots = 3
            midindex = int(round(float(n_associations)/2.))
            subplot_indexes = [0,midindex,n_associations-1]
    else:
        is_bigplot = False
        n_subplots = 1
        subplot_indexes = [0]

    """Start plot figure."""
    this_figure, these_plotareas = pyplot.subplots(nrows=n_subplots, ncols=1, figsize=(output_size/dpi_val, output_size/dpi_val), dpi=dpi_val)

    if not isinstance(these_plotareas, numpy.ndarray):
        these_plotareas = numpy.asarray([these_plotareas])
    if is_bigplot:
        this_figure.subplots_adjust(hspace=0.3,top=0.915)
        this_figure.suptitle(os.path.basename(stis_spectrum.orig_file))
    else:
        this_figure.subplots_adjust(top=0.85,bottom=0.3,left=0.25,right=0.9)

    for c,i in enumerate(subplot_indexes):
        this_plotarea = these_plotareas[c]
        
        all_wls, all_fls, all_flerrs, all_dqs, title_addendum = specutils.stitch_components(stis_spectrum.associations[i], n_consecutive, flux_scale_factor, fluxerr_scale_factor)

        median_flux, median_fluxerr, fluxerr_95th = specutils.get_flux_stats(all_fls, all_flerrs)

        if is_bigplot:
            this_plotarea.set_title(title_addendum, loc="right", size="small", color="red")
        if n_associations > 1 and is_bigplot:
            this_plotarea.set_title("Association "+str(i+1)+"/"+str(n_associations), loc="center", size="small", color="black")

        """Determine optimal x-axis."""
        x_axis_range = specutils.set_plot_xrange("stis", all_wls, all_fls, all_flerrs, all_dqs, n_consecutive, flux_scale_factor, fluxerr_scale_factor, median_flux, median_fluxerr, fluxerr_95th)
        if all(numpy.isfinite(x_axis_range)):
            """Create COS avoid regions."""
            avoid_regions = generate_stis_avoid_regions()
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

#--------------------

def readspec(input_file):
    """
    Reads in a STIS spectrum FITS file (x1d, sx1) and returns the wavelengths, fluxes, and flux uncertainties for the two (FUV segments) or three (NUV stripes).

    :param input_file: Name of input FITS file.

    :type input_file: str

    :returns: STIS1DSpectrum -- The spectroscopic data (wavelength, flux, flux error, etc):
    """
    with fits.open(input_file) as hdulist:
        """ Read in the number of extensions from the primary header.  This will determine whether this is an association (N > 1) or not (N = 1)."""
        n_associations = hdulist[0].header["NEXTEND"]
        """Create an initially empty list that will contain each extension's (association's) spectrum object."""
        all_association_spectra = []
        for exten in xrange(n_associations):
            exten_data_table = hdulist[exten+1].data
            """How many orders (table rows) in this extension?"""
            n_orders = len(exten_data_table["sporder"])
            """Create a list of STISOrderSpectra for this extension."""
            all_order_spectra = [STISOrderSpectrum(nelem=exten_data_table["nelem"][order], wavelengths=exten_data_table["WAVELENGTH"][order], fluxes=exten_data_table["FLUX"][order], fluxerrs=exten_data_table["ERROR"][order], dqs=exten_data_table["DQ"][order]) for order in xrange(n_orders)]
            """Create a STISExposureSpectrum from the STISOrderSpectrum objects."""
            this_exposure_spectrum = STISExposureSpectrum(order_spectra=all_order_spectra)
            all_association_spectra.append(this_exposure_spectrum)
        return STIS1DSpectrum(association_spectra=all_association_spectra, orig_file=input_file)

#--------------------
