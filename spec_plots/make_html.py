"""
.. module:: make_html
   :synopsis: Creates a webpage of thumbnail and full-size previews for review
       purposes.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

#--------------------
# Built-In Imports
#--------------------
from __future__ import absolute_import
from __future__ import print_function
import argparse
from glob import glob
import os
import sys
from builtins import str
#--------------------
# External Imports
#--------------------
import numpy
#--------------------
# Package Imports
#--------------------
from spec_plots import __version__

#--------------------

def find_orig_preview(ufr, orig_files):
    """
    Given aunique file root (`ufr`, which includes the IPPPSSOOT ID), this
    function matches that to the list of original preview files from CADC to
    determine if an original preview exists for this ufr.

    :param ufr: Unique file root to test (this is equivalent to the IPPPSSOOT ID
        from MAST).

    :type ufr: str

    :param orig_files: Array of original preview file names in which to find a
        match.

    :type orig_files: numpy.ndarray

    :returns: int -- The index in `orig_files` that matches the specified `ufr`,
        if there is no match then returns None.
    """

    # Find all indices where the ufr is a substring of the original CADC
    # preview plot file name.  Note that the ufr is input as <ufr>_<size>, so
    # split on underscore and just look for the first part.
    where_match = [i for i, x in enumerate(orig_files) if
                   ufr.split('_')[0] in x]

    # Return the index if found.
    if len(where_match) == 1:
        return where_match[0]
    elif len(where_match) > 1:
        print("*** WARNING in MAKE_HTML: Found more than one CADC match, using"
              " the first one for file " + ufr)
        return where_match[0]
    else:
        return None

#--------------------

def make_html(idir, odir="html/plot_previews/", ofile="plot_previews",
              orig_dir=None, plot_display_width=512):
    """
    Creates diagnostic HTML page showing the generated preview plots (thumb and
    large) as well as the original MAST version of these plots (the one's from
    CADC, provided they can be located in the expected location).

    :param idir: Input directory containing the set of preview images generated
        by spec_plots.

    :type idir: str

    :param odir: Output path for the HTML files.  Default = html/plot_previews/`

    :type odir: str

    :param ofile: Base file name of the output HTML files.  Default =
        plot_previews`

    :type ofile: str

    :param orig_dir: [Optional] Directory containing the original preview plots
        from CADC.  If provided, these will be included in the HTML pages for
        comparison purposes.

    :type orig_dir: str

    :param plot_display_width: [Optional] The width of the preview plots in the
        HTML tables (specified through the <img> tag), in pixels.  Default =
        512.

    :type plot_display_width: int

    :raises: OSError, IOError, ValueError
    """

    # Ensure output directory is treated as an absolute path, make sure it
    # exists, if not, create the output directory.
    odir = os.path.abspath(odir) + os.path.sep
    if odir != '' and not os.path.isdir(odir):
        try:
            os.mkdir(odir)
        except OSError as this_error:
            if this_error.errno == 13:
                sys.stderr.write("*** MAKE_HTML ERROR: Output directory could"
                                 " not be created, "+repr(this_error.strerror)+
                                 "\n")
                exit(1)
            else:
                raise

    # Get all of the thumb and full-size .png files in this directory.
    if idir is not None:
        if os.path.isdir(idir):
            all_thumb_png_files = numpy.asarray(glob(idir+"/*0128.png"))
            n_thumb_png_files = len(all_thumb_png_files)

            all_large_png_files = numpy.asarray(glob(idir+"/*1024.png"))
            n_large_png_files = len(all_large_png_files)

            # <DEVEL> The code blocks that create the file roots should be
            # generators and not arrays, since the array size might be quite
            # large. </DEVEL>
            # Get the file root of the thumb-size spec_plots preview images,
            # which are the first two parts of the file name separated by
            # underscores.  E.g., la7803fkq_x1d_0128.png -> la7803fkq_x1d.
            if n_thumb_png_files > 0:
                all_thumb_png_files_froots = numpy.asarray(
                    ['_'.join(os.path.basename(x).split('_')[0:2])
                     for x in all_thumb_png_files])
                thumb_sort_indexes = numpy.argsort(all_thumb_png_files_froots)
                all_thumb_png_files_froots = (
                    all_thumb_png_files_froots[thumb_sort_indexes])
                all_thumb_png_files = all_thumb_png_files[thumb_sort_indexes]


            # Get the file root for the large-size spec_plots preview files.
            if n_large_png_files > 0:
                all_large_png_files_froots = numpy.asarray(
                    ['_'.join(os.path.basename(x).split('_')[0:2])
                     for x in all_large_png_files])
                large_sort_indexes = numpy.argsort(all_large_png_files_froots)
                all_large_png_files_froots = (
                    all_large_png_files_froots[large_sort_indexes])
                all_large_png_files = all_large_png_files[large_sort_indexes]

        else:
            raise IOError("Could not find directory " + idir)

        # If requested, get a list of all the original preview plots.
        if orig_dir is not None:
            if os.path.isdir(orig_dir):
                all_orig_files_withpath = numpy.asarray(glob(orig_dir+"/*.jpg"))
                all_orig_files = numpy.asarray(
                    [os.path.basename(x).lower() for x in
                     all_orig_files_withpath])
            else:
                raise IOError("Could not find directory containing original"
                              " preview plots.  Looking for " + orig_dir)

        # Get a list of unique <IPPPSSOOT_filetype> base names from both the
        # all_thumb and all_large arrays of spec_plots file names.  This is
        # because, in principle, there might only be a preview thumbnail or a
        # full-size preview for a given IPPPSSOOT_filetype.  Note that this
        # array should be sorted since the "unique" function returns a sorted
        #  array.
        # <DEVEL> Concatenating the arrays here can be expensive if trying to
        # make webpage for a very large number of preview plots.
        # Consider an alternative way of getting the unique set between the two
        # without creating a new array. </DEVEL>
        if n_thumb_png_files > 0 and n_large_png_files > 0:
            uniq_fileroots = numpy.sort(numpy.unique(
                numpy.concatenate([all_thumb_png_files_froots,
                                   all_large_png_files_froots])))
        elif n_thumb_png_files > 0:
            uniq_fileroots = numpy.sort(numpy.unique(
                all_thumb_png_files_froots))
        elif n_large_png_files > 0:
            uniq_fileroots = numpy.sort(numpy.unique(
                all_large_png_files_froots))
        else:
            raise IOError("No thumb-size or large-size files found."
                          "  Looking in " + orig_dir)

        # How many thumbnails per row?
        n_thumbs_per_row = 8

        # Make the HTML file for the thumbnail-sized preview plots.
        with open(odir+ofile+"_thumbs.html", 'w') as ofi:
            ofi.write('<html><head></head><body>\n')
            ofi.write('<table style="border:1px solid black;border-collapse:'
                      'collapse;width:100%">\n')

            for i, ufr in enumerate(uniq_fileroots):

                # Check if the thumb version of the preview exists.
                if n_thumb_png_files > 0 and ufr in all_thumb_png_files_froots:
                    has_thumb = True
                    where_this_thumb = numpy.where(
                        all_thumb_png_files_froots == ufr)[0]
                    if len(where_this_thumb) > 1:
                        print("Warning: This unique file root appeared more"
                              "than once in the list of all thumbnails: " +
                              ufr)
                    where_this_thumb = where_this_thumb[0]
                else:
                    has_thumb = False

                if has_thumb:
                    if i % n_thumbs_per_row == 0:
                        ofi.write('  <tr>\n')

                    # Write the cell containing the thumbnail preview, (or
                    # just fill it grey if missing).
                    ofi.write('    <td style="border:1px solid black;width:'
                              '135px;vertical-align:top"><div style="width:'
                              '128px;text-align:center"><span style='
                              '"font-weight:bold">'+'_<br>'.join(ufr.split('_'))
                              +'</span></div>')
                    if has_thumb:
                        ofi.write('    <img src="'+os.path.relpath(
                            all_thumb_png_files[where_this_thumb], odir)+
                                  '" width="128px">')
                    else:
                        ofi.write('    <div style="background-color:#86867D;'
                                  'width:128px;height:128px"></div>')
                    ofi.write('    </td>\n')

                    if i % n_thumbs_per_row == n_thumbs_per_row-1:
                        ofi.write('  </tr>\n')

                else:
                    print("Warning: Could not find full-size PNG for"
                          " IPPPSSOOT_filetype = " + ufr)
            ofi.write('</table>/n')
            ofi.write('</body></html>/n')

        # Make the HTML file for the large-sized preview plots.
        with open(odir+ofile+"_large.html", 'w') as ofi:
            ofi.write('<html><head></head><body>\n')
            ofi.write('<table style="border:1px solid black;border-collapse:'
                      'collapse;width:'+str(int(round(2.*plot_display_width)))+
                      'px">\n')
            ofi.write('<tr><th>New Version</th><th>Old Version</th></tr>\n')


            for i, ufr in enumerate(uniq_fileroots):

                # Check if the large version of the preview exists.
                if n_large_png_files > 0 and ufr in all_large_png_files_froots:
                    has_large = True
                    where_this_large = numpy.where(
                        all_large_png_files_froots == ufr)[0]
                    if len(where_this_large) > 1:
                        print("Warning: This unique file root appeared more"
                              " than once in the list of all large previews: "+
                              ufr)
                    where_this_large = where_this_large[0]
                else:
                    has_large = False

                # If it has at least one preview image, then print this table
                # row.
                if has_large:
                    ofi.write('  <tr>\n')

                    # Write the cell containing the large preview, (or just
                    # fill it grey if missing).
                    ofi.write('    <td style="border:1px solid black;width:'
                              '100%">')
                    if has_large:
                        ofi.write('    <img src="'+os.path.relpath(
                            all_large_png_files[where_this_large], odir)+
                                  '" width="'+str(plot_display_width)+'px">')
                    else:
                        ofi.write('    <div style="background-color:#86867D;'
                                  'width:'+str(plot_display_width)+'px;height:'+
                                  str(plot_display_width)+'px"></div>')
                    ofi.write('    </td>\n')

                    # If requested, write the cell containing the original
                    # preview plot.  If requested but missing, just fill it with
                    # grey.
                    if orig_dir is not None:
                        # Check if this file root has a match in the original
                        # directory.
                        orig_preview_index = find_orig_preview(ufr,
                                                               all_orig_files)
                        ofi.write('    <td style="border:1px solid black;width:'
                                  '100%">')
                        if orig_preview_index is not None:
                            ofi.write('    <img src="'+os.path.relpath(
                                all_orig_files_withpath[orig_preview_index],
                                odir)+'" width="'+str(plot_display_width)+
                                      'px">')
                        else:
                            ofi.write('    <div style="background-color:'
                                      '#86867D;width:'+str(plot_display_width)+
                                      'px;height:'+str(plot_display_width)+
                                      'px"></div>')
                        ofi.write('    </td>\n')

                    ofi.write('  </tr>\n')

                else:
                    print("Warning: Could not find full-size PNG for "
                          "IPPPSSOOT_filetype = " + ufr)

            ofi.write('</table>\n')
            ofi.write('</body></html>\n')

    else:
        raise ValueError("No preview plot directory specified.")

#--------------------

if __name__ == "__main__":

    # Create argument parser.
    PARSER = argparse.ArgumentParser(description="Create HTML page of preview"
                                     " plots, given a directory containing the"
                                     " plots.")

    PARSER.add_argument("input_dir", action="store", type=str, help="[Required]"
                        " Full path to directory containing preview plots.")

    PARSER.add_argument("-f", action="store", type=str, dest="output_file",
                        default="plot_previews", help="[Optional] Base file"
                        " name of the output HTML files.  If the files already"
                        " exist, they will be overwritten.  Defualt ="
                        " '%(default)s'.")

    PARSER.add_argument("-o", action="store", type=str, dest="output_dir",
                        default="html/plot_previews/", help="[Optional] Specify"
                        " the full path to the output HTML files.  Defualt ="
                        " '%(default)s'.")

    PARSER.add_argument("-p", action="store", type=str, dest="orig_dir",
                        help="[Optional] Full path to a directory containing"
                        " the original versions of the plot previews from CADC."
                        "  If provided, these will be included in the HTML"
                        " table for comparision purposes.")

    PARSER.add_argument("-w", action="store", type=float, dest=
                        "plot_display_width", default=512., help="[Optional]"
                        " Specify the display width of the large-sized preview"
                        " plots, in pixels, for the HTML table.  Default = "
                        "'%(default)s.'")

    # Parse arguments.
    ARGS = PARSER.parse_args()

    # Make sure the requested display width for the large-sized plots is at
    # least greater than a minimum value.
    MIN_DISPLAY_WIDTH = 128.
    if ARGS.plot_display_width < MIN_DISPLAY_WIDTH:
        print("Warning: Display width for full-size preview plots is very,"
              " very small.  Will use a display width of "+
              str(MIN_DISPLAY_WIDTH)+" px instead of "+
              str(ARGS.plot_display_width)+" px.")
        ARGS.plot_display_width = MIN_DISPLAY_WIDTH

    # Call main method.
    make_html(ARGS.input_dir, ARGS.output_dir, ARGS.output_file, ARGS.orig_dir,
              ARGS.plot_display_width)
#--------------------
