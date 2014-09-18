__version__ = '1.0'

"""
.. module:: make_html
   :synopsis: Creates a webpage of thumbnail and full-size previews for review purposes.
.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

from glob import glob
import os
import argparse
import numpy

def make_html(idir=None, ofile="plot_previews.html"):
    if idir is not None:
        """Get all of the thumb and full-size .png files at this location."""
        if os.path.isdir(idir):
            all_thumb_png_files = numpy.asarray(glob(idir+'/*0128.png'))
            all_large_png_files = numpy.asarray(glob(idir+'/*1024.png'))
            n_thumb_png_files = len(all_thumb_png_files); n_large_png_files = len(all_large_png_files)
            if n_thumb_png_files > 0:
                all_thumb_png_files_froots = ['_'.join(os.path.basename(x).split('_')[0:2]) for x in all_thumb_png_files]
                thumb_sort_indexes = numpy.argsort(all_thumb_png_files_froots)
                all_thumb_png_files = all_thumb_png_files[thumb_sort_indexes]
            if n_large_png_files > 0:
                all_large_png_files_froots = ['_'.join(os.path.basename(x).split('_')[0:2]) for x in all_large_png_files]
                large_sort_indexes = numpy.argsort(all_large_png_files_froots)
                all_large_png_files = all_large_png_files[large_sort_indexes]
        else:
            raise IOError("Could not find directory " + idir)
        """Get a list of unique <IPPPSSOOT_filetype> base names from both the all_thumb and all_large arrays.  This is because, in principle, there might only be a preview thumbnail or a full-size preview for a given IPPPSSOOT_filetype."""
        all_fileroots = numpy.concatenate([all_thumb_png_files, all_large_png_files])
        for i,froot in enumerate(all_fileroots):
            all_fileroots[i] = '_'.join(os.path.basename(froot).split('_')[0:2])
        """Note that this array should be sorted since the "unique" function returns a sorted array."""
        uniq_fileroots = numpy.unique(all_fileroots)
        """Open HTML for writing and begin printing HTML table, where each row is one of the unique fileroots."""
        cur_thumb_index = 0; cur_large_index = 0
        with open(ofile, 'w') as of:
            of.write('<html><head></head><body>\n')
            of.write('<table style="border:1px solid black;border-collapse:collapse;width:1160px">\n')
            for ufr in uniq_fileroots:
                if n_thumb_png_files > 0:
                    if '_'.join(os.path.basename(all_thumb_png_files[cur_thumb_index]).split('_')[0:2]) == ufr:
                        has_thumb = True
                    else:
                        has_thumb = False
                else:
                    has_thumb = False
                if n_large_png_files > 0:
                    if '_'.join(os.path.basename(all_large_png_files[cur_large_index]).split('_')[0:2]) == ufr:
                        has_large = True
                    else:
                        has_large = False
                else:
                    has_large = False

                if has_thumb or has_large:
                    of.write('  <tr>\n')
                    """Write the cell containing the thumbnail preview, (or just fill it grey if missing)."""
                    of.write('    <td style="border:1px solid black;width:135px;vertical-align:top"><div style="width:128px;text-align:center"><span style="font-weight:bold">'+ufr+'</span></div>')
                    if has_thumb:
                        of.write('<img src="'+all_thumb_png_files[cur_thumb_index]+'" width="128px">')
                    else:
                        of.write('<div style="background-color:#86867D;width:128px;height:128px"></div>')
                    of.write('</td>\n')

                    """Write the cell containing the large preview, (or just fill it grey if missing)."""
                    of.write('    <td style="border:1px solid black;width:1030px">')
                    if has_large:
                        of.write('<img src="'+all_large_png_files[cur_large_index]+'" width="1024px">')
                    else:
                        of.write('<div style="background-color:#86867D;width:1024px;height:1024px"></div>')
                    of.write('</td>\n')

                    of.write('  </tr>\n')
                else:
                    print "Warning: Could not find either thumbnail or full-size PNG for IPPPSSOOT_filetype = " + ufr
                    import ipdb; ipdb.set_trace()
                if has_thumb:
                    cur_thumb_index+=1
                if has_large:
                    cur_large_index+=1
            of.write('</table>\n')
            of.write('</body></html>\n')
    else:
        raise ValueError("No preview plot directory specified.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create HTML page of preview plots, given an output directory.")
    parser.add_argument("-d", action="store", type=str, dest="input_dir", default=None, help="[Required] Full path to directory containing preview plots.",metavar='location of plots')
    parser.add_argument("-o", action="store", type=str, dest="output_file", default="plot_previews.html", help='[Optional] Full path and file name of the output HTML file.  If the file already exists, it will be overwritten.  Defualt = "plot_previews.html".')
    args = parser.parse_args()
    make_html(args.input_dir, args.output_file)
