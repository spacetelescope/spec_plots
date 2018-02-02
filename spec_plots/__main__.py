"""
.. module:: __main__

   :synopsis: This script is used when make_hst_spec_previews is called directly
   as an executable from the command line.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

#--------------------
# Built-In Imports
#--------------------
from __future__ import absolute_import
import sys
#--------------------
# Package Imports
#--------------------
import spec_plots.make_hst_spec_previews as mhsp
from spec_plots import __version__

#--------------------

def main():
    """ This is the main function used when invoking make_hst_spec_previews
    directly as an executable. """

    # Collect command-line arguments and package into an ArgParse object.
    parser = mhsp.setup_args()
    args = parser.parse_args(sys.argv[1:])

    # Check arguments and options.
    mhsp.check_input_options(args)

    # Call main function.
    mhsp.make_hst_spec_previews(args.input_file,
                                flux_scale_factor=args.flux_scale_factor,
                                fluxerr_scale_factor=
                                args.fluxerr_scale_factor,
                                n_consecutive=args.n_consecutive,
                                output_path=args.output_path,
                                output_type=args.output_type,
                                dpi_val=args.dpi_val,
                                debug=args.debug,
                                full_ylabels=args.full_ylabels,
                                optimize=not args.nooptimize,
                                verbose=args.verbose)

if __name__ == "__main__":
    main()
