import sys
import spec_plots.make_hst_spec_previews as mhsp
#from spec_plots.make_hst_spec_previews import setup_args

def main():
    """ Collect command-line arguments and package into an ArgParse object. """
    parser = mhsp.setup_args()
    args = parser.parse_args(sys.argv[1:])

    """ Check arguments and options. """
    mhsp.check_input_options(args)

    """ Call main function. """
    mhsp.make_hst_spec_previews(args.input_file, \
                               flux_scale_factor = args.flux_scale_factor, \
                               fluxerr_scale_factor = args.fluxerr_scale_factor, \
                               n_consecutive = args.n_consecutive, \
                               output_path = args.output_path, \
                               output_type = args.output_type, \
                               dpi_val = args.dpi_val, \
                               debug = args.debug, \
                               full_ylabels = args.full_ylabels, \
                               verbose = args.verbose)

if __name__ == "__main__":
    main()
