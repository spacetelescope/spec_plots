#!/bin/csh

## This script updates the version string inside python programs.  Note that it does NOT update the version inside the main README file or other documentation files, it only updates python source files.

## Define the lists of files that need to be updated with the new version.  This includes the setup.py file, and, separately, the python source files (since the version information is recorded in setup.py differently).
set setupfile = "../setup.py"
set files = ( "make_hst_spec_previews.py" "make_html.py" "utils/specutils/dq_has_flag.py" "utils/specutils/rms.py" "utils/specutils/stitch_components.py" "utils/specutils/avoidregion.py" "utils/specutils/edge_trim.py" "utils/specutils/set_plot_xrange.py" "utils/specutils/calc_plot_metrics.py" "utils/specutils/calc_covering_fraction.py" "utils/specutils/get_flux_stats.py" "utils/specutils/set_plot_yrange.py" "utils/specutils/debug_oplot.py" "utils/specutils/is_bad_dq.py" "utils/specutils/specutilserror.py" "utils/specutils_cos/extract_subspec.py" "utils/specutils_cos/readspec.py" "utils/specutils_cos/check_segments.py" "utils/specutils_cos/get_segment_names.py" "utils/specutils_cos/cosspectrum.py" "utils/specutils_cos/plotspec.py" "utils/specutils_stis/plotspec.py" "utils/specutils_stis/stis1dspectrum.py" "utils/specutils_stis/get_association_indices.py" "utils/specutils_stis/readspec.py" "__main__.py")


## Print out the help statement if incorrect number of arguments are present.
if ($#argv != 2) then
    echo "Usage: update_version old_ver_string new_ver_string"
    echo "Example: update_version 1.25 1.30"
    exit 1
endif

set old_ver = $1
set new_ver = $2

## Replace the version inside setup.py file.
echo Running command: sed -i -e s/version=\"$old_ver\"/version=\"$new_ver\"/ $setupfile
sed -i -e s/version=\"$old_ver\"/version=\"$new_ver\"/ $setupfile

## Replace the version inside each of the python source files.
foreach file ($files)
     if (-f $file ) then
	 echo Running command: sed -i -e s/__version__\ =\ \'$old_ver\'/__version__\ =\ \'$new_ver\'/ $file
	 sed -i -e s/__version__\ =\ \'$old_ver\'/__version__\ =\ \'$new_ver\'/ $file
     else
	 echo *** ERROR: Missing file ${file}.
    endif
end
