#!/bin/csh

## This script updates the version string inside python programs.  Note that it does NOT update the version inside the main README file or other documentation files, it only updates python source files.

## Define the lists of files that need to be updated with the new version.  This includes the setup.py file, and, separately, the python source files (since the version information is recorded in setup.py differently).
set setupfile = "setup.py"
set files = ("make_hst_spec_previews.py" "make_html.py" "specutils.py" "specutils_cos.py" "specutils_stis.py")

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
