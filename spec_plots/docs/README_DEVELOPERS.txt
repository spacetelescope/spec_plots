Best Practices for Development Of A New Branch

1.)  Crate a new branch with the name v<x.yy>, for example, if beginning work on branch 1.45.0, name the branch 'v1.45.0'.

2.)  You will need to update the README.md file inside the "spec_plots" folder to point to the new branch if you want the READTHEDOCS icons to be up-to-date.  You will also need to make sure the new branch is set to "active" at readthedocs.org after the first push to the new branch, so that it will auto-build after each git commit moving forwards.

3.)  Update the version number in the __init__.py under spect_plots/.  This version is imported by all other modules, so it defines the version of the software throughout the repo.

4.)  When doing git commits, refer to github Issue numbers when they are fixed in the commit notes so that github will assign this commit to that Issue, e.g., "- Fixes #44".

5.)  Upate the History.rst file in the main directory with a summary of the changes, used by PyPI to track changes between versions.

6.)  Re-build the EGG and WHEEL distributions using the "build_package.sh" script at the top level.

7.)  Keep the various documentation files (mostly located inside the "docs" folder) as up-to-date as possible before final commit to the branch.

8.)  After doing a final commit of the branch, make sure all Github Issues are closed for that branch, or re-assigned if punting to a later build for that Issue.

9.)  After regression testing is finished and all Issues assigned to this branch are closed or moved to a future build, initiate a pull request and merge into master.

10.)  To upload to PyPI with twine: twine upload dist/*, with the -u and -p options.

11.)  After uploading to PyPI, you can build the conda package.  Inside the "conda" folder in the top-level directory, run "conda skeleton pypi spec-plots", then "conda build spec-plots".  You can also make a Python 3.6 version by doing "conda build --python 3.6 spec-plots".  NOTE: On SWF's environment, it's necessary to keep anaconda Python at "2" but the overall Python environment at 3.6 when building that version.

12.)  Test the build using "conda install --user-local spec-plots".  If everything looks good, upload it to anaconda.org using "anaconda upload /Users/fleming/anaconda2/conda-bld/osx-64/spec-plots-1.34.2-py27_0.tar.bz2", where the last part is the name of the built .tar.bz2 file.
