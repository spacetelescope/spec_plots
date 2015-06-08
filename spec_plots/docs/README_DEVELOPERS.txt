Best Practices for Development Of A New Branch

1.)  Crate a new branch with the name v<x.yy>, for example, if beginning work on branch 1.45, name the branch 'v1.45'.

2.)  You will need to update the README.md file inside the "spec_plots" folder to point to the new branch if you want the READTHEDOCS icons to be up-to-date.  You will also need to make sure the new branch is set to "active" at readthedocs.org after the first push to the new branch, so that it will auto-build after each git commit moving forwards.

3.)  Update the version number if python source files using the "update_version.sh" script inside the "spec_plots" folder (one level down from the root of the repo), so that the version numbers in source files are correctly updated.  If you move files around or add new python source files, be sure to update the shell script.

4.)  When doing git commits, refer to github Issue numbers when they are fixed in the commit notes so that github will assign this commit to that Issue, e.g., "- Fixes #44".

5.)  Re-build the EGG and WHELL distributions using the "build_package.sh" script at the top level.

6.)  Keep the various documentation files (mostly located inside the "docs" folder) as up-to-date as possible before final commit to the branch.

7.)  After doing a final commit of the branch, make sure all github Issues are closed for that branch, or re-assigned if punting to a later build for that Issue.

8.)  After regression testing is finished and all Issues assigned to this branch or closed for moved to a future build, initiate a pull request and merge into master.
