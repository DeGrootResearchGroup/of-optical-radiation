# Common Allclean body shared by every test case.
#
# Sourced from each case's Allclean *after* the caller cd's into the
# case directory. Standard cleanup: source OpenFOAM's CleanFunctions,
# run cleanCase, and remove the 0/ time directory. Cases needing
# additional cleanup (postProcessing, multi-region polyMesh, etc.)
# add their own `rm -rf` lines after the source.
#
# Why a sourced fragment rather than a sub-script: the caller has
# already chosen its working directory and shell context, so sharing a
# function/macro is unnecessary -- a sourced fragment runs in the
# caller's shell and is the smallest possible abstraction.

. "$WM_PROJECT_DIR/bin/tools/CleanFunctions"

cleanCase
rm -rf 0
