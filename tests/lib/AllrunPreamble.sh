# Common Allrun preamble. Sourced from each case's Allrun right after
# the shebang line.
#
# Does the three things every case Allrun was doing by hand:
#   1. cd into the case directory (so the script is robust to being
#      invoked via an absolute or relative path).
#   2. Source RunFunctions so the rest of the script can use
#      runApplication / runParallel.
#   3. Pre-clean the case so runApplication does not skip on a stale
#      log.<app> from a previous (possibly aborted) run.
#
# What this preamble deliberately does NOT do: restore 0/ from 0.orig.
# That step varies across cases -- some need `rm -rf 0` first, some
# don't have a 0.orig at all, and the three multi-region cases
# restore per-region 0/ after splitMeshRegions -- so each Allrun
# carries its own restore lines.

cd "${0%/*}" || exit 1

. "$WM_PROJECT_DIR/bin/tools/RunFunctions"

./Allclean > /dev/null 2>&1 || true
