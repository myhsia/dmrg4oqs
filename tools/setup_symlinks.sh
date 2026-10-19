#!/usr/bin/env bash
# Recreate the symlinks to the shared files in the sibling `mkct` checkout.
# These links are gitignored: the public repository ships only the DMRG code,
# and CI copies the shared files from the private repository instead.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"   # ../dmrg
mkct="$(cd "$here/../mkct" && pwd)"

mkdir -p "$here/source/library"
ln -sfn ../../mkct/source/config                 "$here/source/config"
ln -sfn ../../../mkct/source/library/qutils.py   "$here/source/library/qutils.py"
ln -sfn ../../../mkct/source/library/fsutils.py  "$here/source/library/fsutils.py"
ln -sfn ../../../mkct/source/library/pltools.py  "$here/source/library/pltools.py"
ln -sfn ../../mkct/.opencode                     "$here/.opencode"
echo "symlinks -> $mkct"
