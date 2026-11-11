#!/usr/bin/env bash
# The repository's session check, declared in `AGENTS.md` as
# `` Session check: `source/tools/session_check.sh` `` and run by the `upkeep`
# plugin at session start with this repository as the working directory.
#
# `session_check.py` does the work (array headers only: the archives are ~10 MB
# each and the plugin's budget is ten seconds).  This wrapper only finds an
# interpreter that can import numpy and yaml, without pre-probing each one --
# a probe costs as much as the check itself.  Each candidate is tried in turn and
# one that reports "cannot run here" (exit 2) is simply skipped.
#
# Exit codes, per the plugin contract: 0 clean (silent), 1 something to report
# (one line per problem; the plugin keeps the first three), 2 cannot run here.
set -uo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"     # .../<repo>
check="$root/source/tools/session_check.py"

candidates=("${MKCT_PYTHON:-}")
[ -n "${CONDA_PREFIX:-}" ] && candidates+=("$CONDA_PREFIX/bin/python")
condabase="$(conda info --base 2>/dev/null || true)"
[ -n "$condabase" ] && candidates+=("$condabase/envs/mkct/bin/python")
[ -n "$(command -v python3 2>/dev/null)" ] && candidates+=("$(command -v python3)")
[ -n "$(command -v python 2>/dev/null)" ] && candidates+=("$(command -v python)")

for cand in "${candidates[@]}"; do
  [ -n "$cand" ] && [ -x "$cand" ] || continue
  "$cand" "$check" "$root" 2>/dev/null
  status=$?
  [ "$status" -eq 2 ] || exit "$status"
done

if command -v conda >/dev/null 2>&1; then
  out="$(conda run -n mkct python "$check" "$root" 2>/dev/null)"
  status=$?
  [ "$status" -eq 2 ] || { [ -n "$out" ] && printf '%s\n' "$out"; exit "$status"; }
fi

echo "no interpreter with numpy and yaml (set MKCT_PYTHON or activate mkct)"
exit 2
