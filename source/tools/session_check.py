# %%
"""Cheap integrity check of the local ``source/output`` data of this checkout.

Run by ``source/tools/session_check.sh`` (declared in ``AGENTS.md`` as the
repository's session check) with the repository root as the single argument.

It reports what CI cannot see: the gitignored ``output/`` archives and the
reference symlinks of *this machine*.  Every failure it lists is one a figure
driver would hit later, and that a workspace rename produces in bulk
(``AGENTS.md``, "stale data symlinks after a workspace rename").

This file is a **canonical shared piece**: it stays generic over any archive
that pairs ``t_<case>`` with a ``<case>`` array (``moments.*.npz``,
``corrmat_*.npz``, ``corr_*.npz`` alike), so ``dmrg4oqs`` keeps a synced copy
(``tools/sync_from_mkct.sh``) and declares it as its own session check.

Array *headers* only -- the archives are ~10 MB each and the plugin allows ten
seconds.  Exit 0 is clean and silent, 1 lists the problems one per line.
"""
import glob
import os
import sys
import zipfile

try:
    import yaml
    from numpy.lib import format as nf
except ImportError as error:      # the wrapper then tries the next interpreter
    sys.stderr.write("session_check: %s\n" % error)
    sys.exit(2)


def shape_of(archive, name):
    """The shape in an ``.npy`` member, read from its header alone."""
    with archive.open(name) as handle:
        version = nf.read_magic(handle)
        if version == (1, 0):
            return nf.read_array_header_1_0(handle)[0]
        if version == (2, 0):
            return nf.read_array_header_2_0(handle)[0]
        raise ValueError("npy format %r" % (version,))


def check(root):
    bad = []

    def rel(path):
        return os.path.relpath(path, root)

    out = os.path.join(root, "source", "output")
    data = os.path.join(out, "data")

    # 1. Dangling symlinks under data/ (a renamed checkout leaves them).  Only
    #    data/ is walked: output/cache is a multi-GB joblib tree and pycache is
    #    byte-compiled code, neither of which this check is about.
    if os.path.isdir(data):
        for dirpath, dirnames, filenames in os.walk(data):
            for name in dirnames + filenames:
                path = os.path.join(dirpath, name)
                if os.path.islink(path) and not os.path.exists(path):
                    bad.append("dangling symlink %s -> %s" % (rel(path), os.readlink(path)))

    # 2. Every config parses (truncated or mis-indented YAML silently defaults).
    for path in sorted(glob.glob(os.path.join(root, "source", "config", "*.yaml"))):
        try:
            with open(path) as handle:
                yaml.safe_load(handle)
        except Exception as error:
            bad.append("config %s: %s" % (rel(path), error))

    # 3. Every archive opens, every member header reads, and each `t_<case>` axis
    #    agrees in length with its `<case>` value array.
    for path in sorted(glob.glob(os.path.join(out, "data", "*", "*.npz"))):
        try:
            archive = zipfile.ZipFile(path)
        except Exception as error:
            bad.append("archive %s: not a zip (%s)" % (rel(path), error))
            continue
        names = [name for name in archive.namelist() if name.endswith(".npy")]
        keys = set(name[:-4] for name in names)
        shapes = {}
        for name in names:
            try:
                shapes[name[:-4]] = shape_of(archive, name)
            except Exception as error:
                bad.append("archive %s: unreadable %s (%s)" % (rel(path), name, error))
        for key in sorted(keys):
            if not key.startswith("t_"):
                continue
            case = key[2:]
            if case not in keys:
                bad.append("archive %s: %s without %s" % (rel(path), key, case))
            elif shapes.get(key) and shapes.get(case) and shapes[key][0] != shapes[case][0]:
                bad.append(
                    "archive %s: %s has %d points, %s has %d"
                    % (rel(path), key, shapes[key][0], case, shapes[case][0])
                )
    return bad


def main():
    if len(sys.argv) != 2:
        print("usage: session_check.py ROOT")
        return 2
    bad = check(sys.argv[1])
    for line in bad:
        print(line)
    return 0 if not bad else 1


# %%
if __name__ == "__main__":
    sys.exit(main())
