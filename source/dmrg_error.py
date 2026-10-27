# %%
"""Accuracy of the MPS/tDMRG references against HEOM.

The metric is the shared correlation-matrix error of ``library/mre.py``: over
the calculation window ``[0, tf]`` of each config, the global relative Frobenius
error ``E = ||C - C'||_F / ||C'||_F`` of the basis correlation matrix (all 16
complex elements) against the HEOM reference matrix, reported in percent.
``tf`` is the config window (forced, for comparability with the other
references); a run whose series ends earlier is compared on the overlap only
and a warning is emitted.

This is the runner; the algorithm is imported from the shared ``mre`` module.
It writes ``output/dmrg_error.json`` (one value per config/case/label, each row
carrying the window ``tf``) and ``output/dmrg_error.md``, and prints the table.

Usage::

    python dmrg_error.py                                  # RUNLIST x cases
    python dmrg_error.py --cases qhm qah                  # quantum only
    python dmrg_error.py --labels "" star64 star128        # tagged runs
"""
# %%
# 1. Requirements
import argparse, json, os, sys, warnings
sys.path.append(os.path.join(os.path.dirname(__file__), "./library"))
sys.pycache_prefix = os.path.abspath("./output/pycache")
import mre as core
from fsutils import load_config, load_group, ensure_dir
from pltools import parse_xlims


# 2. Run selection
RUNLIST = [
    "A-base_line",
    "B-slow_system",
    "C-no_bias",
    "D-strong_tunneling",
    "E-fast_bath",
    "F-slow_bath",
    "G-lowT",
    "H-highT",
    "I-weak_coupling",
    "J-strong_coupling",
    "K-linear_quad",
]
CASES = ["qhm", "qah", "chm", "cah"]
METHODS = ["dmrg"]
# Suffixes of labelled runs (parameter sweeps of the reference scenario).
LABELS = [""]

# The HEOM reference lives in the canonical mkct4oqs checkout.
REF = os.environ.get("DMRG_REF_DATA", "../../mkct4oqs/source/output/data")
JSON = "./output/dmrg_error.json"
MD = "./output/dmrg_error.md"

def _corrmat(prefix, case, label=""):
    """``(t, C)`` of one tDMRG run, or ``None`` if the archive is absent."""
    d = "./output/data/" + prefix + "/"
    tag = f"dmrg_{case}" + (f"_{label}" if label else "")
    path = d + f"corrmat_{tag}.npz"
    if not os.path.exists(path):
        return None
    return load_group(path)[0][case]


def _ref(prefix):
    """``(t, C)`` of the HEOM reference correlation matrix."""
    path = os.path.join(REF, prefix) + os.sep + "corrmat_heom.npz"
    return load_group(path)[0]["qhm"]


def error(prefix, case, label=""):
    """Relative Frobenius error (percent); None if the run is absent."""
    got = _corrmat(prefix, case, label)
    if got is None:
        return None
    t, C = got
    tf = parse_xlims(load_config(prefix)["xlims"], reorder = True)[0][1]
    th, Ch = _ref(prefix)
    if th[-1] < tf * (1.0 - 1e-6):
        warnings.warn(
            f"[{prefix}] the HEOM reference reaches t={th[-1]:g} but the "
            f"metric window is [0, {tf:g}] (config xlims); the missing range "
            f"cannot be scored.", stacklevel = 2)
    win = t <= tf * (1.0 + 1e-9)
    if t[win][-1] < tf * (1.0 - 1e-6):
        warnings.warn(
            f"[{prefix}] dmrg_{case} reaches t={t[win][-1]:g} but the metric "
            f"window is [0, {tf:g}]; the overlap is scored and the window "
            f"recorded as {tf:g}.", stacklevel = 2)
    t, C = t[win], C[win]
    hw = th <= tf * (1.0 + 1e-9)
    Cr = core._interp_mat(Ch[hw], th[hw], t)
    return 100.0 * core._corrmat_error(C, Cr)


def main():
    ap = argparse.ArgumentParser(description = __doc__.splitlines()[0])
    ap.add_argument("--configs", nargs = "*", default = RUNLIST)
    ap.add_argument("--cases", nargs = "*", default = CASES)
    ap.add_argument("--labels", nargs = "*", default = LABELS)
    args = ap.parse_args()

    ensure_dir("./output")
    raw = {}
    for prefix in args.configs:
        rows = raw.setdefault(prefix, {})
        for label in args.labels:
            for c in args.cases:
                got = error(prefix, c, label)
                if got is None:
                    continue
                key = c if not label else f"{c}[{label}]"
                rows[key] = {"E": float(got),
                             "tf": float(parse_xlims(
                                 load_config(prefix)["xlims"],
                                 reorder = True)[0][1])}
    if os.path.exists(JSON):
        with open(JSON) as f:
            old = json.load(f)
        for p, rows in raw.items():          # merge, do not lose other runs
            old.setdefault(p, {}).update(rows)
        raw = old
    with open(JSON, "w") as f:
        json.dump(raw, f, indent = 2, sort_keys = True)

    keys = [c if not label else f"{c}[{label}]"
            for label in args.labels for c in args.cases]
    lines = ["# MPS/tDMRG accuracy against HEOM", "",
             "Metric: the global relative Frobenius error "
             "`E = ||C - C'||_F / ||C'||_F` over `[0, tf]` of the basis "
             "correlation matrix (all 16 complex elements), in percent.", "",
             "| config | " + " | ".join(keys) + " |",
             "|---|" + "---:|" * len(keys)]
    for prefix in args.configs:
        cells = []
        for k in keys:
            v = raw.get(prefix, {}).get(k)
            cells.append(f"{v['E']:.1f}" if v else "--")
        lines.append(f"| {prefix} | " + " | ".join(cells) + " |")
    with open(MD, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print("WROTE", JSON, "and", MD, flush = True)


if __name__ == "__main__":
    main()

# %%
