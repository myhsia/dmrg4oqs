# %%
"""Accuracy of the MPS/tDMRG references against HEOM.

The metric is the one of the Padé-order scan (``library/mre.py``): over the
calculation window ``[0, tf]`` of each config, the relative ``L2`` error of six
curves

    {Re P0, Re P1, Re C_xx, Re C_xy, Im C_xx, Im C_xy}

against the HEOM reference, averaged over the six and reported in percent,
``(E / n_obs) * 100``.  ``tf`` is the config window (forced, for comparability
with the other references); a run whose series ends earlier is compared on the
overlap only and a warning is emitted.

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
import numpy as np
sys.path.append(os.path.join(os.path.dirname(__file__), "./library"))
sys.pycache_prefix = os.path.abspath("./output/pycache")
import mre as core
from fsutils import load_config, load_npz, ensure_dir
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

FLOOR = 5e-2        # ignore reference curves below this fraction of the largest


def _series(prefix, case, quantity, label=""):
    """``(t, y)`` of one observable for one (case, label)."""
    d = "./output/data/" + prefix + "/"
    tag = f"dmrg_{case}" + (f"_{label}" if label else "")
    return load_npz(d + f"{quantity}_{tag}.npz")


def _heom(prefix, quantity):
    d = os.path.join(REF, prefix) + os.sep
    return load_npz(d + f"{quantity}_heom.npz")


def error(prefix, case, label=""):
    """``(mean relative error %, n_curves)``; None if the run is absent."""
    d = "./output/data/" + prefix + "/"
    tag = f"dmrg_{case}" + (f"_{label}" if label else "")
    if not os.path.exists(d + f"pop_{tag}.npz"):
        return None
    tf = parse_xlims(load_config(prefix)["xlims"], reorder = True)[0][1]
    th, yh = _heom(prefix, "pop")
    _, yc = _heom(prefix, "corr")
    _, yx = _heom(prefix, "xorr")
    if th[-1] < tf * (1.0 - 1e-6):
        warnings.warn(
            f"[{prefix}] the HEOM reference reaches t={th[-1]:g} but the "
            f"metric window is [0, {tf:g}] (config xlims); the missing range "
            f"cannot be scored.", stacklevel = 2)

    t, pop = _series(prefix, case, "pop", label)
    _, c = _series(prefix, case, "corr", label)
    _, x = _series(prefix, case, "xorr", label)
    win = t <= tf * (1.0 + 1e-9)
    if t[win][-1] < tf * (1.0 - 1e-6):
        warnings.warn(
            f"[{prefix}] {tag} reaches t={t[win][-1]:g} but the metric window "
            f"is [0, {tf:g}]; the overlap is scored and the window recorded as "
            f"{tf:g}.", stacklevel = 2)
    t = t[win]
    pop = pop[win]
    c = c[win]
    x = x[win]

    hw = th <= tf * (1.0 + 1e-9)
    cmu = core._interp(yc[hw], th[hw], t)
    cx = core._interp(yx[hw], th[hw], t)
    pairs = [(pop[:, 0].real, np.interp(t, th[hw], yh[hw][:, 0].real)),
             (pop[:, 1].real, np.interp(t, th[hw], yh[hw][:, 1].real)),
             (c.real, cmu.real), (x.real, cx.real),
             (c.imag, cmu.imag), (x.imag, cx.imag)]
    # Curves that vanish by symmetry in the reference (e.g. Im C_xx for a
    # biased system) carry no information and would blow the relative error up
    # arbitrarily; drop those below FLOOR times the largest reference norm and
    # average over the rest.
    norms = [np.linalg.norm(b) for _, b in pairs]
    keep = [np.linalg.norm(b) >= FLOOR * max(norms) for _, b in pairs]
    if not any(keep):
        return None
    e = sum(core._rel(a, b) for (a, b), k in zip(pairs, keep) if k)
    return 100.0 * e / sum(keep), sum(keep)


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
                rows[key] = {"E": float(got[0]), "n_obs": int(got[1]),
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
             "Metric `(E/n_obs)*100 %` over `[0, tf]` for "
             "`{Re P0, Re P1, Re C_xx, Re C_xy, Im C_xx, Im C_xy}` "
             "(`n_obs` after dropping curves that vanish by symmetry).", "",
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
