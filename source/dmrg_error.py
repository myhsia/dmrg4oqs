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
    python dmrg_error.py --labels star_c64dt50 star_c64dt25 --plot   # + figure
"""
# %%
# 1. Requirements
import argparse, json, os, sys, warnings
import numpy as np
sys.path.append(os.path.join(os.path.dirname(__file__), "./library"))
sys.pycache_prefix = os.path.abspath("./output/pycache")
import mre as core
from fsutils import load_config, ensure_dir
from corrmat import load_corrmat
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

def _corrmat(prefix, case, label = ""):
    """``(t, C)`` of one tDMRG run (grouped or per-case), or ``None``."""
    got = load_corrmat("./output/data/" + prefix + "/", "dmrg", case, label)
    return got.get(case)


def _ref(prefix):
    """``(t, C)`` of the HEOM reference correlation matrix."""
    got = load_corrmat(os.path.join(REF, prefix), "heom")
    return got.get("qhm")


_WARNED = set()


def error(prefix, case, label = "", window = None):
    """Relative Frobenius error (percent); None if the run is absent.

    ``window`` overrides the config ``tf`` (the runs may be shorter, or a
    long run may be scored on successive windows ``t <= 2, 3, 4, ...``).
    """
    if case != "qhm" and case not in _WARNED:
        _WARNED.add(case)
        warnings.warn(
            f"[{prefix}] scoring case '{case}' against the quantum-harmonic HEOM "
            "matrix: the bath differs physically, so this number is a physical "
            "difference, NOT the method's accuracy.  Only qhm is an accuracy "
            "statement.", stacklevel = 2)
    got = _corrmat(prefix, case, label)
    if got is None:
        return None
    t, C = got
    tf = (parse_xlims(load_config(prefix)["xlims"], reorder = True)[0][1]
          if window is None else float(window))
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
    ap.add_argument("--window", type = float, default = None,
                    help = "score over [0, window] instead of the config tf")
    ap.add_argument("--windows", nargs = "*", type = float, default = None,
                    help = "score over several windows (e.g. 2 3 4) and print "
                           "one table per window")
    ap.add_argument("--plot", action = "store_true",
                    help = "also draw the run comparison for the first config "
                           "-> output/figures/dmrg_error@<config>.pdf")
    args = ap.parse_args()

    ensure_dir("./output")
    if args.windows:
        # successive windows: report each, then keep the largest for the JSON
        for w in args.windows:
            print(f"\n== window [0, {w:g}]", flush = True)
            raw = {}
            for prefix in args.configs:
                rows = raw.setdefault(prefix, {})
                for label in args.labels:
                    for c in args.cases:
                        got = error(prefix, c, label, window = w)
                        if got is not None:
                            rows[c if not label else f"{c}[{label}]"] = {
                                "E": float(got), "tf": float(w)}
            keys = [c if not label else f"{c}[{label}]"
                    for label in args.labels for c in args.cases]
            print("| config | " + " | ".join(keys) + " |")
            print("|---|" + "---:|" * len(keys))
            for prefix in args.configs:
                print(f"| {prefix} | " + " | ".join(
                    f"{raw[prefix][k]['E']:.2f}" if k in raw.get(prefix, {})
                    else "--" for k in keys) + " |")
        return

    raw = {}
    for prefix in args.configs:
        rows = raw.setdefault(prefix, {})
        for label in args.labels:
            for c in args.cases:
                got = error(prefix, c, label, window = args.window)
                if got is None:
                    continue
                key = c if not label else f"{c}[{label}]"
                rows[key] = {"E": float(got),
                             "tf": (float(parse_xlims(
                                 load_config(prefix)["xlims"],
                                 reorder = True)[0][1])
                                 if args.window is None else float(args.window))}
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
    if args.plot and args.configs:
        path = plot_comparison(args.configs[0], args.labels, args.cases, raw)
        print("WROTE", path, flush = True)


def plot_comparison(prefix, labels, cases, raw):
    """One page per case: the observables of every labelled run vs HEOM.

    The runs are the same method at different ``dt``/``chi``/``nmax``, which
    ``FIG0``/``FIG4`` plot one archive at a time; this panel puts them on one
    axis next to the reference.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from corrmat import observables

    d = "./output/data/" + prefix + "/"
    ref = load_corrmat(os.path.join(REF, prefix), "heom").get("qhm")
    ensure_dir("./output/figures")
    tag = [l or "base" for l in labels]
    if list(cases) != list(CASES):
        tag = list(cases) + tag            # record the cases when they are a subset
    path = ("./output/figures/dmrg_error@" + prefix
            + ("_" + "_".join(tag) if tag else "") + ".pdf")
    styles = ["-", "--", "-.", ":", (0, (3, 1, 1, 1))]
    with PdfPages(path) as pdf:
        for case in cases:
            fig, axes = plt.subplots(1, 3, figsize = (15, 4.2))
            tmax = 0.0
            for i, (label, ls) in enumerate(zip(labels, styles)):
                ser = load_corrmat(d, "dmrg", case = case, label = label)
                if not ser:
                    continue
                t, C = ser[case]
                tmax = max(tmax, float(t[-1]))
                obs = observables(C)
                lab = (label or "base") + f"  E={raw[prefix][case if not label else case + '[' + label + ']']['E']:.2f}%"
                axes[0].plot(t, obs["pop"][:, 0].real, linestyle = ls, lw = 1.6, label = lab)
                axes[1].plot(t, obs["corr"].real, linestyle = ls, lw = 1.6)
            if ref is not None and case == "qhm":
                th, Ch = ref
                o = observables(Ch)
                axes[0].plot(th, o["pop"][:, 0].real, "k-", lw = 1.0,
                             label = "HEOM")
                axes[1].plot(th, o["corr"].real, "k-", lw = 1.0)
            for ax, ylab in zip(axes[:2], (r"$P_0(t)$", r"$C_\mu(t)$")):
                ax.set_xlabel(r"$t$"); ax.set_ylabel(ylab)
                ax.grid(alpha = 0.3)
                ax.set_xlim(0.0, float(tmax))
            axes[0].legend(fontsize = 8, frameon = False)
            axes[0].set_title(f"{prefix} {case}: population")
            axes[1].set_title(f"{prefix} {case}: autocorrelation")
            # third panel: the metric
            vals = [raw[prefix].get(case if not l else f"{case}[{l}]", {}).get("E")
                    for l in labels]
            names = [l or "base" for l in labels]
            axes[2].bar(range(len(names)), [v or 0.0 for v in vals],
                        color = "C0", alpha = 0.75)
            axes[2].set_xticks(range(len(names)))
            axes[2].set_xticklabels(names, rotation = 30, fontsize = 7)
            axes[2].set_ylabel(r"$E$ (%)")
            axes[2].set_title("error vs HEOM")
            axes[2].grid(alpha = 0.3, axis = "y")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
    return path


if __name__ == "__main__":
    main()

# %%
