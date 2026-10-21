# %%
"""Plot the tDMRG reference dynamics for one config.

Reads ``output/data/<prefix>/{pop,cohe,cohem,corr,xorr}_dmrg_<case>.npz`` for
``qhm``/``qah``/``chm``/``cah`` and, when present, the HEOM reference
``*_heom.npz``.  Writes ``output/figures/dmrg@<prefix>.pdf``: the legend bar,
the combined population + Re/Im auto-/cross-correlation body, and the
coherence +/- body.

Run from this directory (the code root):  ``python plot_dmrg.py A-base_line``.
"""
import argparse
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "./library"))
sys.pycache_prefix = os.path.abspath("./output/pycache")

from matplotlib.backends.backend_pdf import PdfPages
from pltools import GridPlt, parse_xlims
from fsutils import load_config, ensure_dir, load_npz


VARIANTS = {
    "qhm":  ("C3", "-", 1.5, 1.0),
    "qah":  ("C0", ":", 1.5, 1.0),
    "chm":  ("C4", "-", 1.8, 0.5),
    "cah":  ("C2", ":", 1.8, 0.5),
    "heom": ("k",  "-", 1.2, 1.0),
}
CASES = ["qhm", "qah", "chm", "cah"]
YLIMS = {"corr": [-1, 1, 2], "xorr": [-1, 1, 2],
         "cohe+": [-.5, .5, 2], "cohe-": [-.5, .5, 2]}


def _suffix(name: str) -> str:
    return "heom" if name == "heom" else "dmrg_" + name


def _load(directory: str, names: list) -> dict:
    """Load each present observable per case (``pop`` as ``(t, y)``)."""
    def path(stem, n):
        return directory + stem + "_" + _suffix(n) + ".npz"

    def has(stem, n):
        return os.path.exists(path(stem, n))

    out = {}
    for n in names:
        obs = {}
        for stem, key in (("xorr", "xorr"), ("cohe", "cohe+"), ("cohem", "cohe-")):
            if has(stem, n):
                t, y = load_npz(path(stem, n))
                obs[key] = (t, y, YLIMS[key])
        if has("corr", n):
            t, y = load_npz(path("corr", n))
            obs["corr"] = (t, y, YLIMS["corr"])
        if has("pop", n):
            obs["pop"] = load_npz(path("pop", n))
        out[n] = obs
    return out


def main(prefix: str, cases = None) -> None:
    d = "./output/data/" + prefix + "/"
    wanted = CASES if not cases else [c for c in cases if c in CASES]
    names = [n for n in wanted if os.path.exists(d + "pop_" + _suffix(n) + ".npz")]
    if os.path.exists(d + "pop_heom.npz"):
        names.append("heom")
    cases = [n for n in names if n != "heom"]
    if not cases:
        raise FileNotFoundError("no DMRG data in " + d)
    obs = _load(d, names)

    def have(key):
        return all(key in obs[n] for n in cases)

    if not any(have(q) for q in ("corr", "xorr", "cohe+", "cohe-")):
        raise FileNotFoundError("no correlation data in " + d)

    params = load_config(prefix)
    xlims = parse_xlims(params["xlims"], reorder=True)

    def row(key):
        out = []
        for n in names:
            if key not in obs[n]:
                continue
            t, y, ylim = obs[n][key]
            out.append((t, y, *VARIANTS[n], ylim, r"\ttfamily " + n.upper()))
        return out

    pop_lines = None
    if have("pop"):
        pop_lines = []
        for n in names:
            if "pop" not in obs[n]:
                continue
            t, P = obs[n]["pop"]
            if n == "heom":
                color, ls, lw, alpha = VARIANTS[n]
                pop_lines.append((t, P[:, 0], color, ls, lw, alpha,
                                  r"\ttfamily " + n.upper()))
                pop_lines.append((t, P[:, 1], color, ls, lw, alpha / 2, ""))
                continue
            p, q = cases[cases.index(n) % 2], cases[2 + cases.index(n) % 2]
            color, ls = VARIANTS[n][0], VARIANTS[n][1]
            pop_lines.append((t, P[:, 0], color, ls, VARIANTS[p][2],
                              VARIANTS[p][3], r"\ttfamily " + n.upper()))
            pop_lines.append((t, P[:, 1], color, ls, VARIANTS[q][2],
                              VARIANTS[q][3], ""))

    rows = [row(q) for q in ("corr", "xorr") if have(q)]
    biy = [r"$\expval*{\sigma_x(t) \sigma_x(0)}$",
           r"$\expval*{\sigma_x(t) \sigma_y(0)}$"][:len(rows)]

    ensure_dir("./output/figures")
    with PdfPages("./output/figures/dmrg@" + prefix + ".pdf") as pdf:
        GridPlt.reim4gridplt(
            biylabel=biy, supxlabel="Time", xlims=xlims,
            pop_lines=pop_lines, plot_data=rows, pdf=pdf, close_fig=True)
        if have("cohe+") and have("cohe-"):
            GridPlt.reim4gridplt(
                biylabel=[r"$\expval*{\sigma_+(t)}$", r"$\expval*{\sigma_-(t)}$"],
                supxlabel="Time", xlims=xlims,
                plot_data=[row("cohe+"), row("cohe-")], pdf=pdf, close_fig=True)


RUNLIST = ["A-base_line"]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot the tDMRG reference.")
    parser.add_argument("prefix", nargs="?", default=None)
    parser.add_argument("-cases", nargs="+", default=None,
                        help="Cases to draw (default: every case with data)")
    args = parser.parse_args()
    for p in ([args.prefix] if args.prefix else RUNLIST):
        main(p, cases=args.cases)

# %%
