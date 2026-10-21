# %%
# 1. Requirements
import argparse
import numpy as np
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "./library"))
sys.pycache_prefix = os.path.abspath("./output/pycache")
from dmrg import chain
from dmrg.tdmrg import run_observables as run_observables_star
from qutils import ohmic_modes
from fsutils import load_config, save_time
from pltools import parse_xlims


# Run selection — comment out entries to skip.  ``chm``/``cah`` use the
# classical surrogate (no wavefunction representation).
CASES = ["qhm", "qah", "chm", "cah"]

# case -> (quantum sampling, anharmonic).
_CASE_MATRIX = {
    "qhm": (True,  False),
    "qah": (True,  True),
    "chm": (False, False),
    "cah": (False, True),
}


def _case_name(quantum: bool, anharm: bool) -> str:
    for name, (q, ah) in _CASE_MATRIX.items():
        if q == quantum and ah == anharm:
            return name
    raise ValueError(f"No case for quantum={quantum}, anharm={anharm}")


def _parse_args():
    parser = argparse.ArgumentParser(
        description = "MPS/tDMRG reference dynamics for the spin-boson model.")
    parser.add_argument("prefix", nargs = "?", default = None,
                        help = "Config prefix (default: run RUNLIST)")
    parser.add_argument("-case", choices = ["qhm", "qah", "chm", "cah", "all"],
                        default = "all", help = "Bath case")
    parser.add_argument("-geometry", choices = ["auto", "chain", "star"],
                        default = "auto",
                        help = "Bath geometry (auto: chain for quantum, star for classical)")
    parser.add_argument("-engine", choices = ["twosite", "onesite"],
                        default = "twosite",
                        help = "TDVP engine (onesite uses Krylov subspace expansion, ~3x faster)")
    parser.add_argument("-N", type = int, default = 20,
                        help = "Number of bath modes")
    parser.add_argument("-nmax", type = int, default = 6,
                        help = "Local Fock-space truncation per mode")
    parser.add_argument("-chi", type = int, default = 128,
                        help = "Maximum MPS bond dimension")
    parser.add_argument("-dt", type = float, default = 0.02,
                        help = "Output / TDVP time step")
    parser.add_argument("-dtau", type = float, default = 0.02,
                        help = "Imaginary-time step for the thermal-state preparation")
    parser.add_argument("-wmax_factor", type = float, default = 10.0,
                        help = "Quadrature cutoff, in units of wC")
    parser.add_argument("-tag", type = str, default = "",
                        help = "Optional suffix appended to the output tag (dmrg_<case>_<tag>)")
    parser.add_argument("-t0", type = float, default = None,
                        help = "Segment start time (default: 0)")
    parser.add_argument("-t1", type = float, default = None,
                        help = "Segment end time (default: the config tf)")
    parser.add_argument("-state", type = str, default = "",
                        help = "Checkpoint file: resume if it exists, save at the end")
    parser.add_argument("-njobs", type = int, default = 1,
                        help = "joblib workers for the two initial states within a case (1 = serial)")
    parser.add_argument("-nproc", type = int, default = 1,
                        help = "joblib processes for independent cases (local full runs)")
    return parser.parse_args()


def main(prefix: str, quantum: bool = True, anharm: bool = False,
         N: int = 20, nmax: int = 6, chi: int = 128, dt: float = 0.02,
         wmax_factor: float = 10.0, geometry: str = "chain", dtau: float = 0.02,
         engine: str = "twosite", tag_suffix: str = "",
         t0: float = None, t1: float = None, state: str = "",
         n_jobs: int = 1):
    # 3. Parameters (physical values; ``scale`` is a numerical device only)
    params = load_config(prefix)
    epsilon= params["epsilon"]
    Delta  = params["Delta"]
    lambd  = params["lambd"]
    wC     = params["wC"]
    beta   = params["beta"]
    a = float(params["anharmon"]) * np.sqrt(7.0 / 6.0) if anharm else 0.0
    xlims = parse_xlims(params["xlims"], reorder = True)
    tf    = xlims[0][1]
    t0    = 0.0 if t0 is None else float(t0)
    t1    = tf if t1 is None else float(t1)
    state = (state or "").strip("'\"")
    t     = np.linspace(t0, t1, int(round((t1 - t0) / dt)) + 1)
    omega, c = ohmic_modes(lambd, wC, N = N, wmax_factor = wmax_factor)

    case = _case_name(quantum, anharm)
    tag = f"dmrg_{case}" + (f"_{tag_suffix}" if tag_suffix else "")
    outdir = "./output/data/" + prefix
    geom = geometry if geometry in ("chain", "star") else (
        "chain" if quantum else "star")
    print(f"[{prefix}] {tag}: geom={geom} engine={engine} N={N} nmax={nmax} "
          f"chi={chi} tf={tf:g} beta={beta} wC={wC} lambd={lambd} a={a:.4f}",
          flush=True)

    # 4. Dynamics
    if geom == "chain":
        res = chain.run_observables(epsilon, Delta, omega, c, beta, t,
                                    quantum = quantum, anharm2 = a ** 2,
                                    nmax = nmax, chi = chi, dt = dt, dtau = dtau,
                                    engine = engine, t0 = t0,
                                    state = state or None,
                                    save_state = state or None,
                                    n_jobs = n_jobs)
    else:
        res = run_observables_star(epsilon, Delta, omega, c, beta, t,
                                   quantum = quantum, anharm2 = a ** 2,
                                   nmax = nmax, chi = chi, dt = dt, t0 = t0,
                                   state = state or None,
                                   save_state = state or None)

    # 5. Save to npz
    save_time(t, res["pop"],
              outdir, "pop_" + tag + ".npz", identifier = "P", compress = True)
    save_time(t, res["cohe_p"],
              outdir, "cohe_" + tag + ".npz", identifier = "S+", compress = True)
    save_time(t, res["cohe_m"],
              outdir, "cohem_" + tag + ".npz", identifier = "S-", compress = True)
    save_time(t, res["Cmu"],
              outdir, "corr_" + tag + ".npz",
              identifier = "⟨μ(t)μ(0)⟩", compress = True)
    save_time(t, res["Cx"],
              outdir, "xorr_" + tag + ".npz",
              identifier = "⟨Sx(t)Sy(0)⟩", compress = True)


RUNLIST = [
    "A-base_line",
]
def run(prefix: str, cases = CASES, N: int = 20, nmax: int = 6,
        chi: int = 128, dt: float = 0.02, wmax_factor: float = 10.0,
        geometry: str = "auto", dtau: float = 0.02, engine: str = "twosite",
        tag_suffix: str = "", t0: float = None, t1: float = None,
        state: str = "", nproc: int = 1, n_jobs: int = 1):
    def one(case):
        quantum, anharm = _CASE_MATRIX[case]
        main(prefix, quantum = quantum, anharm = anharm, N = N, nmax = nmax,
             chi = chi, dt = dt, wmax_factor = wmax_factor, geometry = geometry,
             dtau = dtau, engine = engine, tag_suffix = tag_suffix,
             t0 = t0, t1 = t1, state = state, n_jobs = n_jobs)
    # Case-level parallelism only makes sense for independent full runs: a
    # segmented run shares one checkpoint path and must stay serial.
    if nproc > 1 and len(cases) > 1 and not state and t0 is None and t1 is None:
        from joblib import Parallel, delayed
        Parallel(n_jobs = nproc)(delayed(one)(c) for c in cases)
    else:
        for case in cases:
            one(case)


if __name__ == "__main__":
    args = _parse_args()
    cases = CASES if args.case == "all" else [args.case]
    for prefix in ([args.prefix] if args.prefix else RUNLIST):
        run(prefix, cases = cases, N = args.N, nmax = args.nmax,
            chi = args.chi, dt = args.dt, wmax_factor = args.wmax_factor,
            geometry = args.geometry, dtau = args.dtau, engine = args.engine,
            tag_suffix = args.tag, t0 = args.t0, t1 = args.t1, state = args.state,
            nproc = args.nproc, n_jobs = args.njobs)

# %%
