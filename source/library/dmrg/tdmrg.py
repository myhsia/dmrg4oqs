"""tDMRG time evolution and the basis correlation matrix.

Two MPS evolutions suffice: writing ``|psi_b(t)> = U(t) |b> |TFD>``
(``b = 0, 1``), the basis correlation matrix in the Liouville basis
``|0><0|, |0><1|, |1><0|, |1><1|`` is

    C[:, j, k] = <psi_a| phi_j |psi_b>,   phi_k^dagger = |b><a|,

with ``(a, b) = (0,0), (0,1), (1,0), (1,1)`` for ``k = 0, 1, 2, 3`` (the
off-diagonal columns are the cross terms between the two states).  Every
reported observable is a contraction of ``C`` (``library/corrmat.py``); in
particular ``C_mu = C_{1,1} + C_{2,1}``.
"""
import os, pickle
import numpy as np

from corrmat import observables

from tenpy.algorithms.tdvp import TwoSiteTDVPEngine

from .model import build_model, init_mps

# System basis operators phi_j = |i><k| of the correlation matrix.
_PHI = (
    np.array([[1, 0], [0, 0]], dtype = np.complex128),   # |0><0|
    np.array([[0, 1], [0, 0]], dtype = np.complex128),   # |0><1|
    np.array([[0, 0], [1, 0]], dtype = np.complex128),   # |1><0|
    np.array([[0, 0], [0, 1]], dtype = np.complex128),   # |1><1|
)
# phi_j = c_I I + c_x sigma_x + c_y sigma_y + c_z sigma_z
_PHI_COEF = ((0.5, 0.0, 0.0, 0.5), (0.0, 0.5, 0.5j, 0.0),
             (0.0, 0.5, -0.5j, 0.0), (0.5, 0.0, 0.0, -0.5))
# Column k of C: <phi_a| . |phi_b> with (a, b) below.
_AB = ((0, 0), (0, 1), (1, 0), (1, 1))


def _corrmat(psi0, psi1):
    """Basis correlation matrix ``C_ij(t)`` from the two purified states.

    ``C[:, j, k] = <psi_a| phi_j |psi_b>`` with ``phi_k^dagger = |b><a|`` for
    ``k = 0..3``; the diagonal columns use the two states and the off-diagonal
    ones their cross overlaps (Hermiticity gives the conjugate direction).
    """
    X = {("I", 0, 0): 1.0, ("I", 1, 1): 1.0}
    for a, bra in ((0, psi0), (1, psi1)):
        for ax in "xyz":
            X[(ax, a, a)] = bra.expectation_value("Sigma" + ax, sites = [0])[0]
    for op in ("I", "x", "y", "z"):
        tmp = psi1.copy()
        if op != "I":
            tmp.apply_local_op(0, "Sigma" + op, unitary = False)
        X[(op, 0, 1)] = psi0.overlap(tmp)
        X[(op, 1, 0)] = np.conj(X[(op, 0, 1)])
    C = np.zeros((4, 4), dtype = np.complex128)
    for k, (a, b) in enumerate(_AB):
        for j, (cI, cx, cy, cz) in enumerate(_PHI_COEF):
            C[j, k] = (cI * X[("I", a, b)] + cx * X[("x", a, b)]
                       + cy * X[("y", a, b)] + cz * X[("z", a, b)])
    return C


def _evolve_one(b, psi_b, model_args, opts, t, snapdir):
    """Evolve one purified state, snapshotting it at every output time.

    Runs in a joblib worker; the model is rebuilt from ``model_args`` so only
    the MPS crosses the process boundary.
    """
    epsilon, Delta, omega, c, anharm2, nmax = model_args
    lat, model, H = build_model(epsilon, Delta, omega, c, anharm2=anharm2, nmax=nmax)
    eng = TwoSiteTDVPEngine(psi_b, model, opts)
    snaps = []
    for it in range(len(t)):
        path = os.path.join(snapdir, f"b{b}_{it:05d}.pkl")
        with open(path, "wb") as f:
            pickle.dump(psi_b, f)
        snaps.append(path)
        if it < len(t) - 1:
            eng.run()
    return b, snaps


def run_observables(epsilon, Delta, omega, c, beta, t, quantum=True, anharm2=0.0,
                    nmax=6, chi=128, dt=None, verbose=False,
                    t0=None, state=None, save_state=None, normalize_cmu=None,
                    n_jobs=1):
    """Evolve the spin-boson MPS and return the correlation matrix.

    Returns a dict with ``C`` (nT, 4, 4) and, for convenience, its five
    contractions ``pop`` (nT, 2 real), ``cohe_p``, ``cohe_m``, ``Cmu``
    (normalized) and ``Cx`` (nT complex), via ``library/corrmat.py``.

    Segment support: if ``state`` exists it is a pickled ``{"psi", "time"}``
    checkpoint that is resumed (the time is taken from the file); at the end
    the two MPS are written back to ``save_state``.  ``normalize_cmu`` defaults
    to ``t0 == 0`` so only a first segment rescales ``Cmu``.
    """
    t = np.asarray(t, dtype=float)
    lat, model, H = build_model(epsilon, Delta, omega, c, anharm2=anharm2, nmax=nmax)
    dt = float(t[1] - t[0]) if dt is None else float(dt)
    if state and os.path.exists(state):
        with open(state, "rb") as f:
            payload = pickle.load(f)
        psi = payload["psi"]
        t0 = payload["time"]
    else:
        psi = {b: init_mps(lat, omega, beta, nmax, sys_state=b, quantum=quantum)
               for b in (0, 1)}
    if t0 is None:
        t0 = float(t[0])
    if normalize_cmu is None:
        normalize_cmu = (t0 == 0.0)
    opts = {
        'dt': dt, 'N_steps': 1,
        'trunc_params': {'chi_max': chi, 'svd_min': 1e-13},
        'preserve_norm': True, 'start_time': t0,
        'max_N_sites_per_ring': lat.N_sites,
    }
    nT = len(t)
    carr = np.empty((nT, 4, 4), dtype = complex)

    def store(it, C):
        carr[it] = C

    if n_jobs and n_jobs > 1:
        import shutil
        import tempfile
        from joblib import Parallel, delayed
        snapdir = tempfile.mkdtemp(prefix="dmrg_snap_")
        threads = int(os.environ.get(
            "DMRG_THREADS", max(1, (os.cpu_count() or 1) // max(1, n_jobs))))
        try:
            res = Parallel(n_jobs=n_jobs, inner_max_num_threads=threads)(
                delayed(_evolve_one)(
                    b, psi[b], (epsilon, Delta, omega, c, anharm2, nmax),
                    opts, t, snapdir) for b in (0, 1))
            snaps = {b: s for (b, s) in res}
            for it in range(nT):
                with open(snaps[0][it], "rb") as f:
                    p0 = pickle.load(f)
                with open(snaps[1][it], "rb") as f:
                    p1 = pickle.load(f)
                store(it, _corrmat(p0, p1))
                if verbose:
                    print(f"  t={t[it]:.3f} chi={max(p0.chi)} "
                          f"pop0={carr[it, 0, 0].real:.4f}", flush=True)
            with open(snaps[0][-1], "rb") as f:
                final0 = pickle.load(f)
            with open(snaps[1][-1], "rb") as f:
                final1 = pickle.load(f)
            psi = {0: final0, 1: final1}
        finally:
            shutil.rmtree(snapdir, ignore_errors=True)
    else:
        eng = {b: TwoSiteTDVPEngine(psi[b], model, opts) for b in (0, 1)}
        for it in range(nT):
            store(it, _corrmat(psi[0], psi[1]))
            if verbose:
                print(f"  t={t[it]:.3f} chi={max(psi[0].chi)} "
                      f"pop0={carr[it, 0, 0].real:.4f}", flush=True)
            if it < nT - 1:
                for b in (0, 1):
                    eng[b].run()
    obs = observables(carr, normalize = bool(normalize_cmu))
    if save_state:
        os.makedirs(os.path.dirname(save_state) or ".", exist_ok=True)
        with open(save_state, "wb") as f:
            pickle.dump({"psi": psi, "time": float(t[-1] + dt)}, f)
    return {'C': carr, 'pop': obs['pop'].real, 'cohe_p': obs['cohe+'],
            'cohe_m': obs['cohe-'], 'Cmu': obs['corr'], 'Cx': obs['xorr']}
