"""tDMRG time evolution and observables for the spin-boson model.

All five observables (populations, coherences, and auto-/cross-correlations)
are obtained from just *two* MPS evolutions.  Every initial (unnormalized)
system operator is rank one, ``sigma_x |0><0| = |1><0|`` and
``Y |chi><chi| = |Ychi><chi|``, and
``|chi>``, ``|Ychi>`` are linear combinations of the two system basis states.
Writing ``|phi_b(t)> = U(t) |b> |TFD>`` (``b = 0, 1``), one has

    pop_0   = <phi_0| P_0 |phi_0>,   pop_1 = <phi_0| P_1 |phi_0>,
    cohe_+  = <phi_0| sigma_+ |phi_0>,  cohe_- = <phi_0| sigma_- |phi_0>,
    C_mu(t) = <phi_0| sigma_x |phi_1>          (normalized to C_mu(0) = 1),
    C_x(t)  = <phi_chi| sigma_x |phi_{Ychi}>

and the last reduces to the four overlaps ``M_ab = <phi_a| sigma_x |phi_b>``
with fixed coefficients (see :data:`_CX`).
"""
import os, pickle
import numpy as np

from tenpy.algorithms.tdvp import TwoSiteTDVPEngine

from .model import build_model, init_mps

# Cx = sum_ab coeff_ab M_ab with M_ab = <phi_a| sigmax |phi_b>,
# from |chi> = cos(pi/3)|0> + sin(pi/3)|1> and |Ychi> = -i sin(pi/3)|0> + i cos(pi/3)|1>.
_S = np.sin(np.pi / 3.0)
_C = np.cos(np.pi / 3.0)
_CX = np.array([[-1j * _S * _C, 1j * _C ** 2],
                [-1j * _S ** 2, 1j * _S * _C]], dtype=complex)


def _local_expect(psi):
    return (psi.expectation_value('Sigmaz', sites=[0])[0],
            psi.expectation_value('Sigmax', sites=[0])[0],
            psi.expectation_value('Sigmay', sites=[0])[0])


def run_observables(epsilon, Delta, omega, c, beta, t, quantum=True, anharm2=0.0,
                    nmax=6, chi=128, dt=None, verbose=False,
                    t0=None, state=None, save_state=None, normalize_cmu=None):
    """Evolve the spin-boson MPS and return the five observables.

    Returns a dict with ``pop`` (nT, 2 real), ``cohe_p``, ``cohe_m``,
    ``Cmu`` (normalized) and ``Cx`` (nT complex).

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
    eng = {b: TwoSiteTDVPEngine(psi[b], model, opts) for b in (0, 1)}

    nT = len(t)
    pop = np.empty((nT, 2))
    cop = np.empty(nT, complex)
    com = np.empty(nT, complex)
    cmu = np.empty(nT, complex)
    cx = np.empty(nT, complex)
    for it in range(nT):
        sz0, sx0, sy0 = _local_expect(psi[0])
        sx1 = psi[1].expectation_value('Sigmax', sites=[0])[0]
        pop[it, 0] = (1.0 + sz0.real) / 2.0
        pop[it, 1] = (1.0 - sz0.real) / 2.0
        cop[it] = (sx0 + 1j * sy0) / 2.0
        com[it] = (sx0 - 1j * sy0) / 2.0
        tmp = psi[1].copy()
        tmp.apply_local_op(0, 'Sigmax', unitary=False)
        M01 = psi[0].overlap(tmp)
        M = np.array([[sx0, M01], [np.conj(M01), sx1]], dtype=complex)
        cmu[it] = M01
        cx[it] = np.sum(_CX * M)
        if verbose:
            print(f"  t={t[it]:.3f} chi={max(psi[0].chi)} "
                  f"pop0={pop[it, 0]:.4f}", flush=True)
        if it < nT - 1:
            for b in (0, 1):
                eng[b].run()
    if normalize_cmu and cmu[0] != 0:
        cmu = cmu / cmu[0]
    if save_state:
        with open(save_state, "wb") as f:
            pickle.dump({"psi": psi, "time": float(t[-1] + dt)}, f)
    return {'pop': pop, 'cohe_p': cop, 'cohe_m': com, 'Cmu': cmu, 'Cx': cx}
