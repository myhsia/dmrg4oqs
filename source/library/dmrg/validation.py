"""Regression tests for the tDMRG reference.

Run directly (from ``./source``)::

    python library/dmrg/validation.py

Checks the purified MPS evolution and the two-state observable extraction
against exact diagonalization for (i) a bare electronic system, (ii) a single
harmonic mode, and (iii) a single quartic (anharmonic) mode.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from dmrg import run_observables
from qutils import Pauli, basis, ladder_ops, qubit


# ---------------------------------------------------------------------------
# Exact references (harmonic initial bath; quartic only in the dynamics)
# ---------------------------------------------------------------------------
def exact_obs(epsilon, Delta, omega, c, beta, t, Nf=20, anharm2=0.0):
    sz, sx, sy = Pauli("z"), Pauli("x"), Pauli("y")
    n = len(omega)
    dim = Nf ** n
    H = np.zeros((2 * dim, 2 * dim), complex)
    qops = []
    for j in range(n):
        a, ad = ladder_ops(Nf)
        q = (a + ad) / np.sqrt(2.0)
        M = None
        for k in range(n):
            blk = q if k == j else np.eye(Nf)
            M = blk if M is None else np.kron(M, blk)
        qops.append(M)
        nj = None
        for k in range(n):
            blk = (ad @ a) if k == j else np.eye(Nf)
            nj = blk if nj is None else np.kron(nj, blk)
        H += np.kron(np.eye(2), (omega[j] * nj).astype(complex))
        if anharm2:
            H += np.kron(np.eye(2),
                         (anharm2 * omega[j] / 4.0 * np.linalg.matrix_power(q, 4)).astype(complex))
    Q = sum(c[j] * qops[j] for j in range(n))
    H += np.kron(epsilon * sz + Delta * sx, np.eye(dim)) + np.kron(sz, Q)
    pn = np.zeros(dim)
    for idx in range(dim):
        e = 0.0
        for j in range(n):
            m = (idx // Nf ** (n - 1 - j)) % Nf
            e += omega[j] * m
        pn[idx] = np.exp(-beta * e)
    pn /= pn.sum()
    rhoB = np.diag(pn)
    evals, evecs = np.linalg.eigh(H)
    th = 2 * np.pi / 3
    chi_v = qubit(theta = th).ravel()
    Ychi = Pauli("y") @ chi_v
    rho0 = basis(2, 0, kind = "ketbra")
    Sp, Sm = Pauli("+"), Pauli("-")

    def evolve(Aop, Bop):
        rho = np.kron(Bop, rhoB)
        rt0 = evecs.conj().T @ rho @ evecs
        out = np.empty(len(t), complex)
        for i, tt in enumerate(t):
            rt = rt0 * np.exp(-1j * (evals[:, None] - evals[None, :]) * tt)
            rt = evecs @ rt @ evecs.conj().T
            rs = np.zeros((2, 2), complex)
            for i1 in range(2):
                for i2 in range(2):
                    rs[i1, i2] = np.trace(
                        rt[i1 * dim:(i1 + 1) * dim, i2 * dim:(i2 + 1) * dim])
            out[i] = np.trace(Aop @ rs)
        return out

    pop = evolve(rho0, rho0)
    return dict(pop=np.column_stack((pop.real, 1 - pop.real)),
                cohe_p=evolve(Sp, rho0), cohe_m=evolve(Sm, rho0),
                Cmu=evolve(sx, sx @ rho0),
                Cx=evolve(sx, Ychi[:, None] * chi_v[None, :]))


def _cmp(res, ex, keys=("pop", "cohe_p", "cohe_m", "Cmu", "Cx")):
    worst = 0.0
    for k in keys:
        worst = max(worst, float(np.max(np.abs(res[k] - ex[k]))))
    return worst


def test_bare(tol=0.02, nmax=8, chi=32, dt=0.05):
    """Bare electronic system (uncoupled bath): exact two-state Rabi dynamics."""
    epsilon, Delta = 1.0, 1.0
    omega, c = np.array([4.0]), np.array([0.0])
    t = np.linspace(0.0, 8.0, int(8.0 / dt) + 1)
    ex = exact_obs(epsilon, Delta, omega, c, 0.6, t)
    res = run_observables(epsilon, Delta, omega, c, 0.6, t, nmax=nmax, chi=chi, dt=dt)
    res["Cmu"] = res["Cmu"] / res["Cmu"][0]
    ex["Cmu"] = ex["Cmu"] / ex["Cmu"][0]
    return _cmp(res, ex), tol


def test_single_mode(tol=0.02, nmax=12, chi=64, dt=0.05):
    """Single harmonic mode (thermalised)."""
    epsilon, Delta, beta = 1.0, 1.0, 0.6
    worst = 0.0
    for w0, c0 in [(4.0, 0.6), (4.0, 0.2), (0.5, 0.632)]:
        omega, c = np.array([w0]), np.array([c0])
        t = np.linspace(0.0, 4.0, int(4.0 / dt) + 1)
        ex = exact_obs(epsilon, Delta, omega, c, beta, t, Nf=24)
        res = run_observables(epsilon, Delta, omega, c, beta, t, nmax=nmax, chi=chi, dt=dt)
        res["Cmu"] = res["Cmu"] / res["Cmu"][0]
        ex["Cmu"] = ex["Cmu"] / ex["Cmu"][0]
        worst = max(worst, _cmp(res, ex))
    return worst, tol


def test_single_mode_anharmonic(tol=0.05, nmax=16, chi=64, dt=0.05):
    """Single quartic mode (anharmonic bath)."""
    epsilon, Delta, beta = 1.0, 1.0, 0.6
    a = 0.4 * np.sqrt(7.0 / 6.0)
    worst = 0.0
    for w0, c0 in [(4.0, 0.6), (4.0, 0.4)]:
        omega, c = np.array([w0]), np.array([c0])
        t = np.linspace(0.0, 4.0, int(4.0 / dt) + 1)
        ex = exact_obs(epsilon, Delta, omega, c, beta, t, Nf=28, anharm2=a ** 2)
        res = run_observables(epsilon, Delta, omega, c, beta, t, nmax=nmax, chi=chi,
                              dt=dt, anharm2=a ** 2)
        res["Cmu"] = res["Cmu"] / res["Cmu"][0]
        ex["Cmu"] = ex["Cmu"] / ex["Cmu"][0]
        worst = max(worst, _cmp(res, ex))
    return worst, tol


def test_chain_multi_mode(tol=0.05, nmax=10, chi=64, dt=0.1):
    """Chain-mapped TDVP against exact diagonalization (one/two modes)."""
    from dmrg import chain
    epsilon, Delta, beta = 1.0, 1.0, 0.6
    worst = 0.0
    for omega, c in [(np.array([4.0]), np.array([0.6])),
                     (np.array([2.0, 8.0]), np.array([0.5, 0.3]))]:
        t = np.linspace(0.0, 3.0, int(3.0 / dt) + 1)
        ex = exact_obs(epsilon, Delta, omega, c, beta, t, Nf=18)
        res = chain.run_observables(epsilon, Delta, omega, c, beta, t, nmax=nmax,
                                    chi=chi, dt=dt, dtau=0.005)
        res["Cmu"] = res["Cmu"] / res["Cmu"][0]
        ex["Cmu"] = ex["Cmu"] / ex["Cmu"][0]
        worst = max(worst, _cmp(res, ex))
    return worst, tol


TESTS = [
    ("bare electronic system (two-level Rabi)", test_bare),
    ("single harmonic mode", test_single_mode),
    ("single quartic mode (anharmonic bath)", test_single_mode_anharmonic),
    ("chain-mapped TDVP (one/two modes)", test_chain_multi_mode),
]


def main():
    failed = False
    for name, fn in TESTS:
        err, tol = fn()
        ok = err <= tol
        failed |= not ok
        print(f"[{'PASS' if ok else 'FAIL'}] {name}: max err = {err:.4f} "
              f"(tol {tol})")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
