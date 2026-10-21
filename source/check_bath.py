# %%
"""Validity check of the chain-mapped bath.

The chain represents the bath with ``Q = g_c q_0``, so two quantities have an
exact reference from the discrete modes::

    <Q^2>            = sum_j c_j^2 (1/2) coth(beta w_j / 2)
    Re <Q(t) Q(0)>   = sum_j c_j^2 (1/2) coth(beta w_j / 2) cos(w_j t)

This script prepares the thermal state, measures the static weight, evolves the
bath alone (no system) and compares the correlation.  It is the measured
validity of the tDMRG reference: with a truncated Fock space the weight falls
short, and with a small bond dimension the correlation is over-damped.

Usage::

    python check_bath.py [-N 16] [-nmax 16] [-chi 96] [-tf 1.5] [-dt 0.3]
"""
import argparse
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "./library"))
sys.pycache_prefix = os.path.abspath("./output/pycache")

import numpy as np
from tenpy.algorithms.tdvp import TwoSiteTDVPEngine

from dmrg import chain
from qutils import ohmic_modes

BETA = 0.6


def _weights(om):
    return 0.5 / np.tanh(BETA * om / 2.0)


def check(N=16, nmax=16, chi=96, tf=1.5, dt=0.3, dtau=0.02):
    om, cc = ohmic_modes(0.4, 4.0, N=N, wmax_factor=10)
    Om, kappa, gc = chain.chain_params(om, cc)
    lat, model_real, model_prep, NN = chain.build_models(
        1.0, 1.0, om, cc, nmax=nmax)
    bath = chain.prepare_thermal(lat, model_prep, BETA, chi=chi, dtau=dtau,
                                 N=NN, nmax=nmax)
    w = _weights(om)
    static_exact = float(np.sum(cc ** 2 * w))
    static = gc ** 2 * bath.expectation_value('q q', sites=[1])[0].real

    opts = {'dt': dt, 'N_steps': 1, 'preserve_norm': True,
            'trunc_params': {'chi_max': chi, 'svd_min': 1e-13},
            'start_time': 0.0, 'max_N_sites_per_ring': lat.N_sites}
    psi = bath.copy()
    phi = bath.copy()
    phi.apply_local_op(1, 'q', unitary=False)      # q |Psi(0)>, unnormalised
    eng_p = TwoSiteTDVPEngine(psi, model_prep, opts)
    eng_f = TwoSiteTDVPEngine(phi, model_prep, opts)
    t = np.linspace(0.0, tf, int(round(tf / dt)) + 1)
    corr = np.empty(len(t), complex)
    for it in range(len(t)):
        tmp = phi.copy()
        tmp.apply_local_op(1, 'q', unitary=False)
        corr[it] = gc ** 2 * psi.overlap(tmp)
        if it < len(t) - 1:
            eng_p.run()
            eng_f.run()
    exact = np.array([np.sum(cc ** 2 * w * np.cos(om * tt)) for tt in t])
    err = float(np.max(np.abs(corr.real - exact)) / np.max(np.abs(exact)))
    return static, static_exact, t, corr.real, exact, err


def main():
    p = argparse.ArgumentParser(description="Chain-bath validity check.")
    p.add_argument("-N", type=int, default=16)
    p.add_argument("-nmax", type=int, default=16)
    p.add_argument("-chi", type=int, default=96)
    p.add_argument("-tf", type=float, default=1.5)
    p.add_argument("-dt", type=float, default=0.3)
    a = p.parse_args()

    static, exact0, t, corr, ref, err = check(a.N, a.nmax, a.chi, a.tf, a.dt)
    print(f"N={a.N} nmax={a.nmax} chi={a.chi} tf={a.tf}")
    print(f"  static <Q^2> = {static:.3f} / {exact0:.3f} = {100*static/exact0:.1f}%")
    print(f"  correlation relative error = {err:.3f}")
    np.set_printoptions(precision=3, suppress=True)
    print("  t   ", np.round(t, 2))
    print("  DMRG", np.round(corr, 3))
    print("  exact", np.round(ref, 3))
    if static / exact0 < 0.9:
        print("FAIL: static coupling weight below 90% (raise -nmax)")
        sys.exit(1)
    if err > 0.35:
        print("WARN: correlation deviates; the reference is damped")


if __name__ == "__main__":
    main()

# %%
