"""Spin-boson MPS model and purified (thermofield) initial state.

Site ordering::

    [ system ] [ phys_1 ] [ aux_1 ] [ phys_2 ] [ aux_2 ] ...

The system site is a ``SpinHalfSite``; each bath mode is a truncated
``BosonSite`` (physical) followed by a second, inert ``BosonSite`` (auxiliary)
that purifies the finite-temperature bath.  The auxiliary sites carry no
Hamiltonian; tracing them out gives the physical density matrix, so
``<psi| A_sys |psi> = Tr[A rho]``.

Hamiltonian (matching the project's ``nohalf=True`` convention)::

    H = epsilon sigmaz + Delta sigmax
        + sum_j w_j (b_j^dag b_j + 1/2)
        + sigmaz * sum_j c_j q_j
        + sum_j (a^2 w_j / 4) q_j^4,      q_j = (b_j + b_j^dag)/sqrt(2),

where ``a^2`` is passed as ``anharm2`` (the quartic coefficient squared,
``a = a' * sqrt(7/6)``).
"""
import os, sys
import numpy as np
from scipy.linalg import expm

from tenpy.networks.site import SpinHalfSite, BosonSite
from tenpy.models.lattice import Lattice
from tenpy.models.model import CouplingModel, MPOModel
from tenpy.networks.mps import MPS

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from qutils import thermal_weights, ladder_ops, thermal_variances


def _phys_idx(j):
    """MPS index of the physical site of bath mode ``j``."""
    return 1 + 2 * j


def _anc_idx(j):
    """MPS index of the auxiliary (purification) site of bath mode ``j``."""
    return 2 + 2 * j


def _boson_site(nmax, with_ops=True):
    b = BosonSite(Nmax=nmax, conserve=None)
    if with_ops:
        B = b.B.to_ndarray()
        Bd = b.Bd.to_ndarray()
        q = (B + Bd) / np.sqrt(2.0)
        b.add_op('q', q)
        b.add_op('q4', q @ q @ q @ q)
    return b


def _nmax_list(nmax, N):
    """Per-mode Fock truncation; ``nmax`` may be an int or a length-``N`` sequence."""
    if np.isscalar(nmax):
        return [int(nmax)] * N
    s = [int(x) for x in nmax]
    if len(s) != N:
        raise ValueError(f"nmax sequence must have length N={N}, got {len(s)}")
    return s


def make_sites(omega, nmax):
    nm = _nmax_list(nmax, len(omega))
    sites = [SpinHalfSite(conserve=None)]
    for j, _ in enumerate(omega):
        sites.append(_boson_site(nm[j], with_ops=True))
        sites.append(_boson_site(nm[j], with_ops=False))
    return sites


def build_model(epsilon, Delta, omega, c, anharm2=0.0, nmax=6):
    """Return ``(lat, model, H_MPO)`` for the spin-boson star model."""
    sites = make_sites(omega, nmax)
    lat = Lattice([1], sites, bc_MPS='finite')
    M = CouplingModel(lat)
    M.add_onsite_term(epsilon, 0, 'Sigmaz')
    M.add_onsite_term(Delta,   0, 'Sigmax')
    for j in range(len(omega)):
        ip = _phys_idx(j)
        M.add_onsite_term(omega[j], ip, 'N')
        if anharm2:
            M.add_onsite_term(anharm2 * omega[j] / 4.0, ip, 'q4')
        M.add_coupling_term(c[j], 0, ip, 'Sigmaz', 'q')
    H = M.calc_H_MPO()
    model = MPOModel(lat, H)
    return lat, model, H


def _squeeze_matrix(nmax, r):
    """Squeezing operator ``exp[(r/2)(a^2 - a^dag^2)]`` in the Fock basis."""
    a, ad = ladder_ops(nmax + 1)
    return expm(0.5 * r * (a @ a - ad @ ad))


def _pair_tensor(w, beta, nmax, quantum=True):
    """Two-site (physical x auxiliary) purification tensor of one mode.

    ``quantum``: harmonic Gibbs weights, physical variance ``(1/2)coth(beta w/2)``.
    classical surrogate: squeeze the thermal state so that the physical
    position variance matches the classical value ``1/(beta w)`` (the classical
    bath has no zero-point fluctuations).  The classical variance is always
    below the quantum one, so the squeezing parameter is well defined.
    """
    p = thermal_weights(w, beta, nmax)
    if quantum:
        return np.diag(np.sqrt(p))
    a, ad = ladder_ops(nmax + 1)
    q = (a + ad) / np.sqrt(2.0)

    def varq(r):
        S = _squeeze_matrix(nmax, r)
        B = S @ np.diag(np.sqrt(p))
        rho = B @ B.conj().T
        rho = rho / np.trace(rho)
        return float(np.real(np.trace(q @ q @ rho)))

    target = float(thermal_variances(w, beta, quantum=False))
    grid = np.linspace(0.0, 2.0, 41)
    vals = np.array([varq(r) for r in grid])
    below = np.nonzero(vals <= target)[0]
    if len(below) == 0:                       # cannot reach: clamp to closest
        r = float(grid[np.argmin(np.abs(vals - target))])
    else:
        k = below[0]
        lo, hi = (0.0 if k == 0 else float(grid[k - 1])), float(grid[k])
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if varq(mid) <= target:
                hi = mid
            else:
                lo = mid
        r = 0.5 * (lo + hi)
    B = _squeeze_matrix(nmax, r) @ np.diag(np.sqrt(p))
    return B / np.sqrt(np.trace(B @ B.conj().T))


def init_mps(lat, omega, beta, nmax, sys_state=0, quantum=True):
    """Purified initial state ``|sys_state> (x) TFD(bath)``.

    ``|sys_state>`` is the system basis state (0 or 1).  The bath is a product
    over modes of ``sum_n sqrt(p_n) |n>_phys |n>_aux`` (the thermofield double
    of the (possibly squeezed) thermal state).
    """
    Bflat = []
    bs = np.zeros((2, 1, 1), complex)
    bs[sys_state, 0, 0] = 1.0
    Bflat.append(bs)
    nm = _nmax_list(nmax, len(omega))
    for j in range(len(omega)):
        n_j = nm[j]
        pair = _pair_tensor(omega[j], beta, n_j, quantum=quantum)
        Bp = np.zeros((n_j + 1, 1, n_j + 1), complex)
        Bp[:, 0, :] = pair
        Ba = np.zeros((n_j + 1, n_j + 1, 1), complex)
        for n in range(n_j + 1):
            Ba[n, n, 0] = 1.0
        Bflat.append(Bp)
        Bflat.append(Ba)
    return MPS.from_Bflat(lat.mps_sites(), Bflat, bc='finite',
                          unit_cell_width=lat.mps_unit_cell_width,
                          dtype=complex)
