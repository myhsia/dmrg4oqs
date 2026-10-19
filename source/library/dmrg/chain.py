"""Chain-mapped tDMRG reference (harmonic bath) with TDVP.

The star-shaped spin-boson bath (system coupled to every mode) is rotated to a
one-dimensional chain by a Lanczos tridiagonalisation of the (diagonal) matrix
of squared mode frequencies with the coupling vector as seed.  The system then
couples only to the *first* chain site and the bath becomes nearest-neighbour,
so the state has far less entanglement than in the original star geometry.

The finite-temperature bath is prepared as a thermofield double by
imaginary-time TDVP of the *physical* chain starting from the (truncated)
infinite-temperature purification ``sum_n |n>|n>``; the physical chain alone is
then evolved in real time.

Site ordering is interleaved, ``[system, phys_0, aux_0, phys_1, aux_1, ...]``,
so that each physical mode is adjacent to its purifying ancilla (keeping the
bond dimension small).  The physical chain couplings are therefore not
nearest-neighbour in the MPS geometry, so TDVP (not TEBD) is used throughout.

The quartic (anharmonic) term is nonlocal in the chain basis; it is included
as a *local* chain anharmonicity ``(a^2 Omega_k / 4) q_k^4`` as an
approximation (the chain coordinates are the effective normal modes).
"""
import numpy as np

from tenpy.networks.site import SpinHalfSite
from tenpy.models.lattice import Lattice
from tenpy.models.model import CouplingModel, MPOModel
from tenpy.networks.mps import MPS
from tenpy.algorithms.tdvp import TwoSiteTDVPEngine, SingleSiteTDVPEngine

from .model import _boson_site
from .tdmrg import _CX


def lanczos_chain(omega, c):
    """Lanczos tridiagonalisation of ``diag(omega^2)`` seeded by ``c*sqrt(omega)``.

    Returns ``(alpha, beta, g0)`` with ``alpha`` the tridiagonal diagonal,
    ``beta`` the off-diagonal, and ``g0 = |c*sqrt(omega)|`` the coupling of the
    system to the first chain site (in the unit-mass variables).
    """
    omega = np.asarray(omega, dtype=float)
    c = np.asarray(c, dtype=float)
    m = len(omega)
    Lam = omega ** 2
    g = c * np.sqrt(omega)
    g0 = float(np.linalg.norm(g))
    if g0 <= 0.0:
        raise ValueError("zero total coupling: cannot build a chain")
    V = np.zeros((m, m))
    alpha = np.zeros(m)
    beta = np.zeros(m - 1)
    V[:, 0] = g / g0
    b = 0.0
    vprev = np.zeros(m)
    v = V[:, 0].copy()
    for k in range(m):
        w = Lam * v - b * vprev
        a = float(v @ w)
        alpha[k] = a
        w = w - a * v - b * vprev
        for j in range(k + 1):                 # full reorthogonalisation
            w -= (V[:, j] @ w) * V[:, j]
        if k < m - 1:
            b = float(np.linalg.norm(w))
            beta[k] = b
            vprev = v
            v = w / b
            V[:, k + 1] = v
    return alpha, beta, g0


def chain_params(omega, c):
    """TeNPy chain parameters ``(Om, kappa, gc)``.

    ``Om`` are the chain frequencies, ``kappa`` the nearest-neighbour
    couplings, and ``gc`` the system coupling to the first chain site, all in
    the ``(w/2)(p^2 + q^2)`` convention of :mod:`dmrg.model`.
    """
    alpha, beta, g0 = lanczos_chain(omega, c)
    alpha = np.maximum(alpha, 1e-12)
    Om = np.sqrt(alpha)
    kappa = beta / (alpha[:-1] * alpha[1:]) ** 0.25
    gc = g0 / alpha[0] ** 0.25
    return Om, kappa, gc


def _phys(k):
    return 1 + 2 * k                 # physical chain site k (interleaved)


def _aux(k):
    return 2 + 2 * k                 # its purifying ancilla


def _nmax_list(nmax, N):
    """Per-chain-site Fock truncation; ``nmax`` may be an int or a length-``N`` sequence."""
    if np.isscalar(nmax):
        return [int(nmax)] * N
    s = [int(x) for x in nmax]
    if len(s) != N:
        raise ValueError(f"nmax sequence must have length N={N}, got {len(s)}")
    return s


def _make_lattice(Om, nmax):
    N = len(Om)
    nm = _nmax_list(nmax, N)
    sites = [SpinHalfSite(conserve=None)]
    for k in range(N):
        sites.append(_boson_site(nm[k]))
        sites.append(_boson_site(nm[k]))
    return Lattice([1], sites, bc_MPS='finite'), N


def build_models(epsilon, Delta, omega, c, nmax=6, anharm2=0.0):
    """Return ``(lat, model_real, model_prep, N)`` on one shared lattice."""
    Om, kappa, gc = chain_params(omega, c)
    lat, N = _make_lattice(Om, nmax)

    Mr = CouplingModel(lat)
    Mr.add_onsite_term(epsilon, 0, 'Sigmaz')
    Mr.add_onsite_term(Delta, 0, 'Sigmax')
    for k in range(N):
        Mr.add_onsite_term(Om[k], _phys(k), 'N')
        if anharm2:
            Mr.add_onsite_term(anharm2 * Om[k] / 4.0, _phys(k), 'q4')
    for k in range(N - 1):
        Mr.add_coupling_term(kappa[k], _phys(k), _phys(k + 1), 'q', 'q')
    Mr.add_coupling_term(gc, 0, _phys(0), 'Sigmaz', 'q')
    model_real = MPOModel(lat, Mr.calc_H_MPO())

    Mp = CouplingModel(lat)
    for k in range(N):
        Mp.add_onsite_term(Om[k], _phys(k), 'N')
        if anharm2:
            Mp.add_onsite_term(anharm2 * Om[k] / 4.0, _phys(k), 'q4')
    for k in range(N - 1):
        Mp.add_coupling_term(kappa[k], _phys(k), _phys(k + 1), 'q', 'q')
    model_prep = MPOModel(lat, Mp.calc_H_MPO())
    return lat, model_real, model_prep, N


def _infinite_T_mps(lat, N, nmax):
    """``prod_k sum_n |n>_phys,k |n>_aux,k`` (truncated identity purification)."""
    nm = _nmax_list(nmax, N)
    Bflat = [np.zeros((2, 1, 1), complex)]
    Bflat[0][0, 0, 0] = 1.0
    for k in range(N):
        d = nm[k] + 1
        Bp = np.zeros((d, 1, d), complex)
        for n in range(d):
            Bp[n, 0, n] = 1.0
        Ba = np.zeros((d, d, 1), complex)
        for n in range(d):
            Ba[n, n, 0] = 1.0
        Bflat += [Bp, Ba]
    return MPS.from_Bflat(lat.mps_sites(), Bflat, bc='finite', dtype=complex)


def prepare_thermal(lat, model_prep, beta, chi=64, dtau=0.01, N=1, nmax=6):
    """Thermofield-double bath state by imaginary-time TDVP.

    Starts from ``sum_n |n>|n>`` and evolves the *physical* chain by
    ``beta/2`` in imaginary time, giving ``sum_n exp(-beta E_n / 2)|n>|n>``,
    the purification of the Gibbs state.
    """
    nsites = lat.N_sites
    psi = _infinite_T_mps(lat, N, nmax)
    opts = {'dt': -1j * dtau, 'N_steps': 1, 'preserve_norm': True,
            'trunc_params': {'chi_max': chi, 'svd_min': 1e-13},
            'start_time': 0.0, 'max_N_sites_per_ring': nsites}
    eng = TwoSiteTDVPEngine(psi, model_prep, opts)
    eng.evolve(max(1, int(round(beta / 2.0 / dtau))), -1j * dtau)
    psi.canonical_form()
    psi.norm = 1.0            # canonical_form keeps the stored norm; reset it
    return psi


def run_observables(epsilon, Delta, omega, c, beta, t, quantum=True, anharm2=0.0,
                    nmax=6, chi=64, dt=None, dtau=0.01, verbose=False,
                    engine="twosite", expansion_dim=16):
    """Chain-mapped observables (same dict as :func:`dmrg.tdmrg.run_observables`).

    ``engine`` selects the time-evolution engine: ``"twosite"`` (default) uses
    two-site TDVP; ``"onesite"`` uses single-site TDVP with Krylov subspace
    expansion (cheaper, ~3x faster here, at the cost of some accuracy).
    """
    t = np.asarray(t, dtype=float)
    lat, model, model_prep, N = build_models(epsilon, Delta, omega, c, nmax=nmax,
                                             anharm2=anharm2)
    bath = prepare_thermal(lat, model_prep, beta, chi=chi, dtau=dtau,
                           N=N, nmax=nmax)

    dt = float(t[1] - t[0]) if dt is None else float(dt)
    if engine == "onesite":
        Engine = SingleSiteTDVPEngine
        lanczos = {'expansion_dim': expansion_dim, 'N_max': 20}
    else:
        Engine = TwoSiteTDVPEngine
        lanczos = {}
    opts = {'dt': dt, 'N_steps': 1, 'preserve_norm': True,
            'trunc_params': {'chi_max': chi, 'svd_min': 1e-13},
            'start_time': 0.0, 'max_N_sites_per_ring': lat.N_sites,
            'lanczos_params': lanczos}
    psi = {}
    eng = {}
    for b in (0, 1):
        p = bath.copy()
        if b == 1:
            p.apply_local_op(0, 'Sigmax', unitary=True)
        psi[b] = p
        eng[b] = Engine(psi[b], model, opts)

    nT = len(t)
    pop = np.empty((nT, 2))
    cop = np.empty(nT, complex)
    com = np.empty(nT, complex)
    cmu = np.empty(nT, complex)
    cx = np.empty(nT, complex)
    for it in range(nT):
        sz0 = psi[0].expectation_value('Sigmaz', sites=[0])[0]
        sx0 = psi[0].expectation_value('Sigmax', sites=[0])[0]
        sy0 = psi[0].expectation_value('Sigmay', sites=[0])[0]
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
    if cmu[0] != 0:
        cmu = cmu / cmu[0]
    return {'pop': pop, 'cohe_p': cop, 'cohe_m': com, 'Cmu': cmu, 'Cx': cx}
