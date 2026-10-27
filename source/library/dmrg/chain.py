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
import os, pickle
import numpy as np

from tenpy.networks.site import SpinHalfSite
from tenpy.models.lattice import Lattice
from tenpy.models.model import CouplingModel, MPOModel
from tenpy.networks.mps import MPS
from corrmat import observables
from tenpy.algorithms.tdvp import TwoSiteTDVPEngine, SingleSiteTDVPEngine

from .model import _boson_site
from .tdmrg import _corrmat


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
    return MPS.from_Bflat(lat.mps_sites(), Bflat, bc='finite',
                          unit_cell_width=lat.mps_unit_cell_width,
                          dtype=complex)


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


def _evolve_one(b, psi_b, model_args, opts, t, snapdir, engine):
    """Evolve one of the two purified states, snapshotting it at every time.

    Runs in a joblib worker: it rebuilds the lattice/model (cheap and
    deterministic) so only the MPS has to cross the process boundary.  The
    observables are formed afterwards in the parent from the snapshot pairs
    (``M_01`` needs both states at the same time).
    """
    lat, model, model_prep, N = build_models(*model_args)
    Engine = SingleSiteTDVPEngine if engine == "onesite" else TwoSiteTDVPEngine
    eng = Engine(psi_b, model, opts)
    nT = len(t)
    snaps = []
    for it in range(nT):
        path = os.path.join(snapdir, f"b{b}_{it:05d}.pkl")
        with open(path, "wb") as f:
            pickle.dump(psi_b, f)
        snaps.append(path)
        if it < nT - 1:
            eng.run()
    return b, snaps


def run_observables(epsilon, Delta, omega, c, beta, t, quantum=True, anharm2=0.0,
                    nmax=6, chi=64, dt=None, dtau=0.01, verbose=False,
                    engine="twosite", expansion_dim=16,
                    t0=None, state=None, save_state=None, normalize_cmu=None,
                    n_jobs=1, chi_prep=None):
    """Chain-mapped observables (same dict as :func:`dmrg.tdmrg.run_observables`).

    ``engine`` selects the time-evolution engine: ``"twosite"`` (default) uses
    two-site TDVP; ``"onesite"`` uses single-site TDVP with Krylov subspace
    expansion (cheaper, ~3x faster here, but numerically unreliable).

    ``n_jobs > 1`` evolves the two initial states in parallel joblib workers
    (the two propagations are independent), snapshotting each MPS per output
    time so the cross overlap is formed in a short merge pass.  This roughly
    halves the wall time on a two-core machine.

    Segment support: if ``state`` exists it is a pickled ``{"psi", "time"}``
    checkpoint (the thermal preparation and the initial states are skipped);
    at the end the two MPS are written back to ``save_state``.
    ``normalize_cmu`` defaults to ``t0 == 0``.
    """
    t = np.asarray(t, dtype=float)
    lat, model, model_prep, N = build_models(epsilon, Delta, omega, c, nmax=nmax,
                                             anharm2=anharm2)

    dt = float(t[1] - t[0]) if dt is None else float(dt)
    if engine == "onesite":
        Engine = SingleSiteTDVPEngine
        # expansion_dim = 0 gives pure (fixed-chi) one-site TDVP, as in the
        # reference protocol once the bond dimension is established.
        lanczos = {'N_max': 20}
        if expansion_dim:
            lanczos['expansion_dim'] = expansion_dim
    else:
        Engine = TwoSiteTDVPEngine
        lanczos = {}

    if state and os.path.exists(state):
        with open(state, "rb") as f:
            payload = pickle.load(f)
        psi = payload["psi"]
        t0 = payload["time"]
    else:
        bath = prepare_thermal(lat, model_prep, beta,
                               chi=chi if chi_prep is None else chi_prep,
                               dtau=dtau, N=N, nmax=nmax)
        psi = {}
        for b in (0, 1):
            p = bath.copy()
            if b == 1:
                p.apply_local_op(0, 'Sigmax', unitary=True)
            psi[b] = p
    if t0 is None:
        t0 = float(t[0])
    if normalize_cmu is None:
        normalize_cmu = (t0 == 0.0)
    opts = {'dt': dt, 'N_steps': 1, 'preserve_norm': True,
            'trunc_params': {'chi_max': chi, 'svd_min': 1e-13},
            'start_time': t0, 'max_N_sites_per_ring': lat.N_sites,
            'lanczos_params': lanczos}
    nT = len(t)
    carr = np.empty((nT, 4, 4), dtype = complex)

    def store(it, C):
        carr[it] = C

    if n_jobs and n_jobs > 1:
        import shutil
        import tempfile
        from joblib import Parallel, delayed
        snapdir = tempfile.mkdtemp(prefix="dmrg_snap_")
        # Only two states are independent, so n_jobs <= 2; give each worker its
        # share of the remaining cores so all of them stay busy (BLAS threads).
        threads = int(os.environ.get(
            "DMRG_THREADS", max(1, (os.cpu_count() or 1) // max(1, n_jobs))))
        try:
            res = Parallel(n_jobs=n_jobs, inner_max_num_threads=threads)(
                delayed(_evolve_one)(
                    b, psi[b], (epsilon, Delta, omega, c, nmax, anharm2),
                    opts, t, snapdir, engine) for b in (0, 1))
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
        eng = {b: Engine(psi[b], model, opts) for b in (0, 1)}
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
