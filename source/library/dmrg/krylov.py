"""Krylov (Lanczos) propagator on the MPS manifold -- the "global" counterpart
of the local two-site TDVP step.

Two-site TDVP exponentiates a *local* effective Hamiltonian and then truncates
the two-site tensor by an SVD, once per bond per step.  Every truncation
discards a little weight, i.e. a non-unitary damping of the modes the step has
just excited.  A Krylov step instead builds

    K_m(H, psi) = span{psi, H psi, ..., H^{m-1} psi}

as MPS, each product compressed back to a working bond dimension,
orthonormalises it, exponentiates the ``m x m`` projection of ``H`` exactly, and
truncates **once** when forming the combination.

It is the vector-space analogue of the continued fraction built from the moment
sequence in ``mkct_cayley``: the Lanczos coefficients of ``K_m`` are numerically
the Jacobi/CF coefficients of ``<A, L^n A>``, and the Cayley moments are the same
object in the Möbius-rotated variable.

Compression makes the Krylov vectors not exactly orthonormal, so the step keeps
the Gram matrix ``S_ij = <v_i|v_j>`` and solves the small generalised problem
``A c = lambda S c``; ``||S - 1||`` is returned as the step diagnostic.

The step deliberately allows a *work* bond dimension above the output one: the
Krylov construction is where the entanglement is generated, and truncating it
mid-construction is the failure mode this method is meant to avoid.
"""
import numpy as np


def _apply_H(H_MPO, psi, trunc):
    """``H psi`` (``MPO.apply`` is in place, so copy first).

    ``canonical_form(renormalize=False)`` moves the overall scale into
    :attr:`MPS.norm` instead of discarding it -- with the default
    ``renormalize=True`` the magnitude of ``H psi`` (and hence every Lanczos
    coefficient) would silently become 1.
    """
    out = psi.copy()
    H_MPO.apply(out, {'compression_method': 'SVD', 'trunc_params': trunc})
    out.canonical_form(renormalize = False)
    return out


def _scale(psi, c):
    """Scale in place, keeping the factor in ``norm`` (which ``overlap`` uses)."""
    psi.norm = psi.norm * c
    return psi


def _norm(psi):
    return float(np.sqrt(abs(psi.overlap(psi))))


def _unit(psi):
    """Copy scaled to unit norm."""
    out = psi.copy()
    return _scale(out, 1.0 / _norm(out))


def krylov_step(H_MPO, psi, dt, m = 10, trunc = None, reorth = True,
                return_info = False):
    """One global Krylov step ``psi -> exp(-i H dt) psi``.

    ``trunc`` is a tenpy ``trunc_params`` dict (``chi_max``, ``svd_min``) used
    both for the intermediate products and for the final compression.
    """
    trunc = dict(trunc or {'chi_max': 64, 'svd_min': 1e-13})

    v = [psi.copy()]
    v[0].canonical_form(renormalize = False)
    nrm = _norm(v[0])
    _unit(v[0])                                # |v_0> = psi / ||psi||

    alpha, beta = [], []
    cut = float(trunc.get('svd_min', 1e-13))

    def subtract(x, y, c):
        """``x - c y``, keeping the bond dimension at the working cutoff."""
        z = x.add(y, 1.0, -c, cutoff = cut)
        z.compress_svd(trunc)
        z.canonical_form(renormalize = False)
        return z

    w = _apply_H(H_MPO, v[0], trunc)
    for j in range(m):
        if j > 0:
            w = _apply_H(H_MPO, v[j], trunc)
        alpha.append(float(np.real(v[j].overlap(w))))
        # w <- w - alpha_j v_j - beta_{j-1} v_{j-1}
        w = subtract(w, v[j], alpha[j])
        if j > 0:
            w = subtract(w, v[j - 1], beta[j - 1])
        if reorth:                              # compression breaks orthogonality
            for k in range(j + 1):
                ov = v[k].overlap(w)
                if abs(ov) > 1e-13:
                    w = subtract(w, v[k], ov)
        nw = _norm(w)
        beta.append(nw)
        if nw < 1e-14 or j == m - 1:
            break
        _scale(w, 1.0 / nw)
        v.append(w)

    k = len(alpha)
    S = np.eye(k, dtype = complex)
    A = np.zeros((k, k), dtype = complex)
    for i in range(k):
        for j in range(i, k):
            S[i, j] = S[j, i] = v[i].overlap(v[j])
    for i in range(k):
        A[i, i] = alpha[i]
        if i + 1 < k:
            A[i, i + 1] = beta[i]
            A[i + 1, i] = beta[i]

    from scipy.linalg import eigh, expm
    dev = float(np.max(np.abs(S - np.eye(k))))
    if dev < 1e-10:
        c = expm(-1j * dt * A)[:, 0]
    else:
        lam, U = eigh(A, S)
        c = U @ (np.exp(-1j * dt * lam) * np.conj(U[0, :]))

    acc = _scale(v[0].copy(), c[0])
    for j in range(1, k):
        acc = acc.add(v[j], 1.0, c[j], cutoff = cut)
        acc.compress_svd(trunc)
        acc.canonical_form(renormalize = False)
    acc.canonical_form(renormalize = False)
    _scale(acc, nrm)
    acc.compress_svd(trunc)
    acc.canonical_form(renormalize = False)

    info = {"krylov_dim": int(k), "gram_dev": dev,
            "alpha": [float(a) for a in alpha],
            "beta": [float(b) for b in beta]}
    return (acc, info) if return_info else acc
