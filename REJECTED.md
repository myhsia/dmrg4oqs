# Rejected methods, and what they taught

Two propagator-side experiments were built on top of this reference and are
recorded here rather than kept in the tree.  Both were plausible, both were
measured, and both are dead ends -- but one of them produced a *positive*
result about the existing engine that the validity section in `AGENTS.md`
depends on, so that result is written out in full below.

The rule is the same as in `mkct4oqs`: **a rejected method's code may be
deleted; its lesson may not.**  Nothing that was not authored here is deleted
on the strength of a failure recorded here.

Numbers are the house metric: the global relative Frobenius error
`||C - C'||_F / ||C'||_F` of the basis correlation matrix against the HEOM
reference, in percent, over the stated window.  All were measured on the
reference scenario (`A-base_line`, star geometry, N=16, `nmax=8`, `chi=64`,
`dt=0.025`, `t <= 2`) unless a row says otherwise.

---

## 1. The Cayley (Crank-Nicolson) step as a propagator

**What it was.**  Replace the Krylov exponential of the two-site effective
Hamiltonian by `(1 + i H dt/2)^{-1} (1 - i H dt/2)`: exactly unitary, second
order, A-stable, and the natural discretisation of the projected (Dirac-Frenkel)
Schrodinger equation that TDVP solves.  It was the suggestion that motivated the
whole study.

**Why it failed.**  A Cayley step is a function of `H`, so `n` of them compose to
exactly `exp(-i t H_eff)` with

    H_eff = 2 atan(H dt / 2) / dt

(verified to `2.3e-13` against a direct matrix exponential).  That map is
**nonlinear in the level index**: on a harmonic bath it does not shift the
frequencies uniformly, it anharmonicises them.  The leading shift at the star's
highest mode (`w = 40`) is `0.785 w` at `dt = 0.05` -- 21 % slow -- and matching
the Krylov accuracy needs `dt` roughly six times smaller, each step costing an
iterative linear solve.

**Measured** (free mode at `w = 40`, `t = 2`, error in `<q>`; every stepped run
ends at the same final time, so a correct propagator shows no `dt` dependence):

| dt | `w dt` | exact | Krylov (m=20) | Cayley |
|---:|---:|---:|---:|---:|
| 0.2 | 8.0 | 1.0e-14 | 4.5e-4 | 1.5e-1 |
| 0.1 | 4.0 | 1.0e-14 | 1.4e-3 | 1.6e-3 |
| 0.05 | 2.0 | 1.1e-14 | 2.1e-5 | 1.3 |
| 0.025 | 1.0 | 1.1e-14 | 8.7e-15 | 8.3e-1 |

A thermal two-mode star bath gives the same picture (`1.3e-14` Krylov against
`1.1e-2` Cayley at `dt = 0.05`).

**The positive result this produced, which stands.**  Because any propagator
built from `H` satisfies `U(dt)^n = U(n dt)`, a converged exponential cannot show
a `dt` error at fixed output times -- and the Krylov exponential here is exact to
round-off at `dt = 0.05`.  The documented `dt` dependence of the production
engine is therefore **per-step bond truncation**, not integrator accuracy.
Confirmed on the production pipeline, refining only `dt` (qhm, `nmax=8`,
`chi=64`):

| `dt` | E at `t <= 2` | E at `t <= 4` |
|---:|---:|---:|
| 0.05 | 9.00 % | 26.70 % |
| **0.025** | **3.75 %** | **7.58 %** |
| 0.0125 | 3.77 % | 6.93 % (at `nmax=12`) |

That is what fixed `dt` at 0.025 and raised the reference's accuracy by a factor
2.4; it also showed that the published "chi is not converged" table was a `dt`
artefact (`chi = 64` and `chi = 128` agree to `6e-6` at the converged step).

**Lesson.**  For real time, integrate with the Krylov exponential; Cayley belongs
to **imaginary time**, where the propagator is not unitary and A-stability is
what matters (`prepare_thermal`).

**Where the code is.**  Commit `53aae09`
(`source/library/dmrg/integrator.py`, which also carried the exact/Krylov/Cayley
comparison above).

---

## 2. The global Krylov engine on the MPS manifold

**What it was.**  A "global" alternative to the local two-site TDVP step: build
`K_m(H, psi) = span{psi, H psi, ..., H^{m-1} psi}` as MPS (each product
compressed to a working bond dimension), orthonormalise, exponentiate the
`m x m` projection exactly, and truncate **once** when forming the combination --
instead of truncating twice per bond per step.  The construction is the
vector-space analogue of the continued-fraction/moment object in `mkct4oqs`
(the Lanczos coefficients of the Krylov space are the Jacobi coefficients of the
moment sequence).

**Why it failed: this construction, not the idea.**  The implementation is
correct -- exact against diagonalisation on one mode (`8e-10`), and equal to
TDVP on a two-mode model (`5.7e-6`) where nothing is truncated.  At the reference
coupling it is 8-19x worse than TDVP and gets **worse as `chi` grows**:

| engine | `chi = 16` | `chi = 32` | `chi = 64` |
|---|---:|---:|---:|
| two-site TDVP | 8.26 % | 3.76 % | **3.75 %** |
| global Krylov (`m = 8`) | 67.29 % | 72.44 % | -- |

The naive explanation (truncation applied `m` times per step) predicts
convergence toward TDVP as `chi` grows; the measurement contradicts it.  What is
left is the repeated compression of the Krylov vectors themselves (at every
matvec *and* every addition), whose Gram deviates from the identity (`~1e-4` at
`chi = 16`) and feeds a generalised eigenproblem a poor basis.  The symptom is
artificial relaxation with the norm intact: population `0.095` at `t = 2` against
HEOM's `0.415`.

**Lesson.**  A fair test needs a *work* bond dimension (`chi_work >> chi_out`,
truncate once at the end) or variational compression of the Krylov vectors, at
`~m` times the cost per step.  At equal cost TDVP already reaches 3.76 % at
`chi = 32`, so the onus is on the variant.  Do not cite this result as evidence
against Krylov-space methods in general.

**Where the code is.**  Commit `53aae09` (`source/library/dmrg/krylov.py`, the
`-engine krylov` option of `dmrg.py`, and the two regression tests in
`validation.py` that were removed with it).
