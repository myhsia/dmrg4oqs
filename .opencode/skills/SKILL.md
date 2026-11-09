---
name: heom-dmrg
description: HEOM and MPS/tDMRG reference dynamics for the spin-boson model (harmonic and anharmonic bath, linear and quadratic coupling) — use when generating, validating or interpreting the exact reference correlation functions that the other methods are scored against
license: MIT
compatibility: opencode
metadata:
  audience: maintainers
  workflow: github
---

# HEOM and MPS/tDMRG reference dynamics

Project rules: this repository's `AGENTS.md`.  Shared, system-wide conventions (coding style,
git/PR, LaTeX, tooling) are in the global `repo-rules` skill.

## System and observables

The spin-boson model with the four cases `QHM` (quantum harmonic), `QAH` (quantum anharmonic),
`CHM` (classical harmonic) and `CAH` (classical anharmonic).  The anharmonic quartic coefficient is
`a = anharmon * sqrt(7/6)`.  The Ohmic bath density is
`J(omega) = pi*lambd/wC * omega * exp(-omega/wC)`, discretized with the shared Gauss-Legendre modes.
The **classical** cases use a per-mode squeezed thermal state: the classical bath lies below the
quantum zero-point and is not a Gibbs state.

The whole basis correlation matrix is obtained from just *two* MPS propagations.  With
`|psi_b(t)> = U(t)|b>|TFD>` (`b = 0, 1`) and the Liouville basis
`phi_j = |0><0|, |0><1|, |1><0|, |1><1|`:

    C[:, j, k] = <psi_a| phi_j |psi_b>,   phi_k^dagger = |b><a|,
    (a, b) = (0,0), (0,1), (1,0), (1,1)  for  k = 0, 1, 2, 3

(the diagonal columns are the state-0/1 expectations, the off-diagonal ones their cross
overlaps).  The five reported observables are contractions of `C`
(`source/library/corrmat.py`); in particular `C_mu = C_{1,1} + C_{2,1}`.

## HEOM

- Linear coupling: qutip's `OhmicEnvironment` + `HEOMSolver` (ESPIRA-I fit).
- Quadratic coupling: the extended DEOM in `source/library/extended_deom.py`.
- `source/reference.py` runs each config in its `RUNLIST` and writes only the basis correlation
  matrix `output/data/<prefix>/corrmat_heom.npz` (one `(t, C)` pair under `qhm`, the shared
  `save_group` schema, `C(0) = I`); the five observables are contractions of it
  (`source/library/corrmat.py`).
- **Only the `qhm` matrix is stored.**  It is the quantum-harmonic reference; the other three cases
  differ in the bath, so there is nothing to score them against (see "What the reference
  supports").  `dmrg.py` inherits the QHM-calibrated protocol for all four cases.

## MPS / tDMRG

**Star geometry is the production choice.**  The bath modes are independent, so the thermal state
is an exact product of per-mode thermofield doublets: no imaginary-time preparation, and the
truncation error at fixed `nmax` is the only thermal-state error.  The Hamiltonian is
`H_sys = epsilon*sz + Delta*sx` plus `V = sz` coupled to `Q = sum_j c_j q_j`, with each mode
harmonic (and quartic, in the anharmonic cases) on its own site.

**The chain geometry remains available** (`-geometry chain`): the same modes are rotated to a
one-dimensional chain by a Lanczos tridiagonalization with the coupling vector as seed, so the
system couples only to the first site; the finite-temperature bath is then a purification
prepared by imaginary-time TDVP.  The preparation dominates the cost (2 h+ at N=16, `chi=64`,
`nmax=16`) and the accuracy is no better, so it is kept for comparison, not for production.
`-geometry auto` selects chain for the quantum cases and star for the classical ones; pass
`-geometry star` explicitly for the calibrated protocol.

**Engines.**  Two-site TDVP is the reference engine.  The single-site *with subspace expansion*
variant (`-engine onesite`) is numerically wrong (it oscillates and collapses); do not use it.
Pure one-site TDVP at fixed `chi`, without the expansion, is fine.  `-engine krylov` builds the
global Krylov space `{psi, H psi, ..., H^{m-1} psi}` as MPS, exponentiates the `m x m` projection
exactly and truncates once; the step is exact to 1e-9 against ED on one mode, but at fixed `chi`
it is 8-19x worse than TDVP and *worsens* as `chi` grows (67.29 % / 72.44 % against 8.26 % /
3.76 % at `chi` = 16 / 32), so it is a validated reference implementation, not a production
engine.

**Outputs.**  `source/dmrg.py` writes the per-case `corrmat_dmrg_<case>[_<tag>].npz` and, after a
complete unsegmented `-case all` run, the grouped `corrmat_dmrg[_<tag>].npz` (one matrix per
case, `branch = component = ("00","01","10","11")`, `C(0) = I`).  The synced
`FIG0@corrmat.py` draws the 4x4 `C_ij(t)` pages as `corrmat-dmrg@<prefix>.pdf` and
`FIG4@robustness.py` the `<method>@<prefix>[_<label>].pdf` comparison, overlaying the
reference `corrmat_heom.npz` from the data root `./output/data` (`-ref DIR` for another root).

## What the reference supports

**Only QHM has a reference.**  The four cases differ in the *bath* (quantum/classical x
harmonic/anharmonic), while `corrmat_heom.npz` holds the quantum-harmonic matrix only.  Scoring
QAH/CHM/CAH against it measures the physical difference between the baths, not the numerical
error: no accuracy is claimed for them and their cross-case numbers are not to be quoted.  The
protocol is calibrated on QHM and inherited unchanged by the other three, which feed the
quantum-vs-classical comparison.

Score with the house metric (`dmrg_error.py`, the shared `library/mre.py`): the global relative
Frobenius error `||C - C'||_F / ||C'||_F` over all 16 complex elements, on the stated window, in
percent.

Reference protocol: star, N=16, `wmax_factor=10`, `chi=64`, `nmax=16`, `dt=0.025`, window `t <= 4`:

| `nmax` | E at `t <= 4` | E at `t <= 2` |
|---|---:|---:|
| 8 | 7.58 % | 3.75 % |
| 12 | 6.94 % | 2.55 % |
| **16** | **6.62 %** | **1.77 %** |

The `t <= 4` differences halve (0.64 then 0.32), so the limit is `approx 6.45 %`.  Same metric and
window, method comparison (MKCT archives in `mkct4oqs`, MASH/FSSH in `sh4oqs`): MKCT 0.54 %,
tDMRG 7.58 % (6.94 % at `nmax=12`), MASH 17.16 %, FSSH 37.21 %.

## Conventions

- **Order to certify a window: window -> `dt` -> `chi` -> `nmax` -> N.**
- **`dt` dominates.**  0.05 is 3.5x worse than 0.025 at `t <= 4` (26.70 % vs 7.58 %), while
  0.025 -> 0.0125 moves the score by 0.01 points (6.94 -> 6.93 % at `nmax=12`).  The improvement is
  per-step bond truncation, not integrator accuracy: on the exact state the exponential is exact to
  round-off at `dt = 0.05` and the Krylov error is 2.1e-5 (`library/dmrg/integrator.py`).  A Cayley
  step composes to `exp(-i t H_eff)` with `H_eff = 2 atan(H dt/2)/dt`, which anharmonicises the bath;
  it is the right scheme for imaginary time only.
- **`chi` is cheap**: 32 and 64 agree to 0.02 points at N=16, so 32 is enough there.  Raise `chi`
  only together with N.
- **Uniform `nmax`, never a weight-based per-mode schedule.**  `-nmax` accepts a per-mode sequence,
  but ordering it by static weight gives the two lowest modes (`<n> = 7.4, 1.1`) only two levels and
  corrupts the long-time memory (12.35 % against 3.75 % for uniform 8 at `t <= 2`).  Use uniform
  `nmax = 16` for an accuracy statement, 12 for cheaper production; 8 is badly unconverged.
- **N = 16, `wmax_factor = 10`.**  More modes at fixed `chi` are worse at every window (N=40:
  8.06 % at `t <= 4`), and at `t <= 1` the score is identical for `chi` = 32 and 64, so the cause is
  the discretisation, not truncation.  A long-window run starts from the reference protocol and
  raises N only with `chi` in step, timing a short segment first.
- **Window.**  The discrete bath correlation recurs at `t ≈ 2 pi N / w_max` (N=16 -> 2.5), but the
  measured growth is smooth (3.75 / 5.31 / 7.58 / 10.53 / 12.52 % at `t <= 2 / 3 / 4 / 5 / 6`), so
  `t = 4` is a legitimate reporting window; do not push it by adding modes.
- **Parallelism.**  `-nproc` runs the independent cases in joblib workers; `-njobs` (max 2) evolves
  the two initial states in separate workers, each writing MPS snapshots that the parent merges into
  the same archive (bit-identical to serial).  `DMRG_THREADS` caps the BLAS threads per worker
  (default `cpu_count // njobs`).  Memory, not cores, is usually the limit.
- **Segmentation.**  `-t0`/`-t1` with `-state <path>` resumes and writes back the pickled MPS pair and
  its evolution time, so a long run fits under the CI job cap.  Only a complete unsegmented set is
  grouped into `corrmat_dmrg.npz`.
- **Metric windows.**  A run made with `-t1 4` against a config `tf` of 12 is scored on the overlap
  `[0, 4]` and warns; `dmrg_error.py --windows 2 3 4` scores one long run on successive windows.

## Workflows and validation

```bash
cd source
python reference.py -prefix A-base_line        # HEOM corrmat -> output/data/A-base_line/
python dmrg.py A-base_line -geometry star -N 16 -nmax 16 -chi 64 -dt 0.025 \
       -t0 0 -t1 4 -nproc 4 -njobs 2           # four cases -> corrmat_dmrg_*
python dmrg_error.py --cases qhm --window 4     # Frobenius table -> output/dmrg_error.{json,md}
python FIG0@corrmat.py                          # corrmat-dmrg@<prefix>.pdf
python FIG4@robustness.py                       # dmrg@<prefix>.pdf
python check_bath.py -N 16 -nmax 16 -chi 96     # bath correlator against the analytic curve
python library/dmrg/validation.py               # regression suite (ED for N <= 3)
```

Validation criteria: `validation.py` reproduces exact diagonalization for N <= 3 to `~1e-2` (bare
two-level, single harmonic mode, single quartic mode, chain multimode, Krylov and propagator
checks), which shows the mapping and the observable algebra are sound; `check_bath.py` prints the
static `<Q^2>` fraction and the relative error of `Re<Q(t)Q(0)>` against
`sum_j c_j^2 (1/2) coth(beta w_j/2) cos(w_j t)`, which is where a bad thermal state shows up first.

It is also the guard for a silent contraction bug: `get_expval_ABrho0` must pair `Tr[A phi^dagger]`
with `phi[m, n]` (`einsum("mn,kmn->k", A, basis_list.conj())`).  A flipped index order conjugates
non-real operators and leaves the populations intact — only the coherence phases move — so the
observable suite is what catches it; a diffusion or population check alone does not.
