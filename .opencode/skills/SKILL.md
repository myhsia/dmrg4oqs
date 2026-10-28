---
name: heom-dmrg
description: HEOM and MPS/tDMRG reference dynamics for the spin-boson model (harmonic and anharmonic bath, linear and quadratic coupling)
license: MIT
compatibility: opencode
metadata:
  audience: maintainers
  workflow: github
---

# HEOM and MPS/tDMRG reference dynamics

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
- Quadratic coupling: the extended DEOM in
  `source/library/extended_deom.py`.
- `source/reference.py` runs each config in its `RUNLIST` and writes only the basis correlation
  matrix `output/data/<prefix>/corrmat_heom.npz` (one `(t, C)` pair under `qhm`, the shared
  `save_group` schema, `C(0) = I`); the five observables are contractions of it
  (`source/library/corrmat.py`).

## MPS / tDMRG

- The star-shaped bath is rotated to a one-dimensional chain by a Lanczos tridiagonalization
  (coupling vector as seed), so the system couples only to the first chain site.
- The finite-temperature bath is a thermofield-double purification prepared by imaginary-time TDVP.
- **Two-site TDVP is the correct engine.**  The single-site + subspace-expansion variant
  (`-engine onesite`) is numerically wrong (it oscillates/collapses); do not use it.
- `source/dmrg.py` writes the per-case `corrmat_dmrg_<case>.npz` and, for a complete
  `-case all` run, the grouped `corrmat_dmrg.npz`; the synced `FIG0@corrmat.py` /`
  `FIG4@robustness.py` draw the `corrmat-dmrg@<prefix>.pdf` and `dmrg@<prefix>.pdf` figures,
  overlaying the reference `corrmat_heom.npz` when it is pointed at with `-ref`.
- `library/dmrg/` holds `chain.py` (chain mapping, thermal state, TDVP), `model.py`, `tdmrg.py`
  (observable extraction) and `validation.py`.

## Conventions

- Run from `./source` with the `dmrg` conda environment active.
- The discrete bath correlation recurs at `t ≈ 2*pi*N/w_max`, so matching the continuum over
  `t <= 12` needs `N >= 80`; `nmax >= 8`, `chi >= 64` and `dt <= 0.3` are needed too.
  `-wmax_factor` sets `w_max = wmax_factor * wC`.
- Long runs are sharded with `-t0`/`-t1`/`-state`: a segment resumes the pickled MPS pair (and its
  evolution time) and writes it back, so each job stays under the 350-min GitHub job cap.  The
  two-site sweep is sequential and the two propagations cannot be split (`C_mu` needs both states
  at every output time), so more in-job parallelism does not help.
- Shared, system-wide conventions (repo etiquette, Python style, git/PR, LaTeX, tooling) are in the
  workspace `repo-rules` skill.

## Workflows

```bash
cd source
python reference.py                        # HEOM reference -> output/data/<prefix>/
python dmrg.py A-base_line            # four cases -> output/data/A-base_line/
python FIG4@robustness.py -group dmrg -ref ../../mkct4oqs/source/output/data
python library/dmrg/validation.py     # regression suite
```

CI: `.github/workflows/ci.yaml` (syntax, regression suite, HEOM import) and
`.github/workflows/dmrg.yaml` (manual matrix over `<config> x <case>`, artifacts
`dmrg-<config>-<case>`).
