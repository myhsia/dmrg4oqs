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

All five observables (populations, the two coherences, the auto- and cross-correlation) are obtained
from just *two* MPS propagations, because every time-zero system operator is rank one.  With
`|phi_b(t)> = U(t)|b>|TFD>` (`b = 0, 1`):

    pop_0   = <phi_0| P_0 |phi_0>,      pop_1   = <phi_0| P_1 |phi_0>,
    cohe_+  = <phi_0| sigma_+ |phi_0>,  cohe_-  = <phi_0| sigma_- |phi_0>,
    C_mu(t) = <phi_0| sigma_x |phi_1>          (normalized to C_mu(0) = 1),
    C_x(t)  = <phi_chi| sigma_x |phi_{Ychi}>

where the last reduces to four overlaps with fixed coefficients.

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
- `source/dmrg.py` writes `output/data/<prefix>/*_dmrg_<case>.npz`;
  `source/plot_dmrg.py` writes `output/figures/dmrg@<prefix>.pdf` (legend bar, combined
  body, coherence body) and overlays the reference `corrmat_heom.npz` when present.
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
python plot_dmrg.py A-base_line       # -> output/figures/dmrg@A-base_line.pdf
python library/dmrg/validation.py     # regression suite
```

CI: `.github/workflows/ci.yaml` (syntax, regression suite, HEOM import) and
`.github/workflows/dmrg.yaml` (manual matrix over `<config> x <case>`, artifacts
`dmrg-<config>-<case>`).
