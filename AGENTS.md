# Guard Rules

The following directories, files, and patterns MUST NOT be read or accessed:

- `.git/`
- `**/__pycache__/`

Any tool that reads files (`read`, `glob`, `grep`, `shell` with file-reading commands) must respect
these restrictions.

Scientific background and methodology are described in
[`.opencode/skills/SKILL.md`](./.opencode/skills/SKILL.md); load the `heom-dmrg` skill for the physics.

## Layout

All code is under `source/`, the working directory for the drivers (paths such as `./config` and
`./output` are relative to it):

- `source/dmrg.py` — driver: for one config it runs the bath cases and writes the basis
  correlation matrix per case (`corrmat_dmrg_<case>[-<tag>].npz`) and, after a complete
  `-case all` run, the grouped `corrmat_dmrg.npz` (one matrix per case,
  `branch = component = ("00","01","10","11")`, `C(0) = I`).
- `source/FIG0@corrmat.py`, `source/FIG4@robustness.py` — synced copies of the canonical
  figure scripts: `python FIG0@corrmat.py -ref ../../mkct4oqs/source/output/data` draws the
  4x4 `C_ij(t)` pages of every method found here and
  `python FIG4@robustness.py -ref ../../mkct4oqs/source/output/data` the
  `dmrg@<prefix>.pdf` comparison figure (replacing the retired `plot_dmrg.py`).
- `source/reference.py` — reference dynamics (qutip `HEOMSolver` for linear coupling; the extended
  DEOM, i.e. the dissipaton hierarchy, in `library/extended_deom.py` for quadratic coupling), writing
  only the basis correlation matrix `output/data/<prefix>/corrmat_heom.npz` (one `(t, C)` pair
  under `qhm`); the five observables are contractions of it (`library/corrmat.py`).
- `source/library/extended_deom.py` — extended DEOM (dissipaton) for polynomial coupling.  It is the
  dissipaton equation of motion, not the standard HEOM: the hierarchy is stored in the
  square-root-ladder normalisation so that magnitude filtering (`ferr`) can prune DDOs, the generator
  is assembled from compressed edge lists, and the march is adaptive in chunks of output points.
  `a2 = 0` reduces it to the linear HEOM.
- `source/library/dmrg/` — `chain.py` (star→Lanczos chain, two-site TDVP, thermofield-double bath),
  `model.py` (star MPS, coupling MPO, per-mode squeezed classical state), `tdmrg.py` (shared
  observable extraction), `validation.py` (regression suite).
- `source/environment.yaml`, `source/requirements.txt` — the `dmrg` conda environment.

**Canonical archives** are `corrmat_<method>[_<label>].npz`.  The per-curve
`{pop,cohe,cohem,corr,xorr}_<method>_<case>.npz` files are legacy from the pre-corrmat workflow;
no script reads them, and they are archived out of `output/data/<prefix>/`.

## Method

- **HEOM / extended DEOM.**  Linear coupling uses qutip's `OhmicEnvironment` + `HEOMSolver` (ESPIRA-I fit); quadratic
  coupling uses the extended DEOM (`extended_deom.py`).  `reference.py` runs each config in its
  `RUNLIST` and writes the single `corrmat_heom.npz` (the harmonic reference shared by the four
  cases).
- **MPS / tDMRG.**  The star-shaped Ohmic bath (shared Gauss-Legendre discretization) is rotated to a
  one-dimensional chain by a Lanczos tridiagonalization (coupling vector as seed), so the system
  couples only to the first chain site; the finite-temperature bath is a thermofield-double
  purification prepared by imaginary-time TDVP.  The whole basis correlation matrix follows from
  two propagations `|psi_b(t)> = U|b>|TFD>` (b = 0, 1): the diagonal columns are the state-0/1
  expectations and the off-diagonal ones the cross overlaps (see `tdmrg._corrmat`); the five
  reported observables are its contractions.

## Validity (measured, not assumed)

Score with the house metric (`dmrg_error.py`, the shared `library/mre.py`): the global
relative Frobenius error `||C - C'||_F / ||C'||_F` of the basis correlation matrix (all 16
complex elements) against the HEOM reference matrix over the config window, in percent.
Every DMRG value comes from the corrmat workflow (`dmrg.py` -> `corrmat_dmrg_<case>.npz` ->
`corrmat.observables`); no stored curves are involved.

**Reference scenario (star, N=16, `nmax=8`, chi=64, `dt=0.05`, t <= 2):**

| case | E |
|---|---:|
| QHM | **9.0 %** |
| QAH | 19.9 % |
| CHM | 41.3 % |
| CAH | 67.1 % |

Same window and metric, other references (common grid): MKCT 0.3 / 1.6 / 23.7 / 24.8,
MASH 20.0 / 14.7 / 13.8 / 10.6, FSSH 21.1 / 20.7 / 18.7 / 19.4 (QHM / QAH / CHM / CAH).
So tDMRG beats FSSH on QHM, ties it on QAH, and is poor on the classical surrogate.

**The decisive parameter is `dt`, not chi.**  The star's highest mode is `w_max = 40`, so
`dt=0.3` means `w dt ~ 12`; the per-step truncation then damps the fast modes and fabricates
a "too fast initial decay" that looks like a chi/`nmax`/thermal-state problem.  Measured
(N=16, chi=64, `nmax=8`, t <= 2):

| `dt` | QHM | QAH | CHM | CAH |
|---|---:|---:|---:|---:|
| 0.15 | 15.6 | 85.9 | 81.9 | 86.2 |
| **0.05** | **9.0** | **19.9** | **41.3** | **67.1** |

`dt` is converged over this range and the Fock space is converged (`nmax` 8 -> 12 moves the
score by ~1 point at chi=128), but **chi is not certified**: the chi=128 corrmat run is the
open check.  Quote the reference as the most accurate configuration found, with that caveat.

**Geometry.**  Prefer the **star** with a **per-mode `nmax`** (`-nmax` accepts a sequence): the
star's thermal state is an exact product of per-mode TFDs (no imaginary-time preparation, which
dominated the chain runs at 2 h+), and a schedule such as `[20, 8, 4, 2, ...]` keeps 99.4 % of
the static weight `<Q^2>` at 7.5x lower cost (`sum d^3` 10 466 vs 78 608 for uniform 16).  The
chain remains available but its preparation is expensive and no better.

**Window.**  With N modes the discrete bath correlation recurs at `t ~ 2 pi N / w_max`; for
N=16, `w_max=40` that is ~2.5, so the reference is meaningful only up to `t ~ 2` (the metric
window is still the config's `tf`; `dmrg_error.py` warns when the series ends earlier).

Pushing the recurrence out is **not** the way to extend it: at fixed chi, more modes add
entanglement that is truncated harder, so the curve degrades *before* its own recurrence.
A converged 0-12 reference needs N >~ 80 **and** chi >~ 256 together, which is out of reach
here (an intermediate N=40, chi=128, `dt=0.05` attempt was ~40 h and was abandoned).  The
reference stays a short-time statement.

**Classical cases.**  `chm`/`cah` relax far faster than HEOM at every setting: the per-mode
squeezed thermal state is not a Gibbs state, so they are not a quantitative reference.

**Regression.**  `library/dmrg/validation.py` reproduces exact diagonalization for N <= 3 to
~1e-2, so the mapping and the observable algebra are sound; `check_bath.py` measures the bath
itself against `sum_j c_j^2 (1/2) coth(beta w_j/2) cos(w_j t)`.


## Shared pieces (mirrored from `mkct4oqs`)

The **canonical** `source/reference.py`, `source/config/`, and
`source/library/{qutils,fsutils,pltools,extended_deom,mre,corrmat}.py` live in the private `mkct4oqs`
repository.  This
repository keeps **synced copies** so that a public clone and the CI work without a token.

> **Sync reminder.**  After changing any of them in `mkct4oqs`, run
> `bash tools/sync_from_mkct.sh` here, then commit and push.  `sh4oqs` mirrors the same pieces
> and needs the same sync.

`fsutils.py` also carries the grouped-archive writer/loader (`save_group`/`load_group`) used by
the reference (`corrmat_heom.npz`) and the companion's MKCT results; `corrmat.py` contracts the
matrix into the five observables.

## Conventions

- Run everything from `source/` with the `dmrg` environment active.
- Four cases: `QHM`, `QAH`, `CHM`, `CAH`.  The anharmonic quartic coefficient is
  `a = anharmon * sqrt(7/6)`.
- The **classical** cases use a per-mode squeezed thermal state (the classical bath lies below the
  quantum zero-point and is not a Gibbs state).
- **`engine`**: two-site TDVP is the correct engine.  The single-site + subspace-expansion variant
  (`-engine onesite`) is numerically wrong (it oscillates/collapses) — do not use it.  (Pure
  one-site TDVP at fixed χ, without the expansion, is fine and is what the reference protocol uses.)
- **Bath modes**: the discrete bath correlation recurs at `t ≈ 2πN/w_max`, so matching the continuum
  over `t ≤ 12` needs `N ≳ 80` (N=16 → t ≈ 2.5).  `nmax ≥ 8` (or a per-mode schedule keeping
  ≥ 99 % of `⟨Q²⟩`) and χ ≳ 64 are needed, and crucially **`dt ≲ 0.1`** because the star's highest
  mode is `ω_max = wmax_factor*C` (`ω dt` must stay ≲ 3).
- Coding style: see the maintainer's shared `repo-rules` skill.

## Workflows

```bash
cd source
python reference.py                             # reference corrmat -> output/data/<prefix>/
python dmrg.py A-base_line                 # four cases -> output/data/A-base_line/
python FIG0@corrmat.py -ref ../../mkct4oqs/source/output/data
python FIG4@robustness.py -ref ../../mkct4oqs/source/output/data
python library/dmrg/validation.py          # regression suite
```

## CI

- `.github/workflows/ci.yaml` — syntax check, regression suite, and HEOM import.
- `.github/workflows/dmrg.yaml` — manual matrix over `<config> × <case>`; artifacts
  `dmrg-<config>-<case>`.  Retrieve with `gh run download <id> -D <tmp>`.

No token is needed: this repository is public and self-contained.

## Git style

See the maintainer's shared `repo-rules` skill (`type(scope): summary`; commit only when asked;
never commit secrets or `output/`).
