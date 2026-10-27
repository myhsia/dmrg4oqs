# Guard Rules

The following directories, files, and patterns MUST NOT be read or accessed:

- `.git/`
- `*.log`
- `**/__pycache__/`
- `**/output/`

Any tool that reads files (`read`, `glob`, `grep`, `shell` with file-reading commands) must respect
these restrictions.

Scientific background and methodology are described in
[`.opencode/skills/SKILL.md`](./.opencode/skills/SKILL.md); load the `heom-dmrg` skill for the physics.

## Layout

All code is under `source/`, the working directory for the drivers (paths such as `./config` and
`./output` are relative to it):

- `source/dmrg.py` — driver: for one config it runs the four bath cases and writes the basis
  correlation matrix `output/data/<prefix>/corrmat_dmrg_<case>[-<tag>].npz` via `save_group`
  (key = case, `branch = component = ("00","01","10","11")`, `C(0) = I`).
- `source/plot_dmrg.py` — writes `output/figures/dmrg@<prefix>.pdf` (legend bar, combined body,
  coherence body), contracting each `corrmat_dmrg_<case>.npz` and the reference
  `corrmat_heom.npz` with `library/corrmat.py`.
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

Score with the house metric (`dmrg_error.py`, the shared `library/mre.py`): the global relative
Frobenius error `||C - C'||_F / ||C'||_F` of the basis correlation matrix (all 16 complex
elements) against the HEOM reference matrix over the config window, reported in percent.
The tabulated numbers below predate this metric (they were measured with the previous
32-series definition), so re-run `dmrg_error.py` to refresh them.

**Result (reference scenario, N=16, `nmax=8`, star, t ≤ 2):**

| setting | mean |
|---|---|
| χ=64, `dt=0.05` | **6.5 %** |
| χ=64, `dt=0.075` | 7.7 % |
| χ=128, `dt=0.05` | 10.9 % |
| χ=128, `dt=0.075`, `nmax=12` | 13.3 % |
| χ=128, `dt=0.075`, `nmax=8` | 14.5 % |
| N=40, χ=64, `dt=0.075` (t ≤ 4) | 44.9 % |
| χ=96, `nmax=16`, `dt=0.3` (first attempts) | 76.5 % |

For scale the FSSH reference scores 28 % on the same metric.  `dt` is converged at both bond
dimensions (7.7 → 6.5 for `dt` 0.075 → 0.05 at χ=64; 14.5 → 10.9 at χ=128) and the Fock space is
converged (`nmax` 8 → 12 moves 14.5 → 13.3), but **χ is not**: at fixed `dt=0.05` the score goes
6.5 (χ=64) → 10.9 (χ=128), i.e. **more χ is worse**, so neither value is certified.  Adding modes at
fixed χ is harmful (N=40 → 44.9 %), i.e. χ must grow with N.

**The decisive parameter is `dt`, not χ.**  The star's highest mode is `ω_max ≈ 40`, so `dt=0.3`
means `ω dt ≈ 12`; the per-step truncation then damps the fast modes and fabricates a "too fast
initial decay" that looks like a χ/`nmax`/thermal-state problem.  Measured at N=16, χ=64,
`nmax=8`:

| `dt` | pop₀ at t=0.4 / 0.8 / 1.2 (HEOM 0.866 / 0.631 / 0.481) |
|---|---|
| 0.3 | 0.718 / 0.447 / 0.391 (wrong) |
| 0.15 | 0.824 / 0.629 / 0.474 |
| 0.075 | 0.874 / 0.598 / 0.467 |
| **0.05** | **0.861 / 0.620 / 0.467** |

**Geometry.**  Prefer the **star** with a **per-mode `nmax`** (`-nmax` accepts a sequence): the star's
thermal state is an exact product of per-mode TFDs (no imaginary-time preparation, which dominated
the chain runs at 2 h+), and a schedule such as `[20, 8, 4, 2, …]` keeps 99.4 % of the static weight
`⟨Q²⟩` at 7.5× lower cost (`Σd³` 10 466 vs 78 608 for uniform 16).  The chain remains available but
its preparation is expensive and no better.

**Window.**  With N modes the discrete bath correlation recurs at `t ≈ 2πN/w_max`; for N=16,
`w_max=40` that is ≈2.5, so the reference is meaningful only up to `t ≈ 2` (the metric window is
still the config's `tf`; `dmrg_error.py` warns when the series ends earlier).

Pushing the recurrence out is **not** the way to extend it.  Measured, full window 0→12, same bath
discretization and χ=64 (mean error, house metric):

| run | qhm | qah | chm | cah |
|---|---|---|---|---|
| N=16, `dt=0.05`: 0–2.5 / 0–6 / 0–12 | 9.4 / 35.0 / 66.2 | 18.7 / 24.5 / 39.2 | 38.4 / 72.0 / 122.3 | 67.8 / 71.8 / 77.4 |
| N=40, `dt=0.1`: 0–2.5 / 0–6 / 0–12 | 18.7 / 65.2 / 111.1 | 67.0 / 86.0 / 112.1 | 97.3 / 103.6 / 118.3 | 97.6 / 114.4 / 141.2 |

N=40 is worse in *every* window and case — it collapses at `t≈2.2`, before its own recurrence at
6.3 — because the extra modes add entanglement that χ=64 truncates harder.  So the late-time failure
is **truncation, not the recurrence**: a converged 0–12 curve needs N≳80 *and* χ≳256 together, which
is out of reach (N=40/χ=128/`dt=0.05` alone was ≈40 h).  The reference stays a short-time statement.


**Classical cases.**  `chm`/`cah` collapse by `t ≈ 0.6` at every setting: the per-mode squeezed
thermal state is not a Gibbs state, so they are not a quantitative reference and are **not plotted**.

**Regression.**  `library/dmrg/validation.py` reproduces exact diagonalization for N ≤ 3 to ~1e-2,
so the mapping and the observable algebra are sound; `check_bath.py` measures the bath itself.


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
python plot_dmrg.py A-base_line            # -> output/figures/dmrg@A-base_line.pdf
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
