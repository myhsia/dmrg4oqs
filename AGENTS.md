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

- `source/dmrg.py` — driver: for one config it runs the four bath cases and writes
  `output/data/<prefix>/{pop,cohe,cohem,corr,xorr}_dmrg_<case>[-<tag>].npz` via `fsutils.save_time`
  (`t`, `y`, `identifier`).
- `source/plot_dmrg.py` — writes `output/figures/dmrg@<prefix>.pdf` (legend bar, combined body,
  coherence body); overlays `*_heom.npz` when present.
- `source/heom.py` — HEOM reference dynamics (qutip `HEOMSolver` for linear coupling; the extended
  HEOM/DEOM in `library/extended_heom.py` for quadratic coupling), writing
  `output/data/<prefix>/{pop,cohe,cohem,corr,xorr}_heom.npz`.
- `source/library/extended_heom.py` — extended HEOM for polynomial coupling.
- `source/library/dmrg/` — `chain.py` (star→Lanczos chain, two-site TDVP, thermofield-double bath),
  `model.py` (star MPS, coupling MPO, per-mode squeezed classical state), `tdmrg.py` (shared
  observable extraction), `validation.py` (regression suite).
- `source/environment.yaml`, `source/requirements.txt` — the `dmrg` conda environment.

## Method

- **HEOM.**  Linear coupling uses qutip's `OhmicEnvironment` + `HEOMSolver` (ESPIRA-I fit); quadratic
  coupling uses the extended HEOM/DEOM (`extended_heom.py`).  `heom.py` runs the four cases for each
  config in its `RUNLIST` and writes the five `*_heom.npz` files.
- **MPS / tDMRG.**  The star-shaped Ohmic bath (shared Gauss-Legendre discretization) is rotated to a
  one-dimensional chain by a Lanczos tridiagonalization (coupling vector as seed), so the system
  couples only to the first chain site; the finite-temperature bath is a thermofield-double
  purification prepared by imaginary-time TDVP.  All five observables follow from two propagations
  `|phi_b(t)> = U|b>|TFD>` (b = 0, 1), because every time-zero system operator is rank one.

## Validity (measured, not assumed)

Score with the house metric (`dmrg_error.py`, the shared `library/mre.py`): the relative `L2` error
of `{Re P0, Re P1, Re C_xx, Re C_xy, Im C_xx, Im C_xy}` against HEOM over the config window,
averaged and reported in percent.

**Result (reference scenario, t ≤ 2):** star geometry, N=16, `nmax=8`, χ=64, **`dt=0.075`** gives
**7.7 %** (per curve 5.0–10.8 %) — about four times better than the FSSH reference (28 %).  So the
tDMRG reference *is* quantitative in this window.

**The decisive parameter is `dt`, not χ.**  The star's highest mode is `ω_max ≈ 40`, so `dt=0.3`
means `ω dt ≈ 12`; the per-step truncation then damps the fast modes and fabricates a "too fast
initial decay" that looks like a χ/`nmax`/thermal-state problem.  Measured at N=16, χ=64,
`nmax=8`:

| `dt` | pop₀ at t=0.4 / 0.8 / 1.2 (HEOM 0.866 / 0.631 / 0.481) |
|---|---|
| 0.3 | 0.718 / 0.447 / 0.391 (wrong) |
| 0.15 | 0.824 / 0.629 / 0.474 |
| **0.075** | **0.874 / 0.598 / 0.467** |

**Geometry.**  Prefer the **star** with a **per-mode `nmax`** (`-nmax` accepts a sequence): the star's
thermal state is an exact product of per-mode TFDs (no imaginary-time preparation, which dominated
the chain runs at 2 h+), and a schedule such as `[20, 8, 4, 2, …]` keeps 99.4 % of the static weight
`⟨Q²⟩` at 7.5× lower cost (`Σd³` 10 466 vs 78 608 for uniform 16).  The chain remains available but
its preparation is expensive and no better.

**Window.**  With N modes the discrete bath correlation recurs at `t ≈ 2πN/w_max`; for N=16,
`w_max=40` that is ≈2.5, so the reference is meaningful only up to `t ≈ 2` (the metric window is
still the config's `tf`; `dmrg_error.py` warns when the series ends earlier).  A longer window needs
larger N (N=40 → t ≈ 6).

**Classical cases.**  `chm`/`cah` collapse by `t ≈ 0.6` at every setting: the per-mode squeezed
thermal state is not a Gibbs state, so they are not a quantitative reference and are **not plotted**.

**Regression.**  `library/dmrg/validation.py` reproduces exact diagonalization for N ≤ 3 to ~1e-2,
so the mapping and the observable algebra are sound; `check_bath.py` measures the bath itself.


## Shared pieces (mirrored from `mkct4oqs`)

The **canonical** `source/heom.py`, `source/config/`, and
`source/library/{qutils,fsutils,pltools,extended_heom}.py` live in the private `mkct4oqs`
repository.  This
repository keeps **synced copies** so that a public clone and the CI work without a token.

> **Sync reminder.**  After changing any of them in `mkct4oqs`, run
> `bash tools/sync_from_mkct.sh` here, then commit and push.  `sh4oqs` mirrors the same pieces
> and needs the same sync.

`fsutils.py` also carries the grouped-archive writer/loader (`save_group`/`load_group`) used by
the companion's MKCT results.

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
python heom.py                             # HEOM reference -> output/data/<prefix>/
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
