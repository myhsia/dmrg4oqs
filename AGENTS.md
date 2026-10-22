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

`check_bath.py` (run by the dispatch-only `bath` workflow) compares the chain bath with the exact
discrete modes.  At the reference settings (N=16, `nmax=16`, χ=96):

- static weight `⟨Q²⟩` = **95.2 %** of `Σ c²½coth(βω/2)` (so `nmax=16` is needed; `nmax=8` gives 85 %);
- `Re⟨Q(t)Q(0)⟩` relative error **0.126**, but the curve is **over-damped**: 0.067 vs 0.34 at `t=0.3`
  and ≈0 afterwards, i.e. the real-time *bath* dynamics is faithful only for `t ≲ 0.3`.

Consequences for the reference:

- qhm tracks HEOM within ≈0.08 up to `t≈0.3` and ≈0.12 up to `t≈1.1`, then departs (0.234 over 0→2);
  qah collapses (0.609).  Raising χ (32→96) or `nmax` (8→16) does **not** help: with the correct
  weight the over-damped bath is exposed.  The limit is the bath's real-time MPS evolution, not the
  system-side algebra (which reproduces ED for N≤3 to ~1e-2 in `library/dmrg/validation.py`).
- The classical cases (`chm`/`cah`) collapse by `t≈0.6` at every setting: the per-mode squeezed
  thermal state is not a Gibbs state, so they are not a quantitative reference and are **not plotted**.
- The figure is therefore drawn with `python plot_dmrg.py A_base_line -cases qhm qah -tmax 2`
  (quantum cases only, restricted window), overlaying HEOM.

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
  (`-engine onesite`) is numerically wrong (it oscillates/collapses) — do not use it.
- **Bath modes**: the discrete bath correlation recurs at `t ≈ 2πN/w_max`, so matching the continuum
  over `t ≤ 12` needs `N ≳ 80`; `nmax ≥ 8`, `chi ≳ 64` and `dt ≲ 0.3` are needed too.
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
