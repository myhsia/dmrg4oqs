# Guard Rules

The following directories, files, and patterns MUST NOT be read or accessed:

- `.git/`
- `**/__pycache__/`

Any tool that reads files (`read`, `glob`, `grep`, `shell` with file-reading commands) must respect
these restrictions.

Scientific background and methodology are described in
[`.opencode/skills/SKILL.md`](./.opencode/skills/SKILL.md); load the `heom-dmrg` skill for the physics.

## Shared rules

Repository etiquette (preserve uncommitted changes and staging, commit by pathspec, keep the guides
in sync), Python coding style, the git/PR workflow, LaTeX writing and source conventions,
tooling/upkeep, and local/CI parity live in the global `repo-rules` skill
(`~/.config/opencode/skills/repo-rules/SKILL.md`).

## Layout

All code is under `source/`, the working directory for the drivers (paths such as `./config` and
`./output` are relative to it):

- `source/dmrg.py` — driver: for one config it runs the bath cases and writes the basis
  correlation matrix per case (`corrmat_dmrg_<case>[-<tag>].npz`) and, after a complete
  `-case all` run, the grouped `corrmat_dmrg.npz` (one matrix per case,
  `branch = component = ("00","01","10","11")`, `C(0) = I`).
- `source/FIG0@corrmat.py`, `source/FIG4@robustness.py` — synced copies of the canonical
  figure scripts: `python FIG0@corrmat.py` draws the 4x4 `C_ij(t)` pages of every method found
  here and `python FIG4@robustness.py` the `dmrg@<prefix>.pdf` comparison figure.  Both read the
  reference `corrmat_heom.npz` from the data root `./output/data` by default; this repository
  generates its own with `reference.py`, so no link is needed (`-ref DIR` points at another root).
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
  observable extraction, the star path), `krylov.py` (global Krylov engine: exact on the projected
  Krylov space, worse than TDVP at fixed `chi`), `integrator.py` (propagator comparison on the
  exact state: exact exponential vs Krylov vs Cayley), `validation.py` (regression suite).
- `source/environment.yaml`, `source/requirements.txt` — the `dmrg` conda environment.

**Canonical archives** are `corrmat_<method>[_<label>].npz`.  The per-curve
`{pop,cohe,cohem,corr,xorr}_<method>_<case>.npz` files are legacy from the pre-corrmat workflow;
no script reads them, and they are archived out of `output/data/<prefix>/`.

## Environment

- conda environment `dmrg` (`source/environment.yaml` + `source/requirements.txt`), interpreter
  pinned to `python=3.15`; the CI action `dmrg-env` builds the same environment, so local and CI
  agree.
- Run every command from `source/` with the environment active; the drivers resolve `./config`,
  `./output` and `../../mkct4oqs/...` relative to it.
- Entry points: `dmrg.py` (dynamics), `reference.py` (HEOM / extended DEOM), `dmrg_error.py`
  (score), `FIG0@corrmat.py` / `FIG4@robustness.py` (figures), `check_bath.py` and
  `library/dmrg/validation.py` (validation).  Data land in `output/data/<prefix>/`, figures in
  `output/figures/`; both are gitignored.

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

**Only QHM has a reference.**  The four cases differ in the *bath* (quantum/classical x
harmonic/anharmonic) while `corrmat_heom.npz` holds the quantum-harmonic matrix only.  Scoring
QAH/CHM/CAH against it measures the physical difference between the baths, not the numerical
error, so **no accuracy is claimed for them and their cross-case numbers are not to be
quoted**.  The numerical protocol (N, `wmax_factor`, chi, `nmax`, `dt`) is calibrated on QHM
below and inherited unchanged by the other three cases; those runs feed the quantum-vs-
classical comparison, not an error analysis.

### Accuracy (QHM only)

Star, N=16, `wmax_factor=10`, chi=64, `dt=0.025`, scored on the production window **`t <= 4`**:

| `nmax` | E at `t <= 4` | E at `t <= 2` |
|---|---:|---:|
| 8 | 7.58 % | 3.75 % |
| 12 | 6.94 % | 2.55 % |
| **16** | **6.62 %** | **1.77 %** |

The `t <= 4` differences halve (0.64 then 0.32), so the limit is ~6.45 %; the `t <= 2` metric is
much more sensitive to the Fock space (3.75 -> 1.77, a factor 2), i.e. the early-time coherence
is the part that needs the high levels.  `nmax = 16` is the best measured; `nmax = 12` costs
+5 % error for a much cheaper run.

Window dependence of the same protocol (`nmax=8`, `dt=0.025`): 3.75 / 5.31 / 7.58 / 10.53 /
12.52 % at `t <= 2 / 3 / 4 / 5 / 6`.  The reporting limit is `t = 4`; `t = 6` is available and
labelled.

### How each parameter was shown to converge (QHM, `t <= 4`)

Every entry below changes **one** parameter of the reference protocol and is scored on the same
window, so the ordering is the measured sensitivity:

| parameter | setting | E at `t <= 4` |
|---|---|---:|
| `dt` (`nmax=8`) | 0.05 | 26.70 % |
| | **0.025** | **7.58 %** |
| `dt` (`nmax=12`) | 0.025 | 6.94 % |
| | 0.0125 | 6.93 % |
| chi | 32 | 7.60 % |
| | **64** | 7.58 % |
| `nmax` | 8 | 7.58 % |
| | 12 | 6.94 % |
| | **16** | **6.62 %** |
| N (`wmax_factor=10`) | **16** | **7.58 %** |
| | 40 | 8.06 % |
| N=40, chi 32 -> 64 | | 3.85 % -> 3.85 % at `t <= 1` |

Order to use when certifying a new window: **window -> `dt` -> chi -> `nmax` -> N**.

The `dt=0.05` and `chi=32` scan runs are **not retained**: their numbers in the table are the
record and are reproducible from the settings given here.  The N=32/40 runs *are* kept, because
they are the only evidence for the non-monotonicity in N; the other archives under
`output/data/<prefix>/` are production runs.

* **`dt` dominates, and 0.025 is converged from below.** 0.05 is 3.5x worse than 0.025 at
  `t <= 4`, while halving again (0.025 -> 0.0125, at `nmax=12`) moves the score by 0.01 points
  (6.94 -> 6.93 %).  The improvement is per-step bond truncation (see the propagator section),
  not the exponential, which is exact at 0.05 already.
* **chi is cheap to satisfy**: 32 and 64 agree to 0.02 points, so 32 is enough at N=16.
* **`nmax = 8` is badly unconverged**: 12 buys 8 % at `t <= 4` and 32 % at `t <= 2`, and 16 buys
  a further 5 % / 31 %.  This was invisible while `dt = 0.05` truncation dominated, and it is the
  one parameter still moving at the affordable end of the sweep.
* **agreement with HEOM is not monotone in N**: N=40 is worse than N=16 at *every* window,
  including `t <= 1` where aliasing is irrelevant, and chi=32 and chi=64 give the *same* 3.85 %
  there -- so it is the discretisation, not truncation.  With the static weight identical for
  `wmax_factor=10` at every N (`sum c^2 eta` = 2.15982), the likely reading is that HEOM's own
  bath is a fit (ESPIRA-I) and a finer Gauss-Legendre representation of the true continuum can
  move away from it.  Keep N=16.

### Method comparison, QHM only, `t <= 4`

Same metric, same window (MKCT archives in `mkct4oqs`, MASH/FSSH in `sh4oqs`):

| method | E |
|---|---:|
| MKCT | **0.54 %** |
| tDMRG (this reference) | 7.58 % (6.94 % at `nmax=12`) |
| MASH | 17.16 % |
| FSSH | 37.21 % |

Note the metric window: runs are made with `-t1 4` while the config `tf` is 12, so
`dmrg_error.py` scores the overlap `[0, 4]` and warns; use `--windows 2 3 4` to score a long
run on successive windows.

**The propagator is not the limitation (measured; `library/dmrg/integrator.py`).**  On the exact
many-body state, with every stepped run ending at the same final time so that a perfect
propagator must show no `dt` dependence:

| dt | w_max dt | exact | Krylov (m=20) | Cayley |
|---:|---:|---:|---:|---:|
| 0.2 | 8.0 | 1e-14 | 4.5e-4 | 1.5e-1 |
| 0.1 | 4.0 | 1e-14 | 1.4e-3 | 1.6e-3 |
| **0.05** | **2.0** | **1.1e-14** | **2.1e-5** | 1.3e0 |
| 0.025 | 1.0 | 1.1e-14 | 8.7e-15 | 8.3e-1 |

(free mode at `w = 40`, `t = 2`, error in `<q>`; the thermal 2-mode bath gives the same
picture, 1e-14 Krylov vs 1e-2 Cayley at `dt = 0.05`).  So at `dt = 0.05` the propagator is
exact to round-off, yet that setting is 3.5x worse at `t <= 4` (26.70 % vs 7.58 %).  The `dt`
dependence is therefore **bond truncation** per step, not integrator accuracy: `U(dt)^n =
U(n dt)` for any propagator built from `H`, so a converged exponential cannot show a `dt` error
at fixed output times.

A **Cayley** (Crank-Nicolson) step is therefore not an alternative: it is exactly unitary but
composes to `exp(-i t H_eff)` with `H_eff = 2 atan(H dt/2)/dt`
(`validation.py::test_propagators`, 2e-13), which is *nonlinear* in the level index, i.e. it
anharmonicises a harmonic bath instead of just shifting it.  Its `n = 1` shift at the fastest
mode is 21 % at `dt = 0.05`, and matching Krylov needs ~6x smaller steps, each an iterative
linear solve.  Cayley is the right scheme for **imaginary time** (`prepare_thermal`), where the
propagator is not unitary and A-stability is what matters.

**Geometry.**  Prefer the **star**: its thermal state is an exact product of per-mode TFDs, with
no imaginary-time preparation (which dominated the chain runs at 2 h+).  The chain remains
available but its preparation is expensive and no better.

**Per-mode `nmax`: measured, and the weight-based schedule is *wrong*.**  `-nmax` accepts a
per-mode sequence (e.g. `-nmax 20,8,4,2`), and a schedule assigning a large Fock space to the
modes with the largest static weight looks attractive -- `[2,2,8,20,4,2,...]` costs
`sum d^3 = 10 466` against 11 664 for uniform `nmax = 8`.  Measured at the converged step
(qhm, star, N=16, chi=64, `dt = 0.025`, `t <= 2`; the schedule run is t<=2 only, scored on the
same window as the uniform-8 benchmark):

| `nmax` | E |
|---|---:|
| uniform 8 | **3.75 %** |
| `[2,2,8,20,4,2,...]` (weight-ordered) | 12.35 % |

The schedule gives `nmax = 2` to the two lowest modes (`w = 0.21, 1.11`), which are the ones
with the largest thermal occupation (`<n> = 7.4, 1.1`): their *coupling* weight `c^2 eta` is
small, but they carry the long-time memory, so two levels corrupt the dynamics.  Static weight
is therefore not a valid criterion for the Fock truncation; use **uniform `nmax`**.  Note that
uniform 8 is itself not converged: at `dt = 0.025`, `nmax = 12` buys 8 % at `t <= 4` and 32 % at
`t <= 2` (see the parameter table above), and 16 improves it again.  Use **uniform
`nmax = 16`** for the accuracy statement and 12 for cheaper production runs.

The earlier claim that a schedule keeps "99.4 % of the static weight" conflated two different
truncations: 99.44 % is the weight of keeping 10 of the 16 *modes*, while 10 466/78 608 is the
*cost* ratio of the Fock schedule -- the schedule itself keeps 88.9 % (star), not 99.4 %.

**Window.**  The discrete bath correlation recurs at `t ~ 2 pi N / w_max` (N=16 -> 2.5), but
that bound is **soft**: the measured series above (3.75 -> 12.52 % from `t <= 2` to `t <= 6`)
grows smoothly with no breakdown, so `t = 4` (1.6x the recurrence) is a legitimate window.  Do
not push it by adding modes: at fixed chi, more modes add entanglement that is truncated
harder, and the measured N=32/N=40 runs are worse at every window.  A long-window run should
start from `chi = 64`, `dt = 0.025`, N=16 and raise N only with chi raised in step, timing a
short segment first.

**Other cases.**  `qah`/`chm`/`cah` have no reference (see "Only QHM has a reference" above):
their baths differ physically, so their dynamics differ from HEOM by construction.  They run
with the QHM-calibrated protocol and are reported only in the quantum-vs-classical comparison.
The classical cases additionally use a per-mode squeezed thermal state, which is not a Gibbs
state.

**Regression.**  `library/dmrg/validation.py` reproduces exact diagonalization for N <= 3 to
~1e-2, so the mapping and the observable algebra are sound; `check_bath.py` measures the bath
itself against `sum_j c_j^2 (1/2) coth(beta w_j/2) cos(w_j t)`.


## Shared pieces (mirrored from `mkct4oqs`)

The **canonical** `source/reference.py`, `source/config/`, and
`source/library/{qutils,fsutils,pltools,extended_deom,mre,corrmat}.py` live in the private `mkct4oqs`
repository.  This
repository keeps **synced copies** so that a public clone and the CI work without a token.
The latest sync adds `GridPlt._pin_xlabel` to the synced `pltools.py` (a baseline-pinned
x-label, used by the canonical Padé figure); it is unused in this repository.

> **Sync reminder.**  After changing any of them in `mkct4oqs`, run
> `bash tools/sync_from_mkct.sh` here, then commit and push.  `sh4oqs` mirrors the same pieces
> and needs the same sync.  The synced `tools/link_heom.sh` is the no-op case here: this
> repository generates the reference itself with `reference.py`, so it holds the real
> `corrmat_heom.npz` archives and `FIG0`/`FIG4` read them from `./output/data`.

`fsutils.py` also carries the grouped-archive writer/loader (`save_group`/`load_group`) used by
the reference (`corrmat_heom.npz`) and the companion's MKCT results: one entry per case under its
case name, each with its **own** time axis stored as `t_<case>` (there is no shared `t`), plus
optional `identifier`/`component`/`branch` metadata arrays; `corrmat.py` contracts the matrix into
the five observables.

## Conventions

- Set `sys.pycache_prefix = os.path.abspath("./output/pycache")` **after** third-party imports and
  **before** the first project import in every run script.
- Run everything from `source/` with the `dmrg` environment active.
- Four cases: `QHM`, `QAH`, `CHM`, `CAH`.  The anharmonic quartic coefficient is
  `a = anharmon * sqrt(7/6)`.
- The **classical** cases use a per-mode squeezed thermal state (the classical bath lies below the
  quantum zero-point and is not a Gibbs state).
- **`engine`**: two-site TDVP is the reference engine.  The single-site + subspace-expansion
  variant (`-engine onesite`) is numerically wrong (it oscillates/collapses) — do not use it.
  (Pure one-site TDVP at fixed χ, without the expansion, is fine and is what the reference
  protocol uses.)  `-engine krylov` builds the global Krylov space
  `{psi, H psi, ..., H^{m-1} psi}` as MPS (`m = 10`), exponentiates the `m x m` projection
  exactly and truncates once when forming the combination.  It is the vector-space analogue
  of the exact moment expansion, it agrees with ED to 1e-9 (one mode) and with
  TDVP to 6e-6 on the same model (`validation.py`) -- and it is **much worse than TDVP at fixed
  chi, and does not improve with chi**.  Measured (qhm, star, N=16, nmax=8, `dt = 0.025`,
  `t <= 2`):

  | engine | chi = 16 | chi = 32 | chi = 64 |
  |---|---:|---:|---:|
  | two-site TDVP | 8.26 % | 3.76 % | **3.75 %** |
  | global Krylov (`m = 8`) | 67.29 % | 72.44 % | -- |

  Two facts matter.  TDVP converges cleanly in chi (8.26 -> 3.76 -> 3.75), so chi = 32 is
  enough.  Krylov is 8-19x worse and gets *worse* as chi grows -- which contradicts the naive
  "truncation multiplicity" explanation (that picture predicted convergence toward TDVP at
  large chi).  What is left is the repeated compression itself: the Krylov vectors are
  re-compressed at every matvec and addition, the Gram deviates from the identity
  (`~1e-4` at chi = 16), and the generalised eigenvalue solve then produces a poor
  combination; at chi = 16 the population collapses to 0.095 at `t = 2` (HEOM 0.415) with the
  norm still 1, i.e. artificial relaxation.

  So the conclusion is about **this construction**, not about the Krylov idea: a fair test
  needs a work bond dimension (`chi_work >> chi_out`, truncate once at the end) or a
  variational compression of the Krylov vectors, which is `~m` times the cost per step.  Until
  that exists, use TDVP; the Krylov engine stays as a validated reference implementation of
  the moment/CF analogy and as the regression that the global step is exact when nothing is
  truncated.
- **Bath modes**: the discrete bath correlation recurs at `t ≈ 2πN/w_max`, so matching the continuum
  over `t ≤ 12` needs `N ≳ 80` (N=16 → t ≈ 2.5).  **uniform `nmax = 16`** (8 and 12 are not
  converged; the weight-based per-mode schedule is not valid, see above) and χ ≥ 32 are needed.  **`dt ≲ 0.025`** at χ = 64 (`0.05 -> 9.0 %`,
  `0.025 -> 3.75 %`, `0.0125 -> 3.77 %`); the limit is per-step bond truncation, *not* the
  exponential, which is exact to round-off at `dt = 0.05` (see the propagator section).
- **`A-base_line_tf6`** is a **testing-only** config, absent from every `RUNLIST`: it inherits
  `A-base_line` and cuts `xlims` to `[[0, 6], 2]`.  It serves the window/convergence checks, not
  the demo figures, because the DMRG reference is not stable over the longer window (see
  Validity).  It is synced from `mkct4oqs` with the rest of `source/config/`.
- Coding style: see the maintainer's shared `repo-rules` skill.

## Workflows

```bash
cd source
python reference.py -prefix A-base_line          # HEOM corrmat -> output/data/A-base_line/
python dmrg.py A-base_line -geometry star -N 16 -nmax 16 -chi 64 -dt 0.025 \
       -t0 0 -t1 4 -nproc 4 -njobs 2             # four cases -> corrmat_dmrg_*
python dmrg_error.py --cases qhm --window 4       # Frobenius table -> output/dmrg_error.{json,md}
python FIG0@corrmat.py                            # corrmat-dmrg@<prefix>.pdf
python FIG4@robustness.py                         # dmrg@<prefix>.pdf
python check_bath.py -N 16 -nmax 16 -chi 96       # bath correlator against the analytic curve
python library/dmrg/validation.py                 # regression suite
python library/dmrg/integrator.py                 # propagator study -> output/integrator_report.{json,md}
```

The flags above are the calibrated protocol; `python dmrg.py A-base_line` with no flags runs the
defaults (`N=20`, `nmax=6`, `chi=128`, `dt=0.02`), which is **not** the protocol the validity table
was measured with.

## Validation

Two checks, both runnable locally without the CI:

```bash
cd source
python library/dmrg/validation.py     # exact diagonalization for N <= 3, tolerance ~1e-2
python check_bath.py                  # static <Q^2> fraction and the bath correlator
```

- `library/dmrg/validation.py` covers the bare two-level system, one harmonic mode, one quartic
  mode, a chain multimode case, the Krylov engine and the propagators.  It reproduces exact
  diagonalization for N <= 3 to `~1e-2`, which is the evidence that the mapping and the observable
  algebra are right.
- `check_bath.py` evolves the bath alone and prints the static `<Q^2>` fraction and the relative
  error of `Re<Q(t)Q(0)>` against `sum_j c_j^2 (1/2) coth(beta w_j/2) cos(w_j t)`; a thermal state
  that is too mixed shows up here first.
- The scientific statement (which parameter was shown to converge, and to what) is the "Validity
  (measured, not assumed)" section above; this section only covers the commands and criteria.

## CI
**Dependencies: the bot proposes, CI decides** (the policy and its reasoning are in the
machine-global `repo-rules` skill, "Local and CI parity").  `.github/dependabot.yml` watches
`github-actions` daily and, for the `conda` ecosystem, ignores **patch** releases of `python` -- no
functional change, and the class most likely to break the compiled extensions.  The pin in
`source/environment.yaml` moves deliberately, in step with the sibling repositories, and only once
`build_env` resolves the environment and the suite is green.  `requirements.txt` is intentionally
unpinned, so the libraries float to latest at install time.

- `.github/workflows/ci.yaml` — syntax check, regression suite, and HEOM import; runs on push and
  on the daily schedule (`cron: '41 9 * * *'`).  The environment cache key carries the UTC date, so
  the scheduled run re-solves the dependency set once a day.
- `.github/workflows/dmrg.yaml` — manual matrix over `<config> × <case>`; artifacts
  `dmrg-<config>-<case>`.  Retrieve with `gh run download <id> -D <tmp>`.
- `.github/workflows/bath.yaml` — manual bath-validity run of `check_bath.py` at the given
  `N` / `nmax` / `chi` / `tf` / `dt`; kept separate from `ci.yaml` so a push cannot cancel a long
  validity run.
- `.github/actions/dmrg-env/action.yml` — composite action shared by the three workflows: caches
  the conda package dir and pip (keyed on `source/environment.yaml` / `source/requirements.txt`),
  restores the `dmrg` environment, and builds/saves it when the `build` input is true.

No token is needed: this repository is public and self-contained.

## Git style

See the maintainer's shared `repo-rules` skill (`type(scope): summary`; commit only when asked;
never commit secrets or `output/`).

