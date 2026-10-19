# Guard Rules

The following directories, files, and patterns MUST NOT be read or accessed:

- `.git/`
- `*.log`
- `**/__pycache__/`
- `**/output/`

Any tool that reads files (`read`, `glob`, `grep`, `shell` with file-reading commands) must respect
these restrictions.

## Layout

All code is under `source/`, the working directory for the drivers (paths such as `./config` and
`./output` are relative to it):

- `source/dmrg.py` — driver: for one config it runs the four bath cases and writes
  `output/data/<prefix>/{pop,cohe,cohem,corr,xorr}_dmrg_<case>[-<tag>].npz` via `fsutils.save_time`
  (`t`, `y`, `identifier`).
- `source/plot_dmrg.py` — writes `output/figures/dmrg@<prefix>.pdf` (legend bar, combined body,
  coherence body); overlays `*_heom.npz` when present.
- `source/library/dmrg/` — `chain.py` (star→Lanczos chain, two-site TDVP, thermofield-double bath),
  `model.py` (star MPS, coupling MPO, per-mode squeezed classical state), `tdmrg.py` (shared
  observable extraction), `validation.py` (regression suite).
- `source/environment.yaml`, `source/requirements.txt` — the `dmrg` conda environment.

## Shared files (private companion repository)

`source/config/` and `source/library/{qutils,fsutils,pltools}.py` are **not** in this repository.
Locally they are symlinked by `tools/setup_symlinks.sh` (gitignored); CI checks them out from the
private companion repo with the secret `MKCT_REPO_TOKEN` and copies them into `source/`.

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
- Coding style: spaces around `=` in keyword arguments and decorator parameters, one parameter per
  line in multi-line signatures, `) -> type` on the closing line, and a `# %%` marker at the top of a
  module.

## Workflows

```bash
cd source
python dmrg.py A-base_line && python plot_dmrg.py A-base_line
python library/dmrg/validation.py
```

## CI

- `.github/workflows/ci.yaml` — syntax check + regression suite.
- `.github/workflows/dmrg.yaml` — manual matrix over `<config> × <case>`; artifacts
  `dmrg-<config>-<case>`.  Retrieve with `gh run download <id> -D <tmp>`.

Both need `secrets.MKCT_REPO_TOKEN` to fetch the shared files.

## Git style

Commits follow `type(scope): short summary` (`feat`, `fix`, `tune`, `refactor`, `docs`, `chore`),
imperative and concise.  Only commit when explicitly requested; never commit secrets or `output/`
artifacts.
