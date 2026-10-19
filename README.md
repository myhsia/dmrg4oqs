[![CI](https://github.com/myhsia/heomXdmrg-spin_boson/actions/workflows/ci.yaml/badge.svg?branch=main)](https://github.com/myhsia/heomXdmrg-spin_boson/actions)
[![License](https://img.shields.io/github/license/myhsia/heomXdmrg-spin_boson)](LICENSE)

Reference dynamics for the spin-boson system: HEOM and MPS/tDMRG
===============================================================

Two independent numerical references for the spin-boson model, used to
cross-check a companion memory-kernel project:

- **HEOM** — hierarchical equations of motion for linear coupling (`qutip`'s
  `OhmicEnvironment` + `HEOMSolver`) and an extended hierarchy for quadratic
  coupling (`extended_heom.py`).
- **MPS / tDMRG** — the joint system-plus-bath state evolved as a matrix
  product state (two-site TDVP) and reduced to populations, coherences, and
  auto-/cross-correlation functions.

This repository is the single source of truth for the **non-sensitive shared
pieces** as well: the quantum primitives (`qutils.py`), the config/data IO
helpers (`fsutils.py`), the plotting grid (`pltools.py`), and the parameter
files under `source/config/`.  A private companion repository (the
memory-kernel code) links to these rather than duplicating them.

Bath and model
--------------

The Ohmic spectral density is

    J(w) = pi * lambd / wC * w * exp(-w / wC)

discretized with Gauss-Legendre modes.  For tDMRG the quantum modes are rotated
to a one-dimensional chain by a Lanczos tridiagonalization and the
finite-temperature bath is a thermofield-double purification.  The four bath
cases are `QHM` (quantum harmonic), `QAH` (quantum anharmonic), `CHM` (classical
harmonic) and `CAH` (classical anharmonic); the classical cases use a per-mode
squeezed thermal state (a classical bath is not a quantum Gibbs state).

Layout
------

All code lives under `source/`, the working directory for every driver (paths
such as `./config` and `./output` are relative to it):

    source/
      heom.py                    # HEOM / extended-HEOM reference dynamics
      dmrg.py                    # MPS/tDMRG reference dynamics
      plot_dmrg.py               # dmrg@<prefix>.pdf figures
      environment.yaml           # conda environment (`dmrg`)
      requirements.txt
      config/                    # parameter files (shared, non-sensitive)
      library/
        extended_heom.py         # extended HEOM for quadratic coupling
        qutils.py fsutils.py pltools.py
        dmrg/                    # MPS/tDMRG implementation + regression suite

Install
-------

    conda env create -f source/environment.yaml
    conda activate dmrg

(`qutip` has no wheel on the pinned interpreter yet, so it is built from source
on the first environment creation.)

Run
---

From `source/` with the environment active:

```bash
python heom.py                             # HEOM reference -> output/data/<prefix>/
python dmrg.py A-base_line                 # four cases -> output/data/A-base_line/
python plot_dmrg.py A-base_line            # -> output/figures/dmrg@A-base_line.pdf
python library/dmrg/validation.py          # regression suite (ED references)
```

`dmrg.py` flags: `-case qhm|qah|chm|cah|all`, `-geometry auto|chain|star`,
`-N`, `-nmax`, `-chi`, `-dt`, `-dtau`, `-wmax_factor`, `-tag`.  Use the default
two-site engine (the single-site variant is numerically wrong).

Continuous integration
----------------------

- `.github/workflows/ci.yaml` — byte-compile, the tDMRG regression suite, and a
  HEOM import check.
- `.github/workflows/dmrg.yaml` — a manual (`workflow_dispatch`) matrix over
  `<config> × <case>`; each job uploads a `dmrg-<config>-<case>` artifact.

License
-------

MIT; see [`LICENSE`](LICENSE).
