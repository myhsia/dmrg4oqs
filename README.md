[![CI](https://github.com/myhsia/dmrg-spin_boson/actions/workflows/ci.yaml/badge.svg?branch=main)](https://github.com/myhsia/dmrg-spin_boson/actions)
[![License](https://img.shields.io/github/license/myhsia/dmrg-spin_boson)](LICENSE)

MPS / tDMRG reference dynamics for the spin-boson system
========================================================

A temporary, self-contained wavefunction reference for the spin-boson model:
the joint system-plus-bath state is time-evolved as a **matrix product state**
(two-site TDVP) and reduced to populations, coherences, and auto-/cross-
correlation functions.  It is used to cross-check the numerically exact
hierarchy and semiclassical references of a companion project; it is not itself
part of any publication.

The bath is set up in the **star** geometry from a Gauss-Legendre discretization
of the Ohmic spectral density

    J(w) = pi * lambd / wC * w * exp(-w / wC)

and, for the quantum cases, rotated to a **one-dimensional chain** by a Lanczos
tridiagonalization so that the system couples only to the first site.  The
finite-temperature bath is a thermofield-double purification prepared by
imaginary-time TDVP.

Layout
------

All code lives under `source/`, the working directory for every driver (paths
such as `./config` and `./output` are relative to it):

    source/
      dmrg.py                    # driver (one config, four bath cases)
      plot_dmrg.py               # dmrg@<prefix>.pdf figures
      environment.yaml           # conda environment (`dmrg`)
      requirements.txt
      library/dmrg/              # MPS/tDMRG implementation + regression suite

The conda environment is created from `source/environment.yaml`:

    conda env create -f source/environment.yaml
    conda activate dmrg

Shared primitives
-----------------

This repository keeps only the DMRG-specific code.  A few helpers (`qutils.py`,
`fsutils.py`, `pltools.py`) and the parameter files (`source/config/`) live in a
private companion repository and are brought in two ways:

- **Locally** they are symlinked by `tools/setup_symlinks.sh` (the links are
  gitignored, so the DMRG code edits them in place while developing);
- **in CI** they are checked out from the private repository with a fine-grained
  token stored as the secret `MKCT_REPO_TOKEN` (contents: read) and copied into
  `source/` before the run.  Without that secret the `validation` and `dmrg`
  workflows cannot fetch them.

Run
---

From `source/` with the `dmrg` environment active:

```bash
python dmrg.py A-base_line                 # four cases -> output/data/A-base_line/
python plot_dmrg.py A-base_line            # -> output/figures/dmrg@A-base_line.pdf
python library/dmrg/validation.py          # regression suite (ED references)
```

`dmrg.py` flags: `-case qhm|qah|chm|cah|all`, `-geometry auto|chain|star`,
`-N`, `-nmax`, `-chi`, `-dt`, `-dtau`, `-wmax_factor`, `-tag`.

The four bath cases are `QHM` (quantum harmonic), `QAH` (quantum anharmonic),
`CHM` (classical harmonic) and `CAH` (classical anharmonic).  The classical
cases use a per-mode squeezed thermal state matching the classical variance,
because a classical bath is not a quantum Gibbs state.

Continuous integration
----------------------

- `.github/workflows/ci.yaml` — byte-compile check and the regression suite.
- `.github/workflows/dmrg.yaml` — a manual (`workflow_dispatch`) matrix over
  `<config> × <case>`; each job uploads a `dmrg-<config>-<case>` artifact.

Because the shared primitives live in a private repository, both workflows
check it out with `secrets.MKCT_REPO_TOKEN`.

License
-------

MIT; see [`LICENSE`](LICENSE).
