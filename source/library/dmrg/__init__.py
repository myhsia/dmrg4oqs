"""MPS / tDMRG reference dynamics for the spin-boson model.

A temporary, self-contained reference that time-evolves the spin-boson
wavefunction as a matrix product state (tDMRG, i.e. DMRG time evolution) and
extracts the usual observable set: populations, coherences, and the
auto-/cross-correlation functions.

The bath is a set of Harmonic oscillator modes in a *star* geometry
(the two-level system couples to every mode); the finite-temperature bath is
represented by a purification (thermofield-double) MPS.  See
:mod:`dmrg.model` for the Hamiltonian and initial state and :mod:`dmrg.tdmrg`
for the evolution and observables.
"""
from .tdmrg import run_observables  # noqa: F401

__all__ = ["run_observables"]
