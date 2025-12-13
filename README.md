# CaΣoS-Example package

This repository contains tutorials, examples and implementations from paper or textbooks. The purpose is to provide an easy on-boarding with the CaΣoS  toolbox, a nonlinear optimization-oriented sum-of-squares toolbox based on the symbolic framework of CasADi.

### Installation
The example package itself does not not need to be installed. Only CaΣoS and a supported conic solver is needed. Follow the instructions on the [Getting started](https://github.com/ifr-acso/casos/wiki#getting-started) page.

### Folder Structure

```text
Getting Started/	# Tutorials
├── Basics        	# Toolbox basics
├── SOS			  	# How to setup different kinds of SOS problems

Publications/ 					# Code from publications of the iFR-ACSO group
├── Benchmark					# ACC25: Benchmarks of CaΣoS with other state-of-the art SOS toolboxes.
├── Infinitesimal-MPC			# TAC25: Code to synthesize a pair of compatible CBF-CLF to approximate NMPC for Spacecraft Attitude Control.
├── Verify-Manifold-Spacecraft  # ECC26: Stability Verification on Manifolds with Applications in Spacecraft Attitude Control.
```



### CaΣoS Quick links

- [Getting started](https://github.com/ifr-acso/casos/wiki#getting-started)
- Available [conic solvers](https://github.com/ifr-acso/casos/wiki#conic-solvers)
- Convex and nonconvex [sum-of-squares optimization](https://github.com/ifr-acso/casos/wiki/sum%E2%80%90of%E2%80%90squares-optimization)
- Supported [vector, matrix, and polynomial cones](https://github.com/ifr-acso/casos/wiki/cones)
- Some [practical tipps](https://github.com/ifr-acso/casos/wiki/practical-sos-guide) to sum-of-squares
- [Transitioning](https://github.com/ifr-acso/casos/wiki/transitioning-from-other-toolboxes) from other toolboxes
- Example [code snippets](https://github.com/ifr-acso/casos/wiki/numerical-examples)

### Cite us

If you use CaΣoS, please cite us:

> T. Cunis and J. Olucak, ‘CaΣoS: A nonlinear sum-of-squares optimization suite’, in _2025 American Control Conference_, (Boulder, CO), pp. 1659–1666, 2025 doi: [10.23919/ACC63710.2025.11107794](https://doi.org/10.23919/ACC63710.2025.11107794).

<details>

<summary>Bibtex entry</summary>

```bibtex
@inproceedings{Cunis2025acc,
	author = {Cunis, Torbjørn and Olucak, Jan},
	title = {{CaΣoS}: {A} nonlinear sum-of-squares optimization suite},
	booktitle = {2025 American Control Conference},
	address = {Boulder, CO},
	year = {2025},
	pages = {1659--1666},
	doi = {10.23919/ACC63710.2025.11107794},
}
```

</details>

----
