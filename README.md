
<p align="center">
<img
  src="https://github.com/user-attachments/assets/372f3428-89e8-4fa8-b5de-385d922d97da#gh-light-mode-only"
  width="50%" />

<img
  src="https://github.com/user-attachments/assets/372f3428-89e8-4fa8-b5de-385d922d97da#gh-dark-mode-only"
  width="50%" />
</p>


----

# CaΣoS-Example package

This repository contains tutorials, examples and implementations from paper or textbooks. The purpose is to provide an easy on-boarding with the CaΣoS toolbox, a nonlinear optimization-oriented sum-of-squares toolbox based on the symbolic framework of CasADi. New examples are continuously added.

### Installation
The example package itself does not not need to be installed. Only CaΣoS and a supported conic solver is needed. Follow the instructions on the [Getting started](https://github.com/ifr-acso/casos/wiki#getting-started) page.

### Folder Structure

```text
Getting Started/				# Tutorials
├── Basics        				# Toolbox basics
├── SOS			  				# How to setup different kinds of SOS problems
├── Conic			  			# How to setup cone problem using sdpsol interface

Systems and Control/ 			# Code from publications or textbooks
├── Stability					# Stability analysis, e.g., region-of-attraction estiamtion
├── Reachability				# Reachability analysis, e.g., inner-approximation of backward reachable set
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
