# Binary Interaction Approximation (BIA) Code Package

[![DOI](https://zenodo.org)](https://doi.org)

A sequential C/GSL verification and benchmark package for the **Binary Interaction Approximation (BIA)** in $N$-body problems.

---

## Authors and Affiliation

* **Author**: Shun-ichi Oikawa (及川 俊一) — Faculty of Engineering, Hokkaido University, Sapporo 060-8628, Japan
* **Contact**: `lopia312oikawa@gmail.com` | ORCID: [0009-0001-5603-2787](https://orcid.org)

---

## Citation

* **Paper**: Shun-ichi Oikawa, *Communications in Nonlinear Science and Numerical Simulation*, Vol. 32, pp. 273–284 (2016). DOI: [10.1016/j.cnsns.2015.11.011](https://doi.org)
* **Zenodo Archives**: Sequential C/GSL (`10.5281/zenodo.24122562`), 3D Multi-Body Hybrid Fortran/CUDA (`10.5281/zenodo.24131822`)

---

## Background, Requirements & Setup

* **Significance**: Analytical reduction of two-body motion in $(r, \theta)$ via GSL ODE solvers; achieves $\sim 10^{-10}$ precision for $N=28$.
* **Dependencies**: GCC / ANSI C, `gfortran`, and GNU Scientific Library (GSL).
* **Usage**: Configure `nbody.inp` (`28, 3e-7`), place `0000028.ccc`, run `make`, and execute `./gsl_dim`. Full parameter and architectural details can be referenced directly in the project files.

---

## Program Catalog & Benchmarks

* **Core & Hermite Solvers**: Includes `gsl_bia`, trajectory export `gsl_dim_out`, and fixed 100-sub-step Hermite verification solvers (`hmt` in C and `hmt_fort` in Fortran).
* **Performance**: Features bit-level agreement between `hmt` and `hmt_fort` ($1.38 \times 10^{-11}$ error), dimensionless plasma scaling ($\Delta t = 1.0$), and configurable `Correction` flags (`-1` to `2`) pushing relative errors down to $10^{-13}$.

---

## Future Roadmap & License

* **GPGPU Scaling**: Self-contained integration-free BIA loops mapped 1-thread-per-particle for massive parallelization (Zenodo DOI: [10.5281/zenodo.24131822](https://doi.org)).
* **License**: GNU General Public License v3.0 (GPL-3.0).
