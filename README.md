# Binary Interaction Approximation (BIA) Code Package

[![DOI](https://zenodo.org)](https://doi.org)

## SECTION-A: Authors and Affiliation
* **A-1 to A-3:** Shun-ichi Oikawa (Faculty of Engineering, Hokkaido University, Sapporo, Japan; oikawa@khaki.plala.or.jp)

## SECTION-B: Citation
* **B-1 to B-6:** Shun-ichi Oikawa, *A Force Evaluation Free Method to N-Body Problems: Binary Interaction Approximation*, Communications in Nonlinear Science and Numerical Simulation, Vol. 32, pp. 273–284 (2016). DOI: [10.1016/j.cnsns.2015.11.011](https://doi.org)
* **B-7 & B-8:** Software archives available on Zenodo:
  * Sequential C/GSL Package (v1.0.1): [10.5281/zenodo.24122562](https://doi.org)
  * 3D Multi-Body Hybrid Fortran/CUDA Code (v1.0.0): [10.5281/zenodo.24131822](https://doi.org)

## SECTION-C: Background and Academic Significance
* **C-1 & C-2:** Originally developed as a verification tool to compare BIA against the standard Hermite integration method following peer-review recommendations.
* **C-3 to C-5:** Implements an analytical reduction of two-body motion in the \((r, \theta)\) plane using GSL ODE solvers, allowing researchers to study N-body interactions as independent two-body superpositions.
* **C-6 to C-8:** Achieves high relative precision (\(\sim 10^{-10}\)) and fast sequential execution for \(N=28\) without requiring GPU or parallel dependencies.

## SECTION-D: Repository Structure
```text
.
├── All Source Files     # Complete C/Fortran source and header files (.c, .f, .h)
├── Makefile             # Build configuration for compilation
├── nbody.inp            # Configuration input file (Parameters scaled by mass)
├── 0000028.ccc          # Initial condition data file (Example dataset for N=28)
├── LICENSE              # GNU General Public License v3 (GPL v3)
└── README.md            # This documentation framework
```

## SECTION-E: Requirements and Dependencies
* **E-1:** GCC (or ANSI C compiler)
* **E-2:** `gfortran` (for Fortran-linked verification module `hmt_fort`)
* **E-3:** GNU Scientific Library (GSL) with linear algebra support

## SECTION-F: Setup and Configuration
* **F-1:** Configure `nbody.inp` with `[Number of Particles], [Interaction Parameter]` (e.g., `28, 3e-7`).
* **F-2:** Place the matching coordinate file (e.g., `0000028.ccc` for $N=28$) in the working directory.

## SECTION-G: Compilation and Execution
* **G-1:** Build binaries using `make`.
* **G-2:** Run `./gsl_dim` to generate initial high-precision solver boundary data (`*.dim`).

## SECTION-H: Program Catalog and Output Specifications

This package includes several distinct executable variants tailored for specific analytical and verification purposes:

* **H-1: [Category 1] BIA Core Solvers (Standard C/GSL Implementation)**: Includes endpoints (`gsl_bia`), full time-series trajectory exports (`gsl_dim_out`), and achieves $10^{-10}$ relative error even in raw runs.
* **H-2: [Category 2] Verification Hermite Solvers (`hmt` Series)**: Uses custom 4th-order Hermite integration independent of GSL, featuring fixed 100-sub-step schemes to eliminate artificial numerical dissipation. Includes pure C (`hmt`, recommended) and Fortran-linked (`hmt_fort`) versions.

---

## SECTION-I: Algorithmic Insights and Verification Benchmarks

* **I-1 & I-2**: BIA uses fixed sub-stepping (100 subdivisions per $\Delta t$) to bypass adaptive step-size dissipation. A specialized peer-review technique forced the Hermite solver to follow an adaptive RKF solver's exact time-step history to isolate pure trajectory and energy performance.
* **I-3 & I-4**: Cross-platform tests ($N = 28$) between pure C (`hmt`) and Fortran-linked (`hmt_fort`) achieve bit-level agreement (absolute energy error $1.38 \times 10^{-11}$). While `hmt_fort` executes slightly faster due to Fortran array optimizations, `hmt` is the repository default for cross-platform portability.

---

## SECTION-J: Physical Scaling and Dimensionless Normalization

Derived from high-energy fusion plasma regimes ($10\text{ keV}$ protons, $n = 10^{20}\text{ m}^{-3}$, $v_{\text{th}} = 1.0$, $d_{\text{avg}} = 1.0$), this framework executes macroscopic time steps of $\Delta t = 1.0$ while maintaining $\sim 10$ significant digits of precision.
* **J-7 (`eps_th`):** The input parameter in `nbody.inp` (e.g., `3e-7`) is internally doubled due to reduced mass $\mu = m/2$ in the two-body equation ($\varepsilon_{\text{th}} = \pm 6.0 \times 10^{-7}$). Negative values indicate attractive interactions, and positive values indicate repulsive interactions.

---

## SECTION-K: Numerical Properties and Theoretical Considerations

* **K-0 & K-1**: Enforces strict pairwise action-reaction symmetry for zero center-of-mass drift, and avoids time-step collapse during close encounters or binary bound states by utilizing exact analytical inverse solutions via Newton's method.
* **K-2 & K-3**: Long-range interaction limits can be stabilized via Tree-code or hybrid RKF integration. The method achieves $10^{-10}$ relative accuracy without heavy, high-order geometric integrators.
* **K-4 (`Correction` Flag)**: Exposes a toggle (`Correction` = `-1, 0, 1, or 2`). Default `Correction = 2` runs raw uncorrected BIA; `-1, 0, 1` activate experimental conservation algorithms pushing relative error down to $10^{-11} \sim 10^{-13}$.
* **K-5**: Reference verification metrics at $t = 1.0$ ($N = 28, \Delta t = 1.0$) for `gsl_bia.out` (`-1.866446842361e-06`), `gsl_dim.out` (`-1.866457588462e-06`), and `gsl_hmt.out` (`-1.866457581209e-06`).

---

## SECTION-L: Future Roadmap: Extension to Massive GPGPU Parallelism

* **L-1 & L-2**: The integration-free, self-contained BIA loops eliminate thread synchronization bottlenecks, serving as an ideal foundation for scaling.
* **L-3**: Archived on Zenodo ([DOI: 10.5281/zenodo.24131822](https://doi.org)), featuring a 3D Multi-Body Hybrid Fortran/CUDA implementation with true integration-free inverse Newton solvers and 1-thread-per-particle mapping.
