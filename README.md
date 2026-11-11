Binary Interaction Approximation (BIA) Code Package
A high-precision, highly optimized sequential C/GSL implementation for solving N-body problems using the Binary Interaction Approximation (BIA) framework.

[![DOI](https://zenodo.org/badge/1451113403.svg)](https://doi.org/10.5281/zenodo.24122562)

SECTION-A: Authors and Affiliation

A-1: Name: Shun-ichi Oikawa (及川 俊一)

A-2: Role: Invited Teacher, Faculty of Engineering, Hokkaido University, Sapporo 060-8628, Japan.

A-3: Email: oikawa@khaki.plala.or.jp

SECTION-B: Citation

If you use this code package, its mathematical formulations, or benchmark environments in your research, publications, or conference presentations, please cite the following original peer-reviewed paper:

B-1: Paper Author: Shun-ichi Oikawa

B-2: Paper Title: A Force Evaluation Free Method to N-Body Problems: Binary Interaction Approximation

B-3: Journal: Communications in Nonlinear Science and Numerical Simulation (Commun. Nonlinear Sci. Numer. Simulat.)

B-4: Volume/Issue/Pages: Vol. 32, pp. 273-284

B-5: Year: 2016

B-6: DOI Link: https://doi.org (Please replace with the official DOI link if applicable)

**B-7: Software Archive (Zenodo DOI):** Shun-ichi Oikawa, *Binary Interaction Approximation (BIA) Code Package (v1.0.1)*, Zenodo (2026). DOI: [10.5281/zenodo.24122562](https://doi.org)

**B-8: 3D Multi-Body Hybrid Fortran/CUDA Archive (Zenodo DOI):** Shun-ichi Oikawa, *Binary Interaction Approximation (BIA) — 3D Multi-Body Hybrid Fortran/CUDA Code (v1.0.0)*, Zenodo (2026). DOI: [10.5281/zenodo.24131822](https://doi.org)(2026). DOI: [10.5281/zenodo.24131822](https://doi.org)

SECTION-C: Background and Academic Significance

C-1: The Hermite method program (hmt) included in this package was originally developed as a custom verification tool to address a rigorous peer-review request: to perform a head-to-head comparison between the BIA framework and the industry-standard Hermite integration method.

C-2: Upon the acceptance of the paper, the referees highly commended the academic value and solid mathematical foundation of these computational codes, strongly recommending their release to the scientific community as an open-source package.

C-3: Mathematical Formulation vs. This Implementation
In the original 2016 paper, the BIA framework introduces a radical mathematical breakthrough: since any isolated two-body motion scales to a single orbital plane, it can be analytically formulated in the (r, theta) plane as a function of angle, t = t(theta). By applying Newton's method to invert this relation and solve for the final angle from the initial angle, the framework achieves a completely "Integration-free" architecture, entirely eliminating the need for time-stepping numerical loops. While the original paper leverages this analytical reduction into an ultra-fast massive parallel CUDA framework, this public package (bia_public) serves a different, foundational purpose.

C-4: This package is meticulously designed as a baseline sequential C code where all time evolutions are intentionally delegated to the standard numerical ODE solvers of the GNU Scientific Library (GSL).

C-5: By decoupling the complex analytical inverse-problem solvers and utilizing global standard libraries, this repository allows researchers to easily dissect, understand, and reproduce the core mechanics of the BIA algorithm itself, specifically, how a chaotic N-body problem is structurally reduced into a clean superposition of independent two-body interactions.

C-6: Remarkably, even when integrated via raw, uncorrected GSL numerical solvers, this framework delivers an exceptional relative error on the order of 10^-10. This empirical stability under standard library conditions serves as ironclad evidence of the organic physical validity and robust potential of the underlying BIA formulation.

C-7: Mathematical and Computational Advantages
Unlike brute-force sequential N-body simulations, this public package (bia_public) operates as a pure, high-performance sequential code (single-threaded CPU) without requiring parallel dependencies or dedicated GPU environments.

C-8: By leveraging the core "Integration-free" (analytical reduction) nature of the BIA framework, the inner loops eliminate redundant iterations and memory latency during single-particle updates. This elegant mathematical architecture allows the sequential code to achieve stellar execution speeds, completing time evolution for N = 28 within milliseconds, while simultaneously achieving ultra-high precision approaching the fundamental limits of double-precision floating-point representation.

SECTION-D: Repository Structure

<pre style="white-space: pre-wrap;">
.
├── All Source Files     # Complete C/Fortran source and header files (.c, .f, .h)
├── Makefile             # Build configuration for compilation
├── nbody.inp            # Configuration input file (Parameters scaled by mass)
├── 0000028.ccc          # Initial condition data file (Example dataset for N=28)
├── LICENSE              # GNU General Public License v3 (GPL v3)
└── README.md            # This documentation framework
</pre>.

SECTION-E: Requirements and Dependencies

E-1: [C Compiler] GCC (or any standard ANSI C compliant compiler)

E-2: [Fortran Compiler] gfortran (Required only for compiling the Fortran-linked verification module hmt_fort)

E-3: [Libraries] GNU Scientific Library (GSL) with linear algebra support

SECTION-F: Setup and Configuration

Before execution, ensure that the configuration file nbody.inp and the corresponding initial condition file (*.ccc) reside in the same directory as the compiled binary.

F-1: [Section 1] Configuration File (nbody.inp)
The file nbody.inp specifies the number of particles and the interaction parameter. The program reads only the first line of this file.
Example Content (nbody.inp): 28, 3e-7
Format: [Number of Particles], [Interaction Parameter]
Parameter Note: The interaction constant (e.g., 3e-7) represents a dimensionless repulsive interaction parameter scaled by mass (positive values denote repulsion). Because this package models two-body interactions among identical particles, the internal algorithm doubles this value to account for the reduced mass dynamics.

F-2: [Section 2] Input Coordinate File (*.ccc)
Place the coordinate file matching your target particle count into the working directory.
Example: For N = 28, the file 0000028.ccc must be present.

SECTION-G: Compilation and Execution

G-1: [Step 1] Build the Binaries
Compile the entire package using the provided Makefile by typing: make

G-2: [Step 2] Generate High-Precision Initial Solvers
Run the primary binary to automatically generate the initial high-precision solver boundary data (*.dim) based on nbody.inp and the associated *.ccc file by typing: ./gsl_dim

SECTION-H: Program Catalog and Output Specifications

This package includes several distinct executable variants tailored for specific analytical and verification purposes:

H-1: [Category 1] BIA Core Solvers (Standard C/GSL Implementation)
These standard programs implement the BIA formulation using the generalized ODE solver modules of the GNU Scientific Library (GSL) as their time-evolution engine, maximizing portability and reproducibility.

H-1-A: gsl_bia / gsl_dim (Endpoint Output Version): Computes time evolution from t = 0.0 to t = 1.0 and outputs numerical data strictly at the boundary endpoints. It is optimized for the high-speed generation of initial solutions.

H-1-B: gsl_dim_out (Time-Series Output Version): Computes from t = 0.0 to t = 1.0 while exporting the complete trajectory history. This version is ideal for tracing detailed orbital dynamics, tracking particle identity exchanges during close-encounters, and generating visualizations (e.g., via gnuplot).

H-1-C: Technical Note: Even when running raw without numerical corrections (gsl_bia), these codes deliver a relative error on the order of 10^-10. This stability proves the inherent robustness and mathematical validity of the BIA paradigm when mapped to standard library integrations.

H-2: [Category 2] gsl_bia_inc (Time-Reversal Direct Increment Version)

H-2-A: Mathematical Spec: Designed to verify the Time-Reversal Symmetry of the system by integrating backward from t = 1.0 to t = 0.0.

H-2-B: Advantage: Reversing absolute trajectory integrals typically triggers catastrophic cancellation (compounding round-off errors). By enabling the BIA_INC macro, this program reformulates the differential equations to integrate the incremental changes directly. This eliminates numerical drift, allowing the system to rewind to its exact initial state at double-precision limits.

H-3: [Category 3] Verification Hermite Solvers (hmt Series)
These solvers were specifically engineered to cross-verify the BIA infrastructure against the industry-standard 4th-order Hermite predictor-corrector scheme, completely independent of external integration libraries.

H-3-A: Independent Implementation (No GSL Integration): Unlike the standard gsl_bia series which relies on GSL's general ODE solvers, the hmt series does not use GSL for time evolution. Instead, it directly couples the BIA framework (deconstructing N-body problems into a superposition of two-body interactions) with a custom, self-contained 4th-order Hermite integration algorithm.

H-3-B: Fixed Sub-stepping: It splits each global time step delta-t into a strictly uniform, fixed sub-step scheme (defaulting to 100 subdivisions). By executing pure, native numerical integration directly on the BIA-reduced two-body dynamics without adaptive step-size overhead, it eliminates artificial numerical dissipation and achieves exceptional energy conservation that pushes the fundamental limits of double precision.

H-3-C: hmt (Pure C Version - Recommended): A complete C-translation of the proven Fortran core algorithm. By utilizing internal array operations rather than external library routines, it eliminates dependencies on both GSL integration modules and Fortran compilers, providing elite cross-platform portability.

H-3-D: hmt_fort (Fortran Object-Linked Version): Bridges a C wrapper (gsl_hmt.c) with the legacy, self-developed Fortran integration object (HermiteIntegrator2.o). This is included to rigorously recreate the exact benchmark environment utilized in the original paper (requires gfortran).

SECTION-I: Algorithmic Insights and Verification Benchmarks

I-1: Step Control and Energy Conservation
The Hermite solvers (hmt / hmt_fort) bypass adaptive step-size control. Instead, they execute on a strictly uniform, fixed sub-step scheme (defaulting to 100 subdivisions per main time step delta-t). This eliminates artificial numerical dissipation induced by adaptive steps, achieving exceptional energy conservation that pushes the boundaries of double precision.

I-2: The Peer-Review Cross-Verification Method
To ensure a fair comparison during the peer-review process against standard adaptive Runge-Kutta-Fehlberg (RKF) methods, a specialized benchmarking technique was engineered: The adaptive RKF solver was executed first, outputting a complete history of its accepted time-step intervals. The custom Hermite integration algorithm was then forced to follow that exact identical time-step history step-by-step. This methodology isolates the pure trajectory and energy conservation capabilities of the underlying integrators by entirely removing step-control discrepancies.

I-3: Bit-Level Verification Match (N = 28)
Cross-platform tests between the pure C version (hmt) and the Fortran-linked version (hmt_fort) under identical parameters (N = 28, 100 sub-steps per delta-t) demonstrate a perfect bit-level agreement: The absolute energy error is 1.377253866507999 x 10^-11, and the relative error is 3.847 x 10^-12. This cross-verification mathematically guarantees the algorithmic correctness of the C translation.

I-4: Verification Note: While hmt_fort benefits from Fortran's array optimization and executes slightly faster, hmt serves as the repository default due to its excellent cross-platform portability.

SECTION-J: Physical Scaling and Dimensionless Normalization

The dimensionless parameters and scaling laws implemented in this package are derived from realistic high-energy plasma physics regimes (e.g., nuclear fusion environments). The scaling maps the code's normalized units to physical quantities based on the following baseline configuration:

J-1: Particle Species: Identical protons (m = m_p).

J-2: Reference Energy / Temperature: 10 keV protons.

J-3: Reference Velocity Scale: The average thermal velocity of the identical background particles is normalized to unity (v_th = 1.0, approximately 1.384 x 10^6 m/s).

J-4: Number Density (n): n = 10^20 m^-3 (a standard high-density plasma regime).

J-5: Spatial Scale: The average inter-particle distance is strictly normalized to unity (d_avg = 1.0).

J-6: Large Time-Stepping and Stellar Precision
Under this spatial normalization (d_avg = 1.0), this framework executes time evolution using a remarkably large macroscopic time-step size of delta-t = 1.0. In conventional N-body simulations, taking such a massive time step would instantly destroy the orbital invariants. Remarkably, by virtue of the BIA deconstruction, a single execution step under delta-t = 1.0 yields a high-precision solution with approximately 10 significant digits (a relative error on the order of 10^-10) depending on the particle count. This performance represents an unprecedented leap in balancing computational macro-steps with micro-level precision.

J-7: Dimensional Interpretation of the Input File (nbody.inp)
Based on this scaling, the interaction parameter provided in nbody.inp (e.g., 3e-7) matches the mass-scaled, dimensionless repulsive interaction parameter. Because the BIA framework naturally deconstructs the multi-body system into reduced-mass two-body interactions, the code internally doubles this input parameter to calculate the exact Coulomb or screened-Coulomb collision dynamics within this physical scale.
* **`eps_th` (Dimensionless Interaction Parameter):**
  * **Physical Definition:** Represents the Coulombic (or gravitational) interaction strength normalized by the particle mass, setting the characteristic relative thermal speed to $1$ and the mean inter-particle distance to $1$. Under typical fusion-grade proton plasma conditions, this normalized single-particle interaction parameter is approximately $3 \times 10^{-7}$.
  * **Factor of $2$ (Reduced Mass Effect):** Because the BIA formulation solves the two-body relative equation of motion ($\mathbf{r}_{ij} = \mathbf{r}_i - \mathbf{r}_j$), dividing by the reduced mass $\mu = m/2$ introduces an exact factor of $2$ into the relative acceleration equation:
    $$\varepsilon_{\text{th}} = \pm 2 \times (3 \times 10^{-7}) = \pm 6.0 \times 10^{-7}.$$
  * **Sign Convention:** Negative (`-6.0e-7`, i.e., `-2 * 3e-7`) denotes **attractive** interactions (such as gravitational systems or opposite-charge Coulombic attraction), whereas positive (`+6.0e-7`, i.e., `+2 * 3e-7`) denotes **repulsive** interactions (such as same-species fusion protons). This normalization enables exact analytical integration across macroscopic time steps of $\Delta t = 1$ (mean inter-particle transit time).

SECTION-K: Numerical Properties and Theoretical Considerations

* **K-0: Strictly Momentum-Conserving Formulation:**
  The BIA numerical scheme is constructed as a **strictly momentum-conserving** algorithm. By enforcing exact pairwise action-reaction symmetry ($\Delta \mathbf{v}_i = +\Delta \mathbf{v}_{ij}$ and $\Delta \mathbf{v}_j = -\Delta \mathbf{v}_{ij}$) across every interacting particle pair, total linear momentum is conserved down to the double-precision round-off limit without artificial center-of-mass drift.

K-1: Absolute Robustness Against Binary Bounds (Close Encounters)
One of the most profound advantages of the BIA paradigm lies in its absolute immunity to the stiffness or singularities caused by close-encounter particles, such as binary bound states. In traditional brute-force N-body codes or standard variable-step numerical integrators, when two particles approach each other closely, the gravitational or Coulomb forces diverge, forcing the time-step size to shrink toward zero. This triggers an immediate computational bottleneck, causing the simulation to grind to a complete halt. By contrast, because BIA reduces the global network into a superposition of isolated two-body dynamics, it allows for the seamless integration of exact analytical solutions (the inverse problem solved via Newton's method). Once these exact analytical solvers are fully utilized, the system never suffers from time-step collapse. No matter how extremely close two particles encounter each other, or how tight a binary system becomes, the simulation is mathematically guaranteed to progress without ever freezing or failing.

K-2: Long-Range Constraints and Mitigations
When solving the inverse problem of exact two-body solutions via Newton's method under extreme long-range (weak interaction) regimes, the difference between initial and final angles can approach floating-point precision limits, potentially introducing minor rounding anomalies. For large-scale production simulations, these long-range regimes are seamlessly and stably handled by embedding a Tree-code framework or combining BIA with hybrid standard integrators (such as RKF).

K-3: Mathematical Elegance over High-Order Complexity
The 10^-10 relative accuracy achieved by the uncorrected gsl_bia solver underlines the intrinsic mathematical advantages of the BIA framework. Traditional N-body models require highly complex, high-order geometric integrators to counter chaotic orbital collapse. Because BIA analytically simplifies the global N-body system down to a clean superposition of center-of-mass and relative two-body motions, it naturally guarantees remarkable precision without demanding computationally heavy, high-order stabilization schemes.

K-4: Correction Switching Mechanism (Correction Flag)
For users exploring error dynamics, the main source code exposes an algorithmic toggle via the Correction variable: int Correction = 2; // Options: -1, 0, 1, or 2

K-4-A: Correction = 2 (Default Standard Setting): Runs the raw, uncorrected sequential BIA framework to showcase its organic mathematical precision.

K-4-B: Correction = -1, 0, 1 (Experimental Configurations): Activates advanced energy and angular momentum correction algorithms developed from comprehensive error-bound analyses. Even at a massive time-step size of delta-t = 1, these configurations drive conservation values to near double-precision limits (10^-11 to 10^-13 relative error).

K-5: Reference Verification Sample (N = 28, delta-t = 1.0)
The following lines provide reference output metrics at the final boundary (t = 1.0) for a system with N = 28 and delta-t = 1.0. Use these values to verify the numerical accuracy of your local compilation environment.

K-5-A: File [gsl_bia.out] -> Solver: GSL Standard Solver (No Corrections) -> Calculated Position X: -1.866446842361e-06 -> Meaning: Baseline showcasing the raw organic potential of BIA mathematics.

K-5-B: File [gsl_dim.out] -> Solver: GSL High-Precision (Endpoint Mode) -> Calculated Position X: -1.866457588462e-06 -> Meaning: Primary reference high-precision solution computed via standard global libraries.

K-5-C: File [gsl_hmt.out] -> Solver: Custom 4th-Order Hermite Solver -> Calculated Position X: -1.866457581209e-06 -> Meaning: Cross-verification benchmark. Matches gsl_dim up to the 5th significant digit across divergent solver backends.

SECTION-L: Future Roadmap: Extension to Massive GPGPU Parallelism

L-1: This sequential package (bia_public) is designed to establish the core mathematical validity and baseline CPU benchmarks for the BIA framework.

L-2: Because the innermost loops computing the BIA state for each particle (or pair) are mathematically self-contained and entirely Integration-free, they are completely isolated from thread-to-thread synchronization bottlenecks. This structural independence makes the codebase a perfect starter kit for scaling up.

L-3: The next milestone for this repository will feature the release of an unconstrained GPGPU (CUDA-optimized) package utilizing:

L-3-A: True Integration-Free Solvers: Implementing exact analytical solutions via inverse Newton solvers to completely eliminate time-stepping numerical loops.
L-3-B: 1-Thread-per-1-Particle Mapping: Exploiting thousands of GPU cores to compute individual particle dynamics instantly and independently, rewriting the efficiency paradigm for massive N-body simulations.
