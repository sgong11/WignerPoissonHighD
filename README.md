# The Wigner–Poisson Model

This repository is being developed as a home for numerical experiments and reproducible results for the Wigner–Poisson (WP) system and adaptive-rank solvers for quantum kinetic plasma dynamics.  The current focus is the high-dimensional, structure-preserving adaptive-rank Wigner–Poisson solver in hierarchical Tucker / HTACA form.  Over time, this repository will also collect results and reproducibility material for our 1D1V adaptive-rank work, KEEN-wave simulations, conservative low-rank methods, and related Wigner–Poisson test problems.

A separate README will describe the details of the HTACA implementation and code workflow.  This README is intended to explain the model, why it is useful, where it appears in the literature, and why it is computationally challenging.

---

## 1. What is the Wigner–Poisson system?

The Wigner–Poisson system is a deterministic phase-space model for quantum kinetic electron dynamics.  It evolves a Wigner distribution

```math
f = f(x,v,t),
```

where $x\in\mathbb{R}^d$ is position, $v\in\mathbb{R}^d$ is velocity, and $t$ is time.  The Wigner distribution plays a role analogous to the distribution function in classical kinetic theory, but it is not a classical probability density: it may take negative values.  Those negative regions encode genuinely quantum behavior such as interference and non-classical phase-space structure.

The distribution is coupled to an electrostatic potential $\Phi(x,t)$, which is determined self-consistently through Poisson's equation.  In nondimensional form, a common periodic Wigner–Poisson model is

```math
\partial_t f + v\cdot \nabla_x f
=
-\frac{i}{(2\pi)^d H^{d+1}}
\int_{\mathbb{R}^d}\int_{\mathbb{R}^d}
\exp\!\left(\frac{i(v'-v)\cdot y}{H}\right)
\left[\Phi\!\left(x+\frac{y}{2}\right)-\Phi\!\left(x-\frac{y}{2}\right)\right]
f(x,v',t)\,dy\,dv',
```

```math
-\Delta_x \Phi = \int_{\mathbb{R}^d} f(x,v,t)\,dv - 1.
```

Here $H$ is the nondimensional quantum parameter.  With the usual plasma scaling by the Debye length $\lambda_D$ and plasma frequency $\omega_{pe}$,

```math
H = \frac{\hbar}{m_e \lambda_D^2 \omega_{pe}}.
```

The parameter $H$ measures the strength of quantum diffraction and uncertainty effects.  Formally, as $H\to 0$, the Wigner–Poisson model approaches the classical Vlasov–Poisson model.  For finite $H$, the model retains the phase-space structure of kinetic theory while adding quantum wave effects.

---

## 2. Why use Wigner–Poisson instead of a classical kinetic model?

Classical Vlasov–Poisson dynamics describes collisionless electrostatic plasma behavior using a local force term $E\cdot\nabla_v f$.  Wigner–Poisson replaces that local force term by a nonlocal pseudodifferential Wigner operator.  This nonlocal operator allows the model to describe effects that are absent from classical kinetic theory, including

- quantum diffraction,
- wave-packet spreading,
- tunneling-like behavior,
- interference fringes in phase space,
- sign-changing Wigner structures,
- quantum modifications of wave-particle interaction.

This makes WP useful in regimes where electrons are hot enough that a fully degenerate Fermi model may not be required, but dense or cold enough that the electron de Broglie wavelength is not negligible compared with the kinetic or electrostatic length scales.

Examples include warm dense matter, dense plasmas, resonant tunneling and semiconductor devices, superlattices, 2D electron systems, driven plasmonic systems, and quantum corrections to nonlinear electrostatic plasma waves.

---

## 3. Conserved quantities

On periodic domains, the Wigner–Poisson system has the same basic global invariants one expects from a collisionless electrostatic kinetic model.  With

```math
\rho(x,t)=\int f(x,v,t)\,dv,
\qquad
J(x,t)=\int v f(x,v,t)\,dv,
```

the global invariants are

```math
M[f] = \int \rho(x,t)\,dx,
```

```math
P[f] = \int J(x,t)\,dx,
```

and the self-consistent total energy

```math
\mathcal{E}[f]
=
\frac12 \int\!\int |v|^2 f(x,v,t)\,dv\,dx
+
\frac12 \int |\nabla_x \Phi(x,t)|^2\,dx.
```

For numerical simulation, preserving these quantities is important.  Small conservation defects can accumulate over long times and contaminate wave dynamics, phase-space structures, and electric-energy diagnostics.  This is especially important for low-rank and adaptive methods, where compression and recompression can destroy invariants even when the corresponding full-grid method is conservative.

---

## 4. Fourier form of the Wigner operator

A key reason spectral methods are natural for WP is that the Wigner potential operator becomes simple in velocity-Fourier variables.  If $\widehat f(x,k_v,t)$ denotes the Fourier transform of $f$ in velocity, then the Wigner update can be written as

```math
\partial_t \widehat f(x,k_v,t)
=
\frac{i}{H}
\left[
\Phi\!\left(x+\frac{H k_v}{2}\right)
-
\Phi\!\left(x-\frac{H k_v}{2}\right)
\right]
\widehat f(x,k_v,t).
```

With $\Phi$ frozen over a time step, this equation has an analytic exponential update.  This observation is central in many Wigner solvers, including the methods used in our WP papers.

However, this Fourier formulation introduces an additional structure-preservation requirement: because the physical Wigner distribution should be real-valued, the velocity-Fourier representation must preserve Hermitian symmetry.  Low-rank compression and adaptive sampling do not automatically preserve that symmetry, so the numerical method must enforce it deliberately.

---

## 5. Applications of the model

### Quantum transport and semiconductor devices

The Wigner equation and Wigner–Poisson system have long been used for quantum transport in resonant tunneling diodes, semiconductor devices, and superlattices.  In this setting, WP provides a phase-space description of electron transport where nonlocal quantum effects influence tunneling, transient response, current oscillations, and boundary-driven behavior.

Useful entry points include Frensley's Wigner-function model for resonant tunneling devices, Ringhofer's spectral Wigner and Wigner–Poisson methods, and the semiconductor quantum-transport literature of Markowich, Ringhofer, Schmeiser, Arnold, and collaborators.

### Quantum plasmas and warm dense matter

In warm dense matter and high-energy-density plasma regimes, the electron de Broglie wavelength can become comparable to collective plasma scales.  Wigner–Poisson is a reduced mean-field model that retains kinetic phase-space information while including quantum diffraction and uncertainty.  This makes it useful for studying quantum modifications to electrostatic waves, dielectric response, and stopping-power-related physics.

The high-dimensional simulations in this repository are motivated in part by the need for kinetic models of electron dynamics in regimes relevant to inertial-confinement fusion and alpha-particle stopping-power studies.

### KEEN waves and driven nonlinear electrostatic dynamics

Kinetic Electrostatic Electron Nonlinear (KEEN) waves are nonlinear, driven, subplasma-frequency electrostatic structures that depend on trapping and multiharmonic phase-space organization.  Our KEEN-wave WP study uses a 1D1V Wigner–Poisson solver to examine how quantum diffraction changes this classical mechanism.  As $H$ increases, quantum diffraction weakens trapping, damps higher harmonics, diffuses trapped-electron vortices, and drives the system toward a lower long-time electrostatic-energy state.

This makes KEEN waves a useful reduced test problem for understanding how quantum kinetic effects alter driven, nonlinear, nonequilibrium plasma dynamics.

### Collisional extensions

The collisionless WP model is also a starting point for Wigner–Poisson–BGK and other collisional quantum kinetic models.  These extensions are relevant when relaxation, diffusion limits, and dissipative transport must be included.  The adaptive-rank and conservation ideas in this repository are intended to be compatible with future extensions of this type.

---

## 6. Why Wigner–Poisson is computationally hard

Wigner–Poisson is appealing because it keeps a kinetic phase-space description while adding quantum effects.  That same structure creates several major numerical challenges.

### 6.1 Full phase space is large

A $d$-dimensional physical problem lives in a $2d$-dimensional phase space.  Thus,

- 1D1V is a two-dimensional array,
- 2D2V is a four-dimensional tensor,
- 3D3V is a six-dimensional tensor.

A full tensor with $N$ grid points per coordinate has $N^{2d}$ entries.  This becomes prohibitive in 3D3V.

### 6.2 The Wigner potential is nonlocal and oscillatory

The Wigner operator couples the potential at shifted spatial points with the distribution over velocity.  Unlike the Vlasov force term, it is not simply a local derivative in velocity.  Direct evaluation is expensive, communication-heavy, and difficult to scale in parallel.

### 6.3 The solution is oscillatory and sign-indefinite

The Wigner distribution can be negative.  Interference fringes and phase-space oscillations are not numerical artifacts; they are part of the model.  Any numerical method must distinguish physical oscillations from unresolved noise.  Smaller $H$ typically requires finer resolution and/or higher rank because the dynamics become more classical and filamentary.

### 6.4 Fourier updates must preserve Hermitian symmetry

The analytic velocity-Fourier update is efficient, but the discrete representation must preserve conjugate symmetry so that the inverse transform remains real-valued.  This is straightforward in a full tensor when all paired modes are available, but it is subtle under adaptive sampling and low-rank recompression.

### 6.5 Compression can break conservation

Low-rank methods reduce memory and runtime by approximating the solution in a compressed representation.  But truncation can break mass, momentum, and total-energy conservation.  A useful low-rank WP method must therefore preserve both the numerical efficiency and the physical structure of the model.

### 6.6 Self-consistent field coupling matters

The potential $\Phi$ is computed from the density through Poisson's equation.  Errors in density affect the field, and errors in the field feed back into the Wigner update.  The kinetic energy and electrostatic field energy must be treated together when enforcing total energy.

---

## 7. Why adaptive rank and HTACA are used here

The central observation behind this repository is that many WP solutions have exploitable low-rank structure, especially for finite quantum parameter $H$.  Rather than assembling the entire phase-space tensor, the high-dimensional solver uses hierarchical Tucker adaptive cross approximation (HTACA) to sample selected tensor entries and build a compressed hierarchical Tucker representation.

At a high level, HTACA helps because it can

- represent high-dimensional phase-space data without forming the full tensor,
- exploit low-rank structure in the solution,
- reduce storage and computational cost when ranks remain moderate,
- enable 2D2V and 3D3V simulations that are otherwise infeasible on a full grid.

For WP, compression alone is not enough.  The method in this repository is designed to be structure-preserving:

1. **Hermitian-symmetry-aware Fourier sampling:** paired velocity-Fourier modes are handled so that the inverse transform remains real-valued up to roundoff.
2. **Global conservation correction:** mass, momentum, and self-consistent total energy are corrected after adaptive compression.
3. **HTACA formulation:** the solver samples and recompresses tensor data in hierarchical Tucker form without assembling the full phase-space tensor.

A separate method README will describe the HTACA algorithm, tensor tree, sampling interface, Fourier update, and conservation correction in code-level detail.

---

## 8. Papers connected to this repository

The repository will grow to include results from several related WP papers.

| Paper | Status | Role in this repository |
|---|---|---|
| Andrew J. Christlieb, Sining Gong, Jing-Mei Qiu, and Nanyi Zheng, **A Sampling-Based Adaptive Rank Approach to the Wigner–Poisson System** | Accepted in *SIAM Journal on Scientific Computing*; preprint: [arXiv:2506.21314](https://arxiv.org/abs/2506.21314) | First adaptive-rank WP paper; 1D1V ACA/SLAR-Fourier formulation and structure-preserving Fourier update. |
| F. Alejandro Padilla-Gomez, Sining Gong, Michael S. Murillo, F. R. Graziani, and Andrew J. Christlieb, **Quantum Kinetic Modeling of KEEN Waves in a Warm-Dense Regime** | Published in *Physics of Plasmas*, 33, 043902 (2026), DOI: [10.1063/5.0308425](https://doi.org/10.1063/5.0308425) | Application paper showing how WP quantum diffraction modifies KEEN-wave trapping, harmonics, and long-time decay. |
| Andrew J. Christlieb, Sining Gong, Jing-Mei Qiu, and Nanyi Zheng, **A Structure-Preserving Adaptive-Rank Approach to the High-Dimensional Wigner–Poisson System** | Preprint on arXiv: 
https://doi.org/10.48550/arXiv.2606.15067 | Main high-dimensional HTACA WP method; 2D2V and 3D3V simulations with symmetry and conservation correction. |
| Andrew J. Christlieb, Sining Gong, F. Alejandro Padilla-Gomez, and Jing-Mei Qiu, **A Conservative Adaptive Low-Rank Method for the Wigner–Poisson System** | Preprint forthcoming on arXiv | Conservative 1D1V low-rank WP method using a macroscopic correction, Fermi–Dirac-type reconstruction, and global moment correction. |

---

## 9. Suggested reading

The following references are useful entry points for understanding the model, applications, and numerical methods.

### Foundations and model background

- E. Wigner, “On the quantum correction for thermodynamic equilibrium,” *Physical Review*, 40, 749 (1932).
- P. A. Markowich, C. A. Ringhofer, and C. Schmeiser, *Semiconductor Equations*, Springer (2012).
- F. Haas, *Quantum Plasmas: An Hydrodynamic Approach*, Springer (2011).
- M. Bonitz, *Quantum Kinetic Theory*, Springer (2016).
- J. Weinbub and D. K. Ferry, “Recent advances in Wigner function approaches,” *Applied Physics Reviews*, 5, 041104 (2018).

### Early and classical numerical Wigner–Poisson methods

- W. R. Frensley, “Wigner-function model of a resonant-tunneling semiconductor device,” *Physical Review B*, 36, 1570 (1987).
- C. Ringhofer, “A spectral method for the numerical simulation of quantum tunneling phenomena,” *SIAM Journal on Numerical Analysis*, 27, 32–50 (1990).
- C. Ringhofer, “A spectral collocation technique for the solution of the Wigner–Poisson problem,” *SIAM Journal on Numerical Analysis*, 29, 679–700 (1992).
- A. Arnold and C. Ringhofer, “Operator splitting methods applied to spectral discretizations of quantum transport equations,” *SIAM Journal on Numerical Analysis*, 32, 1876–1894 (1995).
- A. Arnold and C. Ringhofer, “An operator splitting method for the Wigner–Poisson problem,” *SIAM Journal on Numerical Analysis*, 33, 1622–1643 (1996).
- S. Shao, T. Lu, and W. Cai, “Adaptive conservative cell average spectral element methods for transient Wigner equation in quantum transport,” *Communications in Computational Physics*, 9, 711–739 (2011).
- S. Shao and J. M. Sellier, “Comparison of deterministic and stochastic methods for time-dependent Wigner simulations,” *Journal of Computational Physics*, 300, 167–185 (2015).

### Low-rank and tensor methods relevant to this code

- O. Koch and C. Lubich, “Dynamical low-rank approximation,” *SIAM Journal on Matrix Analysis and Applications*, 29, 434–454 (2007).
- C. Lubich and I. V. Oseledets, “A projector-splitting integrator for dynamical low-rank approximation,” *BIT Numerical Mathematics*, 54, 171–188 (2014).
- I. V. Oseledets, “Tensor-train decomposition,” *SIAM Journal on Scientific Computing*, 33, 2295–2317 (2011).
- W. Hackbusch, *Tensor Spaces and Numerical Tensor Calculus*, Springer (2012).
- J. Ballani and L. Grasedyck, “Tree adaptive approximation in the hierarchical tensor format,” *SIAM Journal on Scientific Computing*, 36, A1415–A1431 (2014).
- W. Guo and J.-M. Qiu, “A conservative low rank tensor method for the Vlasov dynamics,” *SIAM Journal on Scientific Computing*, 46, A232–A263 (2024).
- W. Guo and J.-M. Qiu, “A local macroscopic conservative (LoMaC) low rank tensor method for the Vlasov dynamics,” *Journal of Scientific Computing*, 101, 61 (2024).

### Applications to quantum plasmas, KEEN waves, and warm dense matter

- F. R. Graziani, J. D. Bauer, and M. S. Murillo, “Kinetic theory molecular dynamics and hot dense matter: Theoretical foundations,” *Physical Review E*, 90, 033104 (2014).
- T.-X. Hu, J.-H. Liang, Z.-M. Sheng, and D. Wu, “Kinetic investigations of nonlinear electrostatic excitations in quantum plasmas,” *Physical Review E*, 105, 065203 (2022).
- J. T. Mendonça, “Landau damping and particle trapping in the quantum regime,” *Reviews of Modern Plasma Physics*, 7, 26 (2023).
- B. Afeyan et al., KEEN-wave papers on ponderomotively driven nonlinear electrostatic waves.
- T. W. Johnston, Y. Tyshetskiy, A. Ghizzo, and P. Bertrand, “Persistent subplasma-frequency kinetic electrostatic electron nonlinear waves,” *Physics of Plasmas*, 16 (2009).
- F. Alejandro Padilla-Gomez, Sining Gong, Michael S. Murillo, F. R. Graziani, and Andrew J. Christlieb, “Quantum Kinetic Modeling of KEEN Waves in a Warm-Dense Regime,” *Physics of Plasmas*, 33, 043902 (2026).

---

## 10. What to look for in the numerical results

When reading or reproducing WP simulations, the most important diagnostics are usually

- phase-space slices of $f(x,v,t)$, including sign-changing structures,
- density $\rho(x,t)=\int f\,dv$,
- electric-field energy $\frac12\int |E|^2\,dx$,
- Fourier mode amplitudes of the electric field,
- conservation errors in mass, momentum, and total energy,
- adaptive ranks or hierarchical Tucker ranks,
- imaginary-part diagnostics after inverse velocity-Fourier transforms,
- convergence with respect to grid resolution, time step, and truncation tolerance.

For KEEN-wave studies, harmonic content and wavelet spectra are especially useful.  For high-dimensional adaptive-rank studies, conservation, real-valuedness, rank growth, and empirical scaling are central.

---

## 11. Repository roadmap

Planned additions include

- a method README for the HTACA formulation used in this code,
- scripts for reproducing high-dimensional WP benchmark results,
- result pages for 2D2V and 3D3V two-stream instability and strong Landau damping,
- conservation and Hermitian-symmetry diagnostic notebooks,
- KEEN-wave result pages and plotting tools,
- 1D1V adaptive-rank and conservative low-rank result summaries,
- links to arXiv records as the remaining preprints are posted.

---

## 12. Citation note

If this repository is useful in your work, please cite the relevant paper from the table above.  The high-dimensional HTACA Wigner–Poisson paper is the main reference for the current code base; the 1D1V adaptive-rank paper is the reference for the original ACA/SLAR-Fourier WP formulation; and the KEEN-wave paper is the reference for the warm-dense driven-wave application.
