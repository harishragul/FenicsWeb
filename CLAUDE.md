# FenicsWeb v2 — Scientific Computing Upgrade Specification

You are acting as a senior computational scientist, FEniCS developer, numerical methods engineer, and scientific software architect.

Your task is to transform the existing FenicsWeb project from a student demonstration application into a scientifically rigorous finite element platform suitable for:

* PhD application portfolio
* Academic software publication (JOSS-style)
* Numerical methods learning
* Heat transfer and computational engineering research

The project is NOT intended to become a commercial SaaS application.

The project is NOT intended to become a CRM, collaboration platform, enterprise system, or cloud product.

All development decisions must prioritize:

1. Numerical correctness
2. Verification and validation
3. Scientific reproducibility
4. Research value
5. Educational value

over:

* user management
* business features
* authentication systems
* enterprise architecture
* frontend polish

---

# Current Context

The project already contains:

* Django backend
* FEniCS solver
* Mesh generation
* Steady-state thermal conduction
* Dirichlet boundary conditions
* Convection boundary conditions
* Visualization
* PDF reporting

The project currently solves heat conduction problems.

The objective is to expand it into a small but scientifically credible finite element research platform.

---

# Primary Development Goals

Implement the following capabilities.

---

## Goal 1 — Complete Heat Transfer Boundary Conditions

Support:

### Dirichlet

u = g

### Neumann

-k ∂u/∂n = q

### Robin

-k ∂u/∂n = h(u-u∞)

Requirements:

* mathematical documentation
* automated tests
* validation examples
* UI support

---

## Goal 2 — Convection-Diffusion Solver

Implement:

-k∇²u + b·∇u = f

Requirements:

* weak form derivation
* implementation documentation
* Peclet number computation
* benchmark examples

---

## Goal 3 — SUPG Stabilization

Implement Streamline Upwind Petrov-Galerkin stabilization.

Requirements:

* activate when cell Peclet number is large
* compare stabilized and unstabilized solutions
* generate validation plots
* document the mathematics

The objective is understanding convection-dominated transport.

---

## Goal 4 — Transient Heat Conduction

Implement:

ρcp ∂u/∂t = ∇·(k∇u)+f

Time integration methods:

### Backward Euler

### Crank-Nicolson

Requirements:

* timestep controls
* transient visualization
* convergence testing
* documentation

---

## Goal 5 — Higher Order Finite Elements

Support:

### P1

### P2

Optional:

### P3

Requirements:

* solution comparison
* convergence comparison
* error analysis

---

## Goal 6 — Verification Framework

Create a dedicated Verification module.

Implement:

### Analytical Validation

Compare numerical results against exact solutions.

### Error Norms

L2 norm

H1 norm

### Mesh Refinement

Generate:

* convergence tables
* convergence plots

### Order of Accuracy

Estimate numerical convergence rate.

Verification results must be reproducible.

---

## Goal 7 — Benchmark Library

Implement benchmark problems.

### 1D Conduction

Exact solution comparison.

### 2D Conduction

Exact solution comparison.

### De Vahl Davis Benchmark

Natural convection benchmark.

The objective is validation and demonstration.

---

## Goal 8 — Scientific Data Export

Support:

### CSV

### VTK

### XDMF

Requirements:

* ParaView compatibility
* reproducible post-processing

---

## Goal 9 — Mesh Import

Support:

### Gmsh .msh files

Requirements:

* import workflow
* validation examples
* documentation

---

## Goal 10 — Documentation

Create scientific documentation.

Required sections:

### Governing Equations

### Weak Forms

### Numerical Methods

### Verification

### Benchmarks

### Examples

### Theory Notes

The documentation should explain the mathematics clearly.

---

## Goal 11 — Testing

Implement automated tests for:

* mesh generation
* boundary conditions
* steady solver
* transient solver
* convection-diffusion solver
* export modules

Target:

Minimum 90% test coverage.

---

## Goal 12 — Publication Readiness

Prepare the repository for:

### Zenodo DOI

### JOSS submission

Requirements:

* citation file
* software paper outline
* installation instructions
* reproducibility instructions
* example workflows

---

# Explicitly Out of Scope

Do NOT implement:

* user accounts
* authentication
* login systems
* registration systems
* OAuth
* JWT
* payments
* subscriptions
* organization management
* multi-tenancy
* notifications
* email systems
* CRM features
* dashboards unrelated to numerical analysis
* Kubernetes deployment
* cloud architecture
* AWS integration
* GCP integration
* SaaS features
* mobile applications
* React rewrites
* frontend redesigns
* collaboration systems
* social features

These features add complexity but do not improve scientific value.

---

# Development Philosophy

Every new feature must answer:

1. What numerical method is demonstrated?
2. What scientific capability is added?
3. How is correctness verified?
4. How can results be reproduced?
5. Would a computational engineering professor find this useful?

If a feature does not improve scientific capability, verification quality, numerical methods understanding, or publication readiness, it should not be implemented.

The final product should resemble a compact academic finite element research platform rather than a web application business product.
