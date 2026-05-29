# Phase 7 — Publication Preparation

**Duration:** 2 weeks
**Status:** Not Started
**Depends on:** All previous phases complete (or substantially complete)
**Unlocks:** JOSS submission, arXiv preprint

---

## Goal

Convert the finished platform into a published, citable software artifact. This phase produces two outputs: an arXiv preprint and a JOSS submission.

By the time Phase 7 starts, the repository will have been public for ~21 weeks — well past the 6-month JOSS minimum.

---

## Tasks

### 7.1 Final CI and test audit

Before submitting anywhere, confirm the test suite is complete and CI is green.

Checklist:

- `pytest tests/ -v` passes locally with zero failures and zero warnings
- GitHub Actions CI badge on README is green
- Test coverage includes at least one test per solver module:
  - `tests/test_conduction.py` — steady-state with Dirichlet, Neumann, Robin BCs
  - `tests/test_convection.py` — Galerkin and SUPG solvers, 1D analytical comparison
  - `tests/test_transient.py` — Backward Euler and Crank-Nicolson, temporal convergence
  - `tests/test_inverse.py` — noiseless recovery error < 0.01
  - `tests/test_mesh.py` — all mesh types create successfully

If any test is missing, write it before proceeding.

---

### 7.2 Zenodo DOI for the software archive

Zenodo provides a permanent DOI for a GitHub repository snapshot — required for the JOSS citation block.

Steps:

1. Go to zenodo.org → Log in with GitHub
2. Enable the FenicsWeb repository in Zenodo's GitHub sync settings
3. Create a GitHub Release (e.g., tag `v1.0.0`) — Zenodo automatically archives it
4. Copy the DOI badge and DOI string (format: `10.5281/zenodo.XXXXXXX`)
5. Add the DOI badge to the README
6. Add the BibTeX citation block to the README:

```bibtex
@software{fenics_web_2026,
  author    = {Your Name},
  title     = {FenicsWeb: A Web-Based Finite Element Framework for
               Convection-Diffusion Heat Transfer},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.XXXXXXX},
  url       = {https://doi.org/10.5281/zenodo.XXXXXXX}
}
```

---

### 7.3 Write the JOSS paper.md

JOSS requires a `paper.md` file at the repository root. Word count: **750–1750 words** (strictly enforced). It is a software description paper — not a research paper. Do not describe results in detail; describe the software and why it matters.

**File:** `paper.md`

Structure (follow this exactly — JOSS editors look for these sections):

```markdown
---
title: 'FenicsWeb: A Web-Based Finite Element Framework for
        Convection-Diffusion Heat Transfer'
tags:
  - Python
  - finite element method
  - heat transfer
  - convection-diffusion
  - FEniCS
  - Django
authors:
  - name: Your Full Name
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
  - name: Your University, Department, Country
    index: 1
date: YYYY-MM-DD
bibliography: paper.bib
---

# Summary
[2–3 sentences: what the software does and what problem it solves]

# Statement of Need
[Why does this tool need to exist? What gap does it fill?
Reference AirTrafficSim, OSSCAR as precedents for web-based
simulation tools. Explain that no web-accessible FEM tool for
convection-diffusion with SUPG stabilisation and inverse problem
capability existed before this work.]

# Features and Functionality
[Describe each solver module in 1–2 sentences each.
Mention BC types, element orders, export formats.]

# Verification and Validation
[1 paragraph: h-convergence study result, De Vahl Davis benchmark.
This is the "research use" evidence JOSS requires.]

# Research Application: Inverse Conductivity Recovery
[1 paragraph on Phase 6 result — the identifiability-vs-Peclet study.
This is the novel result that elevates the paper above a pure software note.]

# Acknowledgements
[Funding, supervisor, institution]

# References
```

The full pre-filled outline is in `docs/publication/paper_outline.md`.

---

### 7.4 Write paper.bib

Create `paper.bib` at the repository root with all references cited in `paper.md`.

Minimum entries required:

```bibtex
@article{fenics_2015,
  author  = {Alnaes, M. S. and others},
  title   = {The FEniCS Project Version 1.5},
  journal = {Archive of Numerical Software},
  year    = {2015},
  volume  = {3}
}

@article{dolfin_adjoint_2013,
  author  = {Farrell, P. E. and others},
  title   = {Automated Derivation of the Adjoint of High-Level
             Mathematical Programs},
  journal = {SIAM Journal on Scientific Computing},
  year    = {2013},
  volume  = {35},
  number  = {4}
}

@article{brooks_hughes_1982,
  author  = {Brooks, A. N. and Hughes, T. J. R.},
  title   = {Streamline Upwind/{Petrov-Galerkin} Formulations for
             Convection Dominated Flows},
  journal = {Computer Methods in Applied Mechanics and Engineering},
  year    = {1982},
  volume  = {32}
}

@article{de_vahl_davis_1983,
  author  = {{De Vahl Davis}, G.},
  title   = {Natural Convection of Air in a Square Cavity:
             A Benchmark Numerical Solution},
  journal = {International Journal for Numerical Methods in Fluids},
  year    = {1983},
  volume  = {3}
}

@article{airtrafficsim_2023,
  author  = {Hui, K. and others},
  title   = {{AirTrafficSim}: An Open-Source Web-Based Air Traffic
             Simulation Platform},
  journal = {Journal of Open Source Software},
  year    = {2023},
  volume  = {8},
  number  = {84}
}
```

---

### 7.5 Post the arXiv preprint

Post before JOSS submission. arXiv gives immediate visibility and does not block JOSS review.

**Target categories:**

- Primary: `cs.CE` (Computational Engineering, Finance, and Science)
- Cross-list: `math.NA` (Numerical Analysis)

**Title:** same as `paper.md` title.

**Abstract** (150–200 words): Summarise the software, the physics covered, the validation results, and the inverse conductivity research finding. End with a sentence on availability.

**Submission steps:**

1. Prepare a PDF of `paper.md` (pandoc or Overleaf works)
2. Submit at arxiv.org — expect 1–2 business days for moderation
3. Once the arXiv ID is assigned, add the arXiv badge to the README

---

### 7.6 JOSS submission

**Pre-submission checklist** (from `docs/publication/joss_checklist.md`):

- [ ] Repo has been public for 6+ months
- [ ] License file is OSI-approved (MIT ✓)
- [ ] `paper.md` and `paper.bib` committed to repo root
- [ ] `CONTRIBUTING.md` exists
- [ ] Tests exist and CI is green
- [ ] Software has a Zenodo DOI
- [ ] README has installation instructions
- [ ] `paper.md` is 750–1750 words

**Submission steps:**

1. Go to joss.theoj.org → Submit
2. Fill in: title, repository URL, software version, language (Python), arXiv preprint URL (optional but recommended)
3. Select editor — suggest "Computational Engineering" scope
4. Pay nothing — JOSS is free to publish
5. An editor will assign two reviewers; typical review time is 4–8 weeks
6. Respond to reviewer comments promptly — most reviews are about documentation and usability, not correctness

**What reviewers check:**

- Does the software install and run from the README instructions?
- Do the tests pass?
- Is the `paper.md` accurate and within word count?
- Does the software provide something not already available?
- Is there evidence of research use?

---

## Files to Create

| File | Purpose |
| ---- | ------- |
| `paper.md` | JOSS submission paper (750–1750 words) |
| `paper.bib` | BibTeX references for paper.md |

**Files to modify:**

| File | Change |
| ---- | ------ |
| `README.md` | Add Zenodo DOI badge, arXiv badge, JOSS badge (after acceptance) |
| `.github/workflows/ci.yml` | Ensure all tests pass on clean checkout |

---

## Definition of Done

- [ ] `pytest tests/ -v` passes with zero failures
- [ ] GitHub Release `v1.0.0` tagged and Zenodo DOI obtained
- [ ] `paper.md` committed, 750–1750 words, all sections present
- [ ] `paper.bib` committed with all cited references
- [ ] arXiv preprint submitted and ID obtained
- [ ] arXiv badge added to README
- [ ] JOSS submission form submitted
- [ ] JOSS submission URL saved (share with thesis supervisor)
