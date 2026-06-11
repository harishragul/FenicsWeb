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
  - `tests/test_mesh.py` — all mesh types create successfully
  - `tests/test_export.py`, `tests/test_mesh_import.py` — export and Gmsh import (Phase 5)

If any test is missing, write it before proceeding.

---

### 7.2 CITATION.cff and Zenodo DOI for the software archive

`CITATION.cff` is the GitHub-recognized citation file format (Goal 12) — it
powers the "Cite this repository" button and is what Zenodo reads to
populate the archive's metadata.

```yaml
cff-version: 1.2.0
message: "If you use this software, please cite it as below."
title: "FenicsWeb: A Web-Based Finite Element Framework for Convection-Diffusion Heat Transfer"
authors:
  - family-names: "Your Last Name"
    given-names: "Your First Name"
    orcid: "https://orcid.org/0000-0000-0000-0000"
version: 1.0.0
date-released: "2026-XX-XX"
license: MIT
repository-code: "https://github.com/<user>/FenicsWeb"
```

Zenodo provides a permanent DOI for a GitHub repository snapshot — required for the JOSS citation block.

Steps:

1. Go to zenodo.org → Log in with GitHub
2. Enable the FenicsWeb repository in Zenodo's GitHub sync settings
3. Create a GitHub Release (e.g., tag `v1.0.0`) — Zenodo automatically archives it and reads `CITATION.cff` for metadata
4. Copy the DOI badge and DOI string (format: `10.5281/zenodo.XXXXXXX`)
5. Add the DOI to `CITATION.cff` and the DOI badge to the README
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
simulation tools. Explain that no web-accessible FEM tool combining
convection-diffusion with SUPG stabilisation, an integrated
verification/benchmark suite, and ParaView-compatible export existed
before this work.]

# Features and Functionality
[Describe each solver module in 1–2 sentences each.
Mention BC types, element orders, time integration schemes,
export formats, and Gmsh mesh import.]

# Verification and Validation
[1–2 paragraphs: h- and p-refinement convergence study results
(L2/H1 rates), SUPG vs. Galerkin Peclet study, and the De Vahl Davis
benchmark comparison. This is the "research use" evidence JOSS requires.]

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

**Abstract** (150–200 words): Summarise the software, the physics covered, and the verification/benchmark results. End with a sentence on availability.

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
| `CITATION.cff` | GitHub/Zenodo citation metadata (Goal 12) |
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
- [ ] `CITATION.cff` committed at repository root
- [ ] GitHub Release `v1.0.0` tagged and Zenodo DOI obtained
- [ ] `paper.md` committed, 750–1750 words, all sections present
- [ ] `paper.bib` committed with all cited references
- [ ] arXiv preprint submitted and ID obtained
- [ ] arXiv badge added to README
- [ ] JOSS submission form submitted
- [ ] JOSS submission URL saved (share with thesis supervisor)
