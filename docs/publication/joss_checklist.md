# JOSS Submission Checklist

**Reference:** https://joss.readthedocs.io/en/latest/submitting.html

This checklist must be fully green before submitting to JOSS. Each item maps to an explicit JOSS requirement or a known reviewer concern based on prior FEniCS tool submissions.

---

## Repository Requirements

- [ ] **Repository is public** on GitHub
  - Phase 0 action. Must have been public for 6+ months before submission.

- [ ] **6+ months of genuine commit history**
  - "Genuine" means regular commits throughout development, not a burst before submission.
  - Run `git log --format="%ad" --date=short | tail -1` to check earliest commit date.

- [ ] **OSI-approved open source license**
  - MIT license already in place ✓
  - Confirm `LICENSE` file is at the repo root.

- [ ] **Repository is not archived or read-only**

---

## Software Requirements

- [ ] **Installable from the README instructions**
  - Test on a clean machine or fresh conda environment.
  - If installation fails, a reviewer will desk-reject.

- [ ] **Tests exist and pass**
  - `pytest tests/ -v` must exit with code 0.
  - Tests must be in the repo, not just described.

- [ ] **CI is configured and passing**
  - GitHub Actions badge must be green on the README.
  - CI must run `pytest` on every push to `main`.

- [ ] **Codebase is substantial**
  - JOSS auto-rejects repos under ~300 lines. FenicsWeb will be well above this after Phase 1.
  - Run `find . -name "*.py" | xargs wc -l | tail -1` to check total line count.

- [ ] **Software version is tagged**
  - Create a git tag `v1.0.0` before submission.
  - This tag is what Zenodo archives.

---

## Documentation Requirements

- [ ] **README with installation instructions**
  - Must include: what the software does, how to install it, how to run it.
  - Screenshots strongly encouraged.

- [ ] **`CONTRIBUTING.md` exists**
  - Must explain how to set up a development environment and run tests.

- [ ] **Community guidelines present**
  - Either in CONTRIBUTING.md or a separate `CODE_OF_CONDUCT.md`.

- [ ] **API documentation or docstrings**
  - Every public function must have at least a one-line docstring.

---

## Paper Requirements

- [ ] **`paper.md` at repository root**
  - Must include YAML front-matter: title, tags, authors with ORCID, affiliations, date, bibliography.

- [ ] **Word count: 750–1750 words**
  - Run `pandoc paper.md | wc -w` to check. JOSS editors enforce this strictly.

- [ ] **All required sections present:**
  - Summary
  - Statement of Need
  - (Functionality description — can be merged into Summary or separate)
  - References

- [ ] **`paper.bib` at repository root**
  - All `[@key]` citations in `paper.md` must have a matching entry in `paper.bib`.

- [ ] **Authors have ORCID iDs**
  - Register free at orcid.org if not already done.

- [ ] **References are complete and correct**
  - FEniCS, dolfin-adjoint, SUPG (Brooks & Hughes), and any benchmark paper cited must have full bibliographic entries.

---

## Research Use Requirement

JOSS explicitly rejects tools described only in aspirational terms ("could be used for research"). At least one of the following must be true:

- [ ] **The authors have used it in a published paper or thesis** — your master's thesis counts.
- [ ] **The paper.md describes a concrete research result** — the Phase 6 inverse problem study satisfies this.
- [ ] **Other researchers have cited or used the tool** — optional at submission time; will grow after publication.

The Phase 6 inverse conductivity study (identifiability vs Peclet number) is the primary evidence of research use. Reference it explicitly in the paper's Statement of Need and in a dedicated section.

---

## Zenodo DOI Requirement

- [ ] **Software archived on Zenodo**
  - Zenodo DOI in the README and in the `paper.md` YAML front-matter.
  - Format: `10.5281/zenodo.XXXXXXX`

---

## Known Reviewer Concerns for FEniCS Tools

Based on published reviews of similar tools (`pulse`, `FEniCS-arclength`, `fenicsx-beat`):

1. **Reviewer will try to install and run it.** Make sure Docker and conda paths both work.
2. **Reviewer will run the tests.** Flaky tests that pass locally but fail in CI will block acceptance.
3. **Reviewer will check the Statement of Need carefully.** Be specific: what problem does FenicsWeb solve that existing tools do not? Answer: web accessibility, integrated V&V dashboard, and the inverse problem capability in a convection-diffusion setting.
4. **Reviewer may ask for more examples.** Prepare a Jupyter notebook or a worked example script as a supplement.

---

## Post-Acceptance Actions

After JOSS acceptance:

- [ ] Add JOSS badge to README: `[![DOI](https://joss.theoj.org/papers/XXXX/status.svg)](https://doi.org/10.21105/joss.XXXX)`
- [ ] Update Zenodo record with the JOSS DOI
- [ ] Add JOSS citation to the thesis bibliography
- [ ] Announce on the FEniCS community forum (fenicsproject.discourse.group)
