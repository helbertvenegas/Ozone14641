# Ozone groups of double Ore extensions of type (14641)

**Authors:** James Gómez and Helbert Venegas  
**Release:** v11 — research verification package, 22 November 2026  
**Status:** **partial mathematical certification**; NOT a complete classification of all PI members.

This repository organizes computational checks and mathematical notes for the 26 Zhang–Zhang families A–Z. The catalogue is comprehensive; proofs of generic ozone groups are not.

## Start here

1. Open [`RELEASE_STATUS.md`](RELEASE_STATUS.md) for verified scope and remaining tasks.
2. Open [`ozone_groups_verification.ipynb`](ozone_groups_verification.ipynb) in Jupyter or Colab for the historical table audit and quadratic-relation examples. For the new invariant-factor work, run [`odd_m_verification.ipynb`](odd_m_verification.ipynb) and [`even_m_and_starred_checks.ipynb`](even_m_and_starred_checks.ipynb). These notebooks use sibling `.py` and `.csv` files; upload the complete repository to Colab or run from the repository root.
3. Read [`ODD_M_THEOREM.tex`](ODD_M_THEOREM.tex) and [`EVEN_M_THEOREM.tex`](EVEN_M_THEOREM.tex). These are **conditional consequences** of the cited manuscript tables and theorems, not independent reconstructions of all defining relations.

## Reproduce the tested checks

Requires Python 3.11+ and packages in `requirements.txt`. From repository root:

```bash
python -m pip install -r requirements.txt
python -m unittest test_even_m test_starred_constraints test_qrv_h1 test_audit_fast test_odd_m -v
python audit_all_families.py
python run_certificates.py
```

The unit-test suite above is the **focused quick suite**; `python -m unittest discover -v` also includes longer symbolic tests and is **not certified as passing in full**. GitHub Actions executes the focused tests; a local PASS does not assert remote GitHub success.

## What PASS means

`AUDIT_26_FAMILIES.md` uses PASS for the narrow test assigned to a family, possibly only a catalogue check or historical example. It does **not** mean that `Oz(B)` is classified for all PI parameters. `ozone_table.csv` contains 73 historic examples whose arithmetic consistency was checked; center-fixation and group exhaustivity were not checked in each case.

## Contents

- `family_catalog.csv`, `ozone_table.csv`: family index and historic examples
- `general_quadratic.py`, `degree2_certificates.py`, `source_monomial.py`: degree-two quadratic/PBW and monomial-coordinate checks
- `odd_m_analysis.py`, `two_parameter_even_m.py`: exact finite-group invariant-factor calculations
- `qrv_ozone_h1.py`, `r1_normal_certificate.py`, `v1_normal_certificate.py`: focused h=1 examples Q/R/V
- `audit_all_families.py`: individual audit report generator
- `RELEASE_STATUS.md`, `MANIFEST_SHA256.txt`: scope and file integrity

## Provenance

Mathematical statements and presentation data are based on the Gómez–Venegas draft manuscripts provided for this project and on the Zhang–Zhang classification they cite. The program contains derived family presentations, not a separately audited transcription of all 26 original Sigma matrices. Cross-check primary references before citing this as supplementary material in a journal.

## Upload to GitHub

Extract the archive, then upload the **contents of this folder** to the root of the chosen repository (not a zip file alone if the goal is a browsable Jupyter notebook). `.github/workflows` should be preserved. This package does not modify GitHub.
