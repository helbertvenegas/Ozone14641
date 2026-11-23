# Verification status, November 22 2026

| Family group | Scope | Mathematically verified by this release | Pending |
|---|---|---|---|
| B E J K L M N O P S T U W Z | 14 monomial families | Exact parameter exclusions at one representative specialization each; plane relations under source basis changes; PBW reduction and 20 automorphism-pattern equivalence tests per family; numerical formulas in prior notebook | Independence from ZZ09 matrix source; pointwise center action; full ozone-group exhaustivity for arbitrary PI parameters |
| D G Q V | Nonmonomial families with PI criteria in working draft | Existing G PI-degree spot check from v2 | Source matrices Sigma, center generators, group matrices and exhaustive certificates |
| C F I R | Remaining nonmonomial PI examples/criteria | Table arithmetic consistency only | Generic PI criterion as applicable; Sigma, exact center and ozone calculations |
| A H X Y | No PI members in working classification | Recorded as exclusions only | Cite exact theorem establishing each obstruction |

**Crucial**: The original Sigma matrices in `source_monomial.py` are derived by conjugating the monomial relations in `tab:mono`. This is an algebraic *internal consistency check* and must not be mistaken for independent verification against Section 4 of Zhang–Zhang.

The older manuscript's `cor:ozsist` requires: all six relations, pointwise fixation of a proven generating set of `Z(B)`, and nonvanishing determinant. Finite-field Groebner counts do not automatically certify characteristic-zero completeness.

The table's four starred rows (`C*`, `Q*`, `R*`, `V*`) remain conditional; no PI-rank or exhaustive ozone certificate is asserted for them.

The uniform ozone theorem of the more recent monomial manuscript provides **theorem-dependent** formulas for 14 families; software verification of selected specializations does not independently establish those theorems.

## v4: D/G mathematical upgrade

- **G(p,f,1), n=ord(p)>=3, f admissible**: conditional on the full-center theorem `thm:G` in the supplied manuscript, a complete symbolic **proof** now gives the invariant factors of Oz for every n. See `DG_THEOREM_NOTE.tex`. `dg_exact.py` checks order spectra and abstract group structures for n=3..20 using modular integer arithmetic. The proof, not numerical enumeration, establishes the general theorem.
- **D(p,1)**: a rigorous necessary-condition reduction for every n>=3 yields diagonal candidate matrices. These bounds are **not** a full ozone determination (n=4: 32 candidates versus table order 16); a missing full-center generator is required.
- **D/G with arbitrary h**: no general ozone theorem has been verified in this release. Center computations at h=1 cannot simply be transferred to h≠1.

## v5 decision: stop the D/G closure attempt and proceed to Q/R/V

- `D(p,1)` is **not closed**: for primitive fourth root `p`, the known central constraints leave 32 possible matrices versus the historical table's 16 ozone automorphisms. The missing center constraints are not guessed.
- `G(p,f,1)` has a **conditional general proof**, given the full-center theorem in the working manuscript; `G(p,f,h)` for general root-of-unity twist `h` is **not closed**.
- `Q` and `V`: `qrv_dynamics.py` certifies the characteristic polynomials of their normal-localization matrices, preservation of their affine conics, and exact order four of the Vieta maps, all symbolically over the rationals. These are **not ozone-group certificates**.
- `R`: PI criterion attributed to ZZ08 in the manuscript; a degree-two normal element requires a distinct localization model. The prior tabulated ozone examples are not yet certified in characteristic zero.
- The original four `Q*/R*/V*/C*` rows still carry a rank/center completeness caveat.

## v6: Q(1) theorem; V(1) exact upper bound

- **Q(h=1): PROVED**, conditional only on the sourced Zhang--Zhang/Gómez--Venegas defining presentation: `Oz(Q(1)) = C4`. An exact center computation produces a new central element of bidegree (1,2); this rules out half the eight degree-one candidates. `Ad_x1` is a normal-inner ozone generator of order four, proving exhaustivity. See `QRV_H1_THEOREMS.tex` and `qrv_ozone_h1.py`.
- **V(h=1): upper bound only**, `|Oz(V(1))| <= 4` from defining relations and known central elements. All four candidates preserve central components in tested low bidegrees but this does **not** prove center-pointwise action for arbitrary degree.
- **Q,V with arbitrary root-of-unity h**: still open.
- **R**: still requires exact source relations and independent degree-two normal model.
- Tests: symbolic rational PBW reductions, not finite-field calculations.

## v7: three *specialized* proofs at h = 1

- **Q(1)**: exact group C4 via central element of bidegree (1,2) and normal x1 (`qrv_ozone_h1.py`).
- **R(1)**: exact group C4 via central element z=x1*y1+x1*y2+x2*y1-x2*y2 and nonzero normal u=x2*y1+x1*y2 (`r1_normal_certificate.py`). R defining relations independently transcribed from https://4d-as.ncag.info/families/ore-r/, pending final check with the original ZZ09.
- **V(1)**: exact group C4 via normal w=x1*(y1^2-y2^2)-x2*(y1^2+y1*y2), whose conjugation has order four, meeting a previous upper bound four (`v1_normal_certificate.py`).
- **Q, R, V generic PI twist h**: **not closed**. The table's h=-1, h=i and starred examples are not certified by these three h=1 theorems.
- Focused tests: `python -m unittest test_qrv_h1 -v`; CI workflow `.github/workflows/qrv-exact.yml` is provided but has not been executed on GitHub.
- See `INFORME_AVANCE.md` for a full article-and-repository status assessment.
