# Mathematical progress: v9, 22 Nov 2026

This report distinguishes proved *conditional consequences* of the manuscripts,
computational checks, and open research problems. In particular, "audited" never
means "ozone group generically determined".

## New partial closure: six two-parameter families

For **M, N, O, P, W, Z**, write `m=ord(u^2)`, `m'=ord(v^2)` as in the
revised monomial draft. On the **even-m stratum** (equivalently
`4 | ord(u)`), its Theorem B, Table `tab:split`, and Propositions
`prop:MNP`, `prop:OWZ` imply

    Oz(B) ≅ C_2 × C_{2m'} × C_{2m'}.

Proof: the realized Galois subgroup is trivial, hence ozone = character
group of Gamma/Gamma_K; the relations of that quotient are
`2m'*alpha=2m'*beta=2*delta=0`. This gives the claimed Smith normal
form. A complete proof in manuscript notation is in `EVEN_M_THEOREM.tex`.
`two_parameter_even_m.py` checks 180 integer parameter/family instances
with exact Smith normal forms. The general argument is mathematical, not
an extrapolation from the instances.

**Dependency:** this does **not** independently prove the monomial
classification, the crossed-product restriction theorem, or the claims of
Propositions MNP/OWZ; these must be verified from defining relations.

## Four starred rows

Under the hypothesis that their tabulated ozone groups are correct and
the general ozone order-divisibility theorem applies, **all four**
PI degrees must be multiples of **6**. This is necessary only; the
ranks, center completeness, and group exhaustivity remain unresolved.
See `STARRED_ROWS_OBSTRUCTIONS.md`.

## Families still open

* C, F, I: generic PI criterion and group require a separate classification;
  source manuscript has only special PI examples.
* D: center constraints insufficient for a full group even at h=1.
* G: generic twist h not done; h=1 proof depends on the center theorem.
* Q, R, V: h=1 cyclic four specialized proofs; generic twists not done.
* M, N, O, P, W, Z: odd m parameter strata still require invariant factors,
  especially non-split swap extensions.
* Starred C*, Q*, R*, V*: center and PI degree/exhaustivity still open.

## Release integrity

- Original v8 scripts and manuscript notes preserved.
- New proofs and tests do not silently change earlier proofs.
- `python -m unittest test_even_m test_starred_constraints test_qrv_h1 test_audit_fast -v`
  passed 11 tests in the current environment.
- All 26 families **remain in the catalog**; several are catalog-only.
- Remote GitHub repository remains unchanged.
