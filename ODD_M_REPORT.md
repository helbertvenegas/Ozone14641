# Odd-m ozone groups: conditional classification and exact tests

Authors: James Gómez and Helbert Venegas. Working certificate, November 2026.

**Scope.** Six PI monomial families M, N, O, P, W, Z. This result is conditional on the monomial presentations, restriction criterion and crossed-product action tables in the revised Gómez–Venegas manuscript (`Texto pegado (3).txt`). It independently verifies the finite-group extension arithmetic, not the published status of those foundational lemmas.

Put `m=ord(u^2)` (odd) and `d=ord(v^2)`. Define

- `A = C_2 × C_(2d) × C_(2d)`; order `8d²` (split index-2 extension).
- `B = C_(4d) × C_(2d)`; order `8d²` (non-split index-2 extension).
- `D = C_(4d) × C_(4d)`; order `16d²`.
- `E = C_(2d) × C_(2d)`; order `4d²`.

| Family | ord(u) odd, ord(v) odd | ord(u) odd, ord(v) even | ord(u) ≡2 mod4, ord(v) odd | ord(u) ≡2 mod4, ord(v) even |
|---|---|---|---|---|
| M, P | A | B | B | A when ord(v) ≡2 mod4; otherwise B |
| N | B | A when ord(v) ≡2 mod4; otherwise B | A | B |
| O, W, Z | E | D | A | B |

**Methods.** Rational phases represent roots of unity exactly. For each specialization the script reads the manuscript's action on `r,rho,Xtilde,Ytilde`, compares automorphisms on `K`, computes squares of lift matrices, and reduces their residues in `N_K/2N_K`. The `|G|=2` case requires the field `L=K(c)` and the extra degree relation, addressed mathematically in the theorem note.

**Tests.** `python test_odd_m.py` checks 720 phase instances, with U=ord(u) and V=ord(v). This is **not** an exhaustive proof for all roots of unity; the uniform symbolic argument is in `ODD_M_THEOREM.tex`.

**Known limitations.** Not a revalidation of the full centralizer/Galois argument, the original matrices of Zhang–Zhang, the PI classification for all families, or all earlier GitHub notebook cells. Do not treat this certificate as resolving C, D, F, G, I, Q, R, V.
