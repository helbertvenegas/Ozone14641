# Mathematical status and remaining work — v11

All labels below distinguish a **general statement in the manuscript**, a **conditional theorem**, a **particular certified example**, and **open classification**.

| Families | What is available | What is not yet independently certified |
|---|---|---|
| A, H, X, Y | No-PI obstruction reported in source classification | Independent general verification of exclusion |
| B, E, J, K, L, S, T, U | General ozone **order** stated in monomial manuscript; representative relation checks | Independent proof of complete group structure/center action for every PI parameter |
| M, N, O, P, W, Z | General **orders** from manuscript, invariant factors proposed for even and odd m; exact finite-phase/Smith checks | Primary relation-by-relation verification of crossed-product action tables and uniform center-fixation/exhaustivity |
| D | Examples and necessary conditions for h=1 | Generic ozone classification including h twists; missing central restrictions |
| G | h=1 formula conditional on claimed center theorem | Verify center theorem and generic twists h |
| Q, R, V | h=1 cases described as C4 with normal-element / central upper-bound certificates (subject to source presentations) | Generic PI twists and starred cases |
| C, F, I | Particular PI examples in old source | General PI criteria and ozone classification |

## Four unresolved starred examples

`C* (p=ω, h²=-1/3)`; `Q*, R*, V* (h=ω)` have reported groups but **no certified PI degree/central rank** in the historical table. If those claimed orders (12 or 36) are correct, the general divisibility bound forces 6 to divide their PI degrees; this is necessary, not sufficient.

## Odd-m caveat

`ODD_M_THEOREM.tex` claims classification **conditional on** the manuscript's Theorem B and action/restriction tables. `test_odd_m.py` checks 720 formal root-of-unity order combinations for the finite-group phase computation. Some are excluded by the families' permitted parameter locus; therefore 720 must not be described as 720 valid PI algebras or 720 independent center checks.

## Release readiness

**Ready to upload:** an honestly labelled reproducible *work-in-progress supplementary repository*.

**Not ready to claim:** a complete referee-grade certificate of the 22 families with PI members. Before such a claim: compare all 26 Sigma matrices with ZZ09, complete the 8 open generic groups, check all four starred centers/PI degrees, audit the odd-m action tables, and run full CI/notebooks against pinned environment.
