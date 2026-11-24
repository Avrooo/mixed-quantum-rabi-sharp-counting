---
title: "Supplemental Material for “Exact polynomial eigenstates and sharp counting bounds in the biased mixed quantum Rabi model”"
author: "Sakib Shahrear Avro"
---

# S1. Status of the results

| Statement | Status | Scope |
|---|---|---|
| Termination criterion, energy | Theorem | arbitrary $n$, $\lvert t\rvert<1/2$, $t\neq0$ |
| Limiting roots, Krawtchouk structure | Theorem | arbitrary $n$; spin: spectral identity only |
| Attainment of $B_n$ | Theorem | arbitrary $n$; region depends on $n$, extreme regime |
| Bound $\le B_n$, sector A, $g_2>0$ | Theorem (sharp) | arbitrary $n$, all $\varepsilon$ |
| Bound $\le B_n$, sector A, $g_2<0$ | Theorem | $t\notin\Lambda_n$, all $\varepsilon$ |
| Bound for union of sectors | Corollary | $2B_n$ off $\Lambda_n$; $B_n$ conditional on Problem 2 |
| Common components at $t\in\Lambda_n$ | Open | $n=1$ settled; $n=2$–$4$ numerical, no real $\varepsilon$ found |
| Sector A empty for $g_2<0$ | Open | theorem for $n=1$; exact checks $n\le5$, limit system $n\le8$ |
| Factorization of $\mathrm{Res}(t)$ | Observation | exact for $n\le6$ |

# S2. Inventory of exact checks

All scripts are in `code/`; logs are in `code/logs/`. Checks were written and run by the author; they do not replace independent reading of the proofs.

| Check (script) | Result |
|---|---|
| Closed-form recursion vs. raw operators (`test_exact.py`) | 98/98 exact |
| Four-term recurrence from the top-weight recursion (`audit_appC.py`, part 1) | 1140/1140 residuals vanish exactly |
| Appendix C normalizations as exact identities in $t$ (`audit_appC.py`, part 2) | 12/12 identities hold |
| Product form vs. recursion at $\lambda_+=k/M$ (`audit_appC.py`, part 3, 90 digits) | 87 cases, deviation $2.6\times10^{-72}$; $A_{n+1}=A_{n+2}=0$ |
| Parity, exact stripping, exact degrees, endpoints, $\mathrm{Res}\neq0$ for $t>0$ (`audit_parity_degree.py`) | 120/120 $(n,t)$ pairs, $n\le12$ |
| Limit system: degrees, top-form resultant, $B_n$ rational zeros, distinctness, Jacobians (`audit_limit_bezout.py`) | all pass, $n\le24$ |
| Zero squared-spin eigenvalue; no common zero on $w=0$ (`audit_zero_spin.py`) | pass, $\ell,n\le24$ |
| Krawtchouk identity (`verify_krawtchouk.py`) | $\ell\le12$, deviation $1.5\times10^{-71}$ |
| Spin identity $\det(w(S^2-K^2)-8K^2)$ and sector dimensions, exact rationals (`verify_prop12.py`) | pass, $\ell\le12$; for even $\ell$ the required sector excludes the zero eigenvalue (dimension $\ell/2$) |
| Highest-weight recursion vs. full constraint polynomials (`topw.py`, `logs/topw.log`) | 21/21 configurations ($n=1,\dots,7$; $t=1/10,-1/7,2/5$) |
| Weighted-degree law and leading coefficients (`degree_law.py`) | pass, $n\le14$ |
| Reyes-Bustos--Wakayama consistency at $g_1\to0$ | *no script or log is archived*; stated in the manuscript as a computational observation only, without a case count |
| $n=1$: symbolic $A,B,C,D$ and factorization of $AD-BC$ (`audit_symbolic_n1.py`) | identical; $-6QV$ with all-negative coefficients |
| $n=1$, $t<0$ incompatibility (`audit_n1_negative.py`) | 0 feasible of 2000 rational points |
| $n=1$, $t=t_p$: exact interval enclosure (`audit_tp_certified.py`) | $a_1>0$, $b_1>0$, $R<0$ certified; bracket width $4\times10^{-31}$ |
| $n=1$ threshold $t_0$ (`audit_n1_threshold.py`) | one root in $(0,1/2)$, $t_0=0.308000678$ |
| Factorization of $\mathrm{Res}$ (`audit_divisibility.py`) | $n\le6$: only quartics $P_{k/M}$ and powers of $q$ |
| Exceptional $t$ by Sturm counts (`exceptional.py`) | 1, 2, 4, 5, 8, 10 roots in $(-1/2,0)$; none in $(0,1/2)$ |
| Resultant near $t_*$ (`audit_tstar.py`) | smooth through $t_*$, $n\le6$ |
| Certified root counts (`data/CERTIFIED_ROOT_COUNTS.csv`) | 187 rows $t>0$: resultant degree $=B_n$, count $\le B_n$ (68 attain $B_n$); 111 certified rows $t<0$ empty |
| Fig. 4 vs. certified rows | 24/24 agree |
| Fock diagonalization vs. exact states (`validation_table.py`) | residuals $\le3\times10^{-15}$; overlaps $1.000000000$ |

# S3. Full-precision parameters of Table V

States at $t=0.1$ ($g_2/\omega=0.1923$, $\omega=1$), refined to 40 digits by Newton iteration in 60-digit decimal arithmetic (`code/supp_table5.py`, log `code/logs/supp_table5.log`). The final column is the absolute value of the two constraint residuals.

| $n$ | $\varepsilon$ | $\Delta$ | $g_1$ | $E_n$ | constraint residual |
|---|---|---|---|---|---|
| 1 | 0 | 2.57023503709286602336922650299738072266 | 0.27391653453590756582702933415046067446 | 2.60884504544446907559299490135225005542 | 1e-59 |
| 1 | 1/2 | 2.02884070375945939783305321113958164633 | 0.25141065329047161618796505104011281388 | 2.12805734131382546368137146235128944062 | 1e-59 |
| 2 | 0 | 2.45223910646018944612338574862951743194 | 0.24537864883386598222478051657623439933 | 3.55600379672791155219652768526343618047 | 1e-57 |
| 3 | 0 | 4.48473964915886093413631381708742040955 | 0.21660337366372345240129234336604979480 | 4.50068291701400368751575791118388290823 | 1e-56 |
| 4 | 1/2 | 4.70922282429516957977530534417509838853 | 0.38601728278787159614560616632993191096 | 4.75785981825723610313662696762919339292 | 2e-55 |
| 5 | 1/2 | 5.84037761027254204887702376599819973424 | 0.17854911137128767521962557619949059327 | 5.87127227217327378098297045095501669688 | 1e-54 |
| 6 | 0 | 4.14085368243659083338122782336981285483 | 0.17973615406680958801114008887777501107 | 7.29365808290091329554069843107438572557 | 2e-54 |

# S4. Reproducibility

`code/reproduce.sh` regenerates all logs and figures with exact rational arithmetic unless a precision is stated. Figures: `make_figures.py`, `make_figures2.py`, `make_workflow.py` (900 dpi PNG and vector PDF).

## S5. Fig. 8 checks
States n=1 (eps=0.5), 3, 6 at t=0.1, K=160. Norm 1 to 1e-15; max |P_analytic-P_Fock| <= 7e-16; max density difference <= 1.4e-14; 1-overlap <= 6e-16; |E_Fock-E_n| <= 2e-15; Wigner integral 1 within 1e-10. Overlap deficit < 1e-15 at N_cut=32 (n=1,3) and 48 (n=6).
