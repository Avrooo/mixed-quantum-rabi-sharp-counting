# Proof-audit report (v1.0.0)

## 1. What this audit is and is not
- **Auditor:** an AI assistant (Claude) working with the author. This is an **internal, AI-assisted audit**, not an independent human mathematical referee report. **No claim of proof verification is made.** The proofs remain the author's own and unrefereed (as disclosed in the manuscript Acknowledgments).
- **Method:** (i) re-reading selected proof steps for logical gaps and sign errors; (ii) independent re-derivation of selected elementary facts; (iii) re-running all logged computations from a clean directory and environment; (iv) a new independent numerical check of existing claims using raw operators only (`code/audit_independent.py`). No new theorem, simulation or figure was added.
- **Limits:** computation cannot prove statements for all n. The Word equation objects could not be rendered by Word in the sandbox; see `manuscript/README.md`.

## 2. Errors found and corrected
| # | Location | Problem | Correction |
|---|---|---|---|
| 1 | Prop. 14(ii), case B>0>D | C/D<A/B was said to imply AD-BC>0; it implies AD-BC<0 (since BD<0). Other case was right. Conclusion unaffected (both cases contradict sgn(AD-BC)=sgn Q*). | Text corrected (tests: `tests/test_prop14_sign.py`) |
| 2 | Prop. 4 | Presented as a proposition; equivalence with Reyes-Bustos--Wakayama constraint polynomials not proved for all n and the "112 cases" have no archived script/log | Now Observation 4, computational, no count |
| 3 | Appendix D | Minimal-monomial exactness stated for arbitrary n; verified only n<=12 | Restricted; shown not used by Theorem 6 |
| 4 | Eq. (11)-(12), Sec. VI rescaling | missing overbars (highest-weight parts) | corrected |
| 5 | Appendix C | corrupted `\mathcal\bar A_y` typesetting | corrected |
| 6 | Thm. 13 proof | analyticity across t=0 and openness only asserted | argued term by term (see §3) |
| 7 | Repository | `make_figures.py` failed in a clean checkout (missing `tables/` directory); old `reproduce.sh` ignored script failures | fixed; clean-environment test passes; failures now exit nonzero |

## 3. Results read and what was concluded
Status words: **read** = argument read, no error found beyond §2; **computational** = supported by exact/numerical checks only; none means "proved".

| Result | Audit finding |
|---|---|
| Thm 1, 2 (degree structure, termination) | read. Leading-order bookkeeping (c2, c2', degree alternatives) is consistent; 98/98 exact comparisons with raw operators |
| Prop 3 | read; closed-form residual checked (n=0 in the 98) |
| Obs 4 | computational only, unverified in repo |
| Prop 5 | read; weights/degrees consistent; degrees d1,d0 and B_n=d1*d0=floor((n+1)^2/4) checked algebraically for n=1..12; exactness n<=12 (120/120) |
| Thm 6 | read. Logic = Bezout + Thm 7; needs exact degrees d_i (from Thm 7) but not minimal-monomial exactness; sector-only bound |
| Thm 7 | read. The reality/sign analysis of lambda_+ = 1/2 -/+ (m/2r) sqrt(q/(q-8t)) was re-derived independently: for t<0 lambda_+ in (0,1/2); for 0<t<t*=1-sqrt(3)/2 lambda_+<0; for t*<t<1/2 non-real. Degenerate cases (x=0, g=0, q=8t) are argued in Appendix C.4; checked computationally (276 exact cases n<=12). The Appendix C proof is long and was NOT independently re-derived in full |
| Thm 9 | read. |Lambda_n| = number of reduced k/M in (0,1/2) with 3<=M<=n+2 reproduces 1,2,4,5,8,10,13,15,20,22 (n=1..10); Sturm counts agree for n<=6 |
| Thm 10 | computational for n<=24 (rational roots, Jacobians); the general-n proof via Krawtchouk structure was read but not independently re-derived |
| Prop 11 | computational l<=12 plus proof read |
| Prop 12 | determinant identity exact l<=12 (`verify_prop12.py`); sector statement added |
| Thm 13 | IFT argument now explicit; supported by exact certification (18 configs) and by the new independent check: 21 states (n=2-4, t=1/20,1/50, eps=-100; plus two moderate-bias cases) each with raw-operator residual <=1.7e-12 relative and Fock overlap 1.0 (`audit_independent.py`) |
| Prop 14 | corrected; certified interval arithmetic for t_p |
| Obs 15 | computational n<=6 |
| Cor 16 | follows from Thm 6 and sector exchange; combined bound is qualified in the text |

## 4. Remaining open items
1. No independent human audit; the author should obtain one (Thms 6, 7, 10, 13 are critical).
2. Observation 4: no archived verification; restore with a script or leave as is.
3. Two negative-coupling configurations timed out (uncertified, not claimed).
4. Equations not rendered by Word here; verify the Word-exported PDF yourself.
5. Reference metadata caveats listed in `docs/Supplemental_Material.md` / release notes.
