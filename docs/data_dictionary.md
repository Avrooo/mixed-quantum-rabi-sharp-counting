# Data dictionary
- `CERTIFIED_ROOT_COUNTS.csv`: n, t (exact rational), g2_over_omega, eps (exact), bezout (B_n), degP/degQ (constraint degrees), res_*_deg and *_sqfree (resultant degrees), pos_v_roots/pos_u_roots (positive roots), count (certified number of admissible (Delta^2,g1^2)), pairing_injective, status (certified, certified_empty_by_Sturm, TIMEOUT_240s), seconds.
- `COLLAPSE_COUNTS.csv`, `EXACT_EXCEPTIONAL_CASES.csv`, `EXCEPTIONAL_EPS_SCAN.csv`, `LIMIT_SYSTEM_COUNTS.csv`: outputs of `exceptional.py`, `exc_scan.py`, `limit_sys.py`.
- `fig8/fig8_summary.csv`: per state n: eps, Delta, g1, E_exact, N_cut (160), norm, weight_up/down, mean photon number, Wigner min and negative volume (integral of max(-W,0)), Wigner integral, analytic-vs-Fock maxima, 1-overlap, |dE|, truncated weight (last 20 Fock levels).
- `fig8/fig8_cutoff_convergence.csv`: 1-overlap and |dE| vs N_cut (8..96).
- `fig8/fig8_state_n*.npz`: `fock_coeffs_analytic`, `fock_coeffs_fockdiag` (length 2*160: spin-up block then spin-down block), `P_analytic`, `P_fockdiag`, `x`, `rho_up`, `rho_down`, `p_grid`, `W_xp` (rows x, columns p).
