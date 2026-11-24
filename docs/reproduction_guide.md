# Reproduction guide
1. `pip install -r requirements.txt` (Python 3.13).
2. `sh reproduce.sh` regenerates every log in `code/logs/`, all figures in `figures/` and `code/data/fig8/`, then runs `pytest`.
3. Compare logs with the claims in `docs/THEOREM_STATUS.csv`. Counts: 98 (`test_exact.log`), 21 (`topw.log`), 1140 (`audit_appC.log`), 120/120 (`audit_parity_degree.log`), 18 (`checks3.log`), 111 certified-empty negative rows (`CERTIFIED_ROOT_COUNTS.csv`), 2000 (`audit_n1_negative.log`).
4. Fig. 8: `code/data/fig8/` contains every number plotted; no manual extraction from the manuscript is required.
Run time and numerical results can differ in the last digits across BLAS builds; checks use tolerances (1e-12) except exact rational arithmetic.
