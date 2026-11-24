#!/bin/sh
# Re-runs every verification script and regenerates all figures; any failure stops the run with a nonzero exit status.
# Run from this directory (code/). Needs Python 3 with numpy, scipy, matplotlib; no LaTeX.
set -eu
mkdir -p logs ../figures
for s in test_exact topw degree_law verify_prop12 checks1 checks2 checks3 checks4 symcheck verify_krawtchouk audit_divisibility audit_monotone audit_tstar degdrop exceptional \
         audit_symbolic_n1 audit_n1_negative audit_tp audit_tp_certified audit_n1_threshold audit_parity_degree audit_appC audit_limit_bezout audit_zero_spin audit_independent \
         validation_table_run supp_table5 make_figures make_figures2 make_figures3 make_workflow; do
  echo "== $s"
  if ! python3 -I -X utf8 "$s.py" > "logs/$s.log" 2>&1; then echo "FAILED: $s (see code/logs/$s.log)"; tail -n 8 "logs/$s.log"; exit 1; fi
done
python3 -I check_logs.py
echo "REPRODUCTION COMPLETED"
