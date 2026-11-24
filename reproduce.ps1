# Windows reproduction (PowerShell). Stops with a nonzero exit code on the first failure. No LaTeX required.
$ErrorActionPreference = "Stop"
Set-Location "$PSScriptRoot\code"
New-Item -ItemType Directory -Force -Path logs | Out-Null
$scripts = "test_exact","topw","degree_law","verify_prop12","checks1","checks2","checks3","checks4","symcheck","verify_krawtchouk","audit_divisibility","audit_monotone","audit_tstar","degdrop","exceptional","audit_symbolic_n1","audit_n1_negative","audit_tp","audit_tp_certified","audit_n1_threshold","audit_parity_degree","audit_appC","audit_limit_bezout","audit_zero_spin","audit_independent","validation_table_run","supp_table5","make_figures","make_figures2","make_figures3","make_workflow"
foreach ($s in $scripts) {
  Write-Host "== $s"
  python -I -X utf8 "$s.py" *> "logs\$s.log"
  if ($LASTEXITCODE -ne 0) { Write-Host "FAILED: $s"; Get-Content "logs\$s.log" -Tail 8; exit 1 }
}
python -I check_logs.py; if ($LASTEXITCODE -ne 0) { exit 1 }
Set-Location ..
python -I -m pytest -q -p no:cacheprovider tests; if ($LASTEXITCODE -ne 0) { exit 1 }
Write-Host "ALL STEPS PASSED"
