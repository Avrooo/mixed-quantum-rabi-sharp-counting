# Changelog
## v1.0.0 (2026-11-24)
- First public release. Word manuscript (authoritative), supplement, verification scripts with logs, full-precision data, figures 1-8.
- Fixed Proposition 14 sign (case B>0>D); demoted Proposition 4 to Observation 4; restricted Appendix D claim to n<=12; clarified Prop. 12; made Thm. 13 analyticity/openness explicit.
- Reproduction scripts now stop with a nonzero exit status on any failed calculation; added `check_logs.py` gate, `audit_independent.py`, `verify_prop12.py`, tests; fixed missing-directory failure in `make_figures.py`; LaTeX no longer required.
