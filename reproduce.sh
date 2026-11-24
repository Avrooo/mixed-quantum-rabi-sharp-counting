#!/bin/sh
# Full reproduction: verification scripts, figures, log gate, tests. Any failure => nonzero exit status. No LaTeX required.
set -eu
cd "$(dirname "$0")/code"
sh ./reproduce.sh
cd ..
python3 -I -m pytest -q -p no:cacheprovider tests
echo "ALL STEPS PASSED"
