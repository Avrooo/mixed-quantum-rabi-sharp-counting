# v1.0.0 release notes
**Tag:** `v1.0.0`. **Title:** Exact polynomial eigenstates and sharp counting bounds in the biased mixed quantum Rabi model - code, data and logs.

Steps for GitHub + Zenodo (to be done by the author; not done in the build sandbox, which has no network):
1. Create the public GitHub repository `mixed-quantum-rabi-sharp-counting`, push this tree, replace `<OWNER>` in `CITATION.cff` and `.zenodo.json`.
2. Enable the repository in Zenodo (Account > GitHub), then create the GitHub release `v1.0.0` from this tree.
3. After Zenodo issues the DOI, add it to `CITATION.cff` (`doi:`), the README, and the manuscript data-availability text; this requires a v1.0.1 metadata-only release or a Zenodo version update.
4. Verify `SHA256SUMS.txt` after cloning: `sha256sum -c SHA256SUMS.txt`.

Known limitations: see `docs/PROOF_AUDIT_REPORT.md` (no independent human audit; Observation 4 unverified; two uncertified timeout rows).
