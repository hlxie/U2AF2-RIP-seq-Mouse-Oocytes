# U2AF2 RIP-seq analysis

Code collection for U2AF2 RIP-seq in mouse oocytes and its association with transcriptomic changes after U2af2 deletion. Reference coordinates: **GRCm38/mm10; Ensembl release 101**.

**Release status: organized draft, not a verified reproduction of the final manuscript.** Six historical analysis scripts were recovered. Target collection/plotting scripts were made portable; new upstream utilities are explicitly identified as reconstructed implementations. Real BAMs, peak files, matrices and final run logs were not available during packaging. No manuscript figures or reported statistics were regenerated here.

## Included modules

| File | Purpose | Origin |
|---|---|---|
| `scripts/01_consensus.py` | ≥2/5 IP support segments, merged consensus BED and SAF | Reconstructed |
| `scripts/02_count_peaks.py` | Paired-end fragment counting with featureCounts | Reconstructed from recorded flags |
| `scripts/02_peak_enrichment.py` | Per-library CPM; descriptive IP/Input and IP/IgG enrichment | Reconstructed; explicit pseudocount |
| `scripts/03_intron_profiles.py` | Exon-excluded peak-masked intron and acceptor-upstream profiles | Reconstructed; definitions in METHODS |
| `scripts/15_collect_target_rip.py` | Target-gene peak evidence and actual IP support | Recovered; CLI/path changes |
| `scripts/32_GO2_RIP_gene_expression.py` | Mouse-level aggregation and expression/RIP association preparation | Recovered; portable reference path |
| `scripts/34_RIP_gene_correlations.py` | Twelve prespecified RIP summaries vs expression change | Recovered |
