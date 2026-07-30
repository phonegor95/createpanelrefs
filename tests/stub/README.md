# Stub-run fixtures

Placeholder inputs for `-profile test_gcnv -stub-run` only.

The optional branches of `GERMLINECNVCALLER_COHORT` (allelic LOH, allelic
rescue, HBA classifier) are selected by `params.allelic_snp_vcf` being non-null,
and the subworkflow resolves it with `checkIfExists: true` along with its `.tbi`.
Under `-stub-run` no process opens either file, so a header-only VCF and an empty
index are enough to exercise the channel wiring — which is the point of the
profile.

These files are deliberately useless for a real analysis. Supply a genuine
common-SNP panel (e.g. gnomAD AF>=0.1, bgzipped and tabix-indexed) via
`--allelic_snp_vcf` for anything other than a stub run.
