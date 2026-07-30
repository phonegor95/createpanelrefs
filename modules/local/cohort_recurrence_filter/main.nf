process COHORT_RECURRENCE_FILTER {
    tag "${meta.id}"
    label 'process_single'

    // python3 + bgzip/tabix in one image (htslib + pysam mulled biocontainer)
    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/htslib_pysam_tabix_pip_variant-extractor:a12ef217eccf6ba8'
        : 'community.wave.seqera.io/library/htslib_pysam_tabix_pip_variant-extractor:a12ef217eccf6ba8'}"

    input:
    tuple val(meta), path(sex_map), path(vcfs)

    output:
    tuple val(meta), path("recurfilt/*.recurfilt.vcf.gz"),          emit: flagged
    tuple val(meta), path("recurfilt/*.recurfilt.vcf.gz.tbi"),      emit: flagged_tbi
    tuple val(meta), path("recurfilt/*.recurfilt.pass.vcf.gz"),     emit: pass
    tuple val(meta), path("recurfilt/*.recurfilt.pass.vcf.gz.tbi"), emit: pass_tbi
    // python3 template ⇒ eval() outputs are not permitted (Bash-only); container is pinned by hash so versions are static.
    tuple val("${task.process}"), val('python'), val('3.12.4'), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('htslib'), val('1.20'),   topic: versions, emit: versions_htslib

    when:
    task.ext.when == null || task.ext.when

    script:
    // Two-arm cohort-recurrence thresholds (overridable via params).
    recip_threshold  = params.recur_recip_threshold  ?: 0.5
    recip_coverage   = params.recur_recip_coverage   ?: 0.5
    overlap_threshold = params.recur_overlap_threshold ?: 0.8
    // Below this many peers per call the frequency is too quantised to act on;
    // the filter degrades to a no-op instead of flagging near-arbitrarily.
    min_cohort       = params.recur_min_cohort ?: 5
    template 'cohort_recurrence_filter.py'

    stub:
    """
    mkdir -p recurfilt
    for v in ${vcfs}; do
        s=\$(basename \$v .filtered_segments.vcf.gz)
        echo | gzip -c > recurfilt/\${s}.recurfilt.vcf.gz
        touch recurfilt/\${s}.recurfilt.vcf.gz.tbi
        echo | gzip -c > recurfilt/\${s}.recurfilt.pass.vcf.gz
        touch recurfilt/\${s}.recurfilt.pass.vcf.gz.tbi
    done
    """
}
