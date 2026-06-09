process CNV_QC_OUTLIER {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/htslib_pysam_tabix_pip_variant-extractor:a12ef217eccf6ba8'
        : 'community.wave.seqera.io/library/htslib_pysam_tabix_pip_variant-extractor:a12ef217eccf6ba8'}"

    input:
    tuple val(meta), path(vcfs)
    path segdup_bed
    path blacklist_bed

    output:
    tuple val(meta), path("cnv_qc_per_sample.tsv"), emit: report
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    mad_k = params.qc_outlier_mad_k ?: 3
    template 'cnv_qc_outlier.py'

    stub:
    """
    echo "sample\ttotal_calls\tDEL\tDUP\tCN0\tLCR_calls\tLCR_frac\tCNV_Mb\tflag\tflag_reason" > cnv_qc_per_sample.tsv
    echo '"${task.process}":' > versions.yml
    echo '    python: 3.12.0' >> versions.yml
    """
}
