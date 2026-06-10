process BCFTOOLS_FILTER_GCNV {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0'
        : 'biocontainers/bcftools:1.20--h8b25389_0'}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.filtered_segments.vcf.gz"),     emit: vcf
    tuple val(meta), path("*.filtered_segments.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | head -n1 | sed 's/^bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    // The sex-aware include expression is supplied via ext.args (see
    // conf/modules/germlinecnvcaller_cohort.config), selected from meta.sex
    // which is derived from the per-sample chrX ploidy call.
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools filter \\
        ${args} \\
        --output-type z \\
        --output ${prefix}.filtered_segments.vcf.gz \\
        ${vcf}

    bcftools index --tbi ${prefix}.filtered_segments.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip -c > ${prefix}.filtered_segments.vcf.gz
    touch ${prefix}.filtered_segments.vcf.gz.tbi
    """
}
