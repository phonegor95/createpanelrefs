process ALLELIC_RESCUE_COHORT {
    tag "${meta.id}"
    label 'process_medium'

    container "${params.gatk_container ?: 'https://depot.galaxyproject.org/singularity/gatk4:4.6.1.0--py310hdfd78af_0'}"

    input:
    // meta carries the PoN name; lists are the collected per-sample outputs.
    tuple val(meta), path(rescue_tsvs), path(hetcounts)

    output:
    tuple val(meta), path("cohort_allelic_contrast.tsv"), emit: tsv
    tuple val(meta), path("rescue_contrast/*.private.vcf.gz"),     emit: private_vcf, optional: true
    tuple val(meta), path("rescue_contrast/*.private.vcf.gz.tbi"), emit: private_tbi, optional: true
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Cohort-contrast: a real private deletion is hemizygous in this sample where
    # the rest of the cohort stays heterozygous (high het-support); systematic
    # segdup/ROH artifacts are homozygous in ~everyone (low support); common CNVs
    # recur across the cohort. het-support is the primary discriminant.
    allelic_cohort_contrast.py \\
        --rescue-glob '*.rescue.tsv' \\
        --counts-glob '*.hetcounts.tsv' --counts-prefix '' --counts-suffix '.hetcounts.tsv' \\
        ${args} \\
        --out-tsv cohort_allelic_contrast.tsv --out-dir rescue_contrast

    # bgzip/index the per-sample private-candidate VCFs (if any were written).
    for v in rescue_contrast/*.private.vcf; do
        [ -e "\$v" ] || continue
        bgzip -f "\$v" && tabix -p vcf -f "\$v.gz"
    done
    """

    stub:
    """
    echo -e "sample\tchrom\tstart\tend\tclass\tconf\tcohort_het_support\tcohort_recurrence\tcontrast_class" > cohort_allelic_contrast.tsv
    mkdir -p rescue_contrast
    """
}
