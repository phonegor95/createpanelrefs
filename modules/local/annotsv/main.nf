process ANNOTSV {
    tag "${meta.id}"
    label 'process_medium'

    // AnnotSV 3.5.10. Defaults to the biocontainer; override with
    // params.annotsv_container to reuse a site-local .sif and avoid a pull.
    container "${params.annotsv_container ?: 'https://depot.galaxyproject.org/singularity/annotsv:3.5.10--py311hdfd78af_0'}"

    input:
    tuple val(meta), path(vcf)
    path annotations_dir

    output:
    tuple val(meta), path("*.annotsv.tsv"), emit: tsv, optional: true
    tuple val("${task.process}"), val('annotsv'), eval("AnnotSV --version 2>&1 | head -n1 | sed 's/^AnnotSV //'"), topic: versions, emit: versions_annotsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: '-genomeBuild GRCh38'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    AnnotSV \\
        -SVinputFile ${vcf} \\
        -annotationsDir ${annotations_dir} \\
        -outputDir . \\
        -outputFile ${prefix}.annotsv.tsv \\
        ${args}

    # AnnotSV emits nothing when the input has 0 SVs; keep the task green.
    if [ ! -s ${prefix}.annotsv.tsv ]; then
        echo "[ANNOTSV] no SVs annotated for ${prefix}" >&2
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.annotsv.tsv
    """
}
