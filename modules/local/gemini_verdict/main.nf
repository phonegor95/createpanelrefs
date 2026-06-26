process GEMINI_VERDICT {
    tag "${meta.id}"
    label 'process_single'

    // python3 (stdlib only) is enough; reuse the GATK image already pulled.
    container "${params.gatk_container ?: 'https://depot.galaxyproject.org/singularity/gatk4:4.6.1.0--py310hdfd78af_0'}"
    // GEMINI_API_KEY is injected from Nextflow secrets (never written to .command.sh):
    //   nextflow secrets set GEMINI_API_KEY "<key>"
    secret 'GEMINI_API_KEY'

    input:
    tuple val(meta), path(yaml)

    output:
    tuple val(meta), path("*.HBA.verdict.md"), emit: verdict

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def model  = params.gemini_model ?: 'gemini-3.5-flash'
    def proxy  = params.gemini_proxy ? "--proxy ${params.gemini_proxy}" : ''
    """
    gemini_verdict.py \\
        --sample ${prefix} \\
        --yaml ${yaml} \\
        --model ${model} ${proxy} \\
        --out ${prefix}.HBA.verdict.md
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '# HBA independent verdict  sample=%s\\nstatus: OK\\nVERDICT: AGREE\\nCONFIDENCE: high\\n' ${prefix} > ${prefix}.HBA.verdict.md
    """
}
