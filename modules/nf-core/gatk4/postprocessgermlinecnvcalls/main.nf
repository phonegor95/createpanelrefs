process GATK4_POSTPROCESSGERMLINECNVCALLS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity'
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta), path(model_shards), path(call_shards), path(ploidy_calls), path(dict), val(sample_index)

    output:
    tuple val(meta), path("*.genotyped_intervals.vcf.gz"), emit: genotyped_intervals
    tuple val(meta), path("*.genotyped_intervals.vcf.gz.tbi"), emit: genotyped_intervals_tbi, optional: true
    tuple val(meta), path("*.genotyped_segments.vcf.gz"), emit: genotyped_segments
    tuple val(meta), path("*.genotyped_segments.vcf.gz.tbi"), emit: genotyped_segments_tbi, optional: true
    tuple val(meta), path("*.denoised_copy_ratios.tsv"), emit: denoised_copy_ratios
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def model_shard_args = model_shards.collect { "--model-shard-path ${it}" }.join(' ')
    def call_shard_args = call_shards.collect { "--calls-shard-path ${it}" }.join(' ')

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK PostprocessGermlineCNVCalls] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    export THEANO_FLAGS="base_compiledir=\$PWD"
    export PYTENSOR_FLAGS="base_compiledir=\$PWD"
    export OMP_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}

    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        PostprocessGermlineCNVCalls \\
        ${model_shard_args} \\
        ${call_shard_args} \\
        --contig-ploidy-calls ${ploidy_calls} \\
        --sample-index ${sample_index} \\
        --sequence-dictionary ${dict} \\
        --output-genotyped-intervals ${prefix}.genotyped_intervals.vcf.gz \\
        --output-genotyped-segments ${prefix}.genotyped_segments.vcf.gz \\
        --output-denoised-copy-ratios ${prefix}.denoised_copy_ratios.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.genotyped_intervals.vcf.gz
    touch ${prefix}.genotyped_segments.vcf.gz
    touch ${prefix}.denoised_copy_ratios.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
