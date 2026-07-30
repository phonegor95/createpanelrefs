process ALLELIC_RESCUE {
    tag "${meta.id}"
    label 'process_medium'

    // GATK (CollectAllelicCounts) + python3 (stdlib only) in one image.
    container "${params.gatk_container ?: (workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/gatk4:4.6.1.0--py310hdfd78af_0'
        : 'biocontainers/gatk4:4.6.1.0--py310hdfd78af_0')}"

    input:
    tuple val(meta), path(cram), path(crai), path(gcnv_vcf)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)
    tuple path(snp_vcf), path(snp_tbi)
    path segdup_bed

    output:
    tuple val(meta), path("*.rescue.tsv"),        emit: tsv
    tuple val(meta), path("*.rescue.vcf.gz"),     emit: vcf
    tuple val(meta), path("*.rescue.vcf.gz.tbi"), emit: tbi
    tuple val(meta), path("*.hetcounts.tsv"),     emit: hetcounts
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sd = segdup_bed ? "--segdup-bed ${segdup_bed}" : ''
    def gv = gcnv_vcf ? "--gcnv-vcf ${gcnv_vcf}" : ''
    def args = task.ext.args ?: ''
    def avail_mem = (task.memory ? (task.memory.giga * 0.85).intValue() : 6)
    """
    # 1. Genome-wide per-SNP allelic counts at high MAPQ on the (recal) CRAM. The
    #    high min-MAPQ is what makes the hemizygosity readable at uniquely-mappable
    #    common SNPs; the SNP panel should be common biallelic SNVs (e.g. gnomAD
    #    AF>=0.1) so the het background is dense.
    gatk --java-options "-Xmx${avail_mem}g" CollectAllelicCounts \\
        -I ${cram} -R ${fasta} -L ${snp_vcf} \\
        --minimum-mapping-quality 30 -O ${prefix}.allelic_counts.tsv

    # 2. Route-A scorer: focal hemizygous-LOH runs in a het background, classed by
    #    local depth ratio + segdup context (annotate/flag, never an auto-call).
    allelic_seg_call.py \\
        --allelic-counts ${prefix}.allelic_counts.tsv --sample ${prefix} \\
        ${sd} ${gv} ${args} \\
        --out-tsv ${prefix}.rescue.tsv --out-vcf ${prefix}.rescue.vcf

    bgzip -f ${prefix}.rescue.vcf && tabix -p vcf -f ${prefix}.rescue.vcf.gz

    # 3. Compact per-sample het-site counts (chrom pos ref alt) for the cohort
    #    contrast step -- a small fraction of the genome-wide counts.
    awk 'BEGIN{OFS="\t"} !/^@/ && \$1!="CONTIG"{ dp=\$3+\$4; if(dp>=8){ b=\$4/dp;
        if(b>=0.25 && b<=0.75) print \$1,\$2,\$3,\$4 } }' \\
        ${prefix}.allelic_counts.tsv > ${prefix}.hetcounts.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "sample\tchrom\tstart\tend\tspan_bp\thom_sites\tmean_depth_ratio\tflank_het\tflank_het_rate\tloh_p\tgcnv_CN\tclass\tconfidence" > ${prefix}.rescue.tsv
    echo | gzip -c > ${prefix}.rescue.vcf.gz; touch ${prefix}.rescue.vcf.gz.tbi
    : > ${prefix}.hetcounts.tsv
    """
}
