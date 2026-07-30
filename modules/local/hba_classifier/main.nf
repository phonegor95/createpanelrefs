process HBA_CLASSIFIER {
    tag "${meta.id}"
    label 'process_single'

    // GATK (CollectAllelicCounts) + samtools + python3 (stdlib) + perl in one image.
    container "${params.gatk_container ?: (workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/gatk4:4.6.1.0--py310hdfd78af_0'
        : 'biocontainers/gatk4:4.6.1.0--py310hdfd78af_0')}"

    input:
    tuple val(meta), path(cram), path(crai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)
    tuple path(snp_vcf), path(snp_tbi)

    output:
    tuple val(meta), path("*.HBA.yaml"),         emit: yaml
    tuple val(meta), path("*.HBA.crosscheck.tsv"), emit: crosscheck
    tuple val("${task.process}"), val('gatk4'),    eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | sed -n '1s/samtools //p'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Optional Sentieon segdup-caller YAML (a REFERENCE, not ground truth) for cross-check.
    def sent = params.hba_sentieon_dir ? "${params.hba_sentieon_dir}/${prefix}/${prefix}.yaml" : ''
    def sy = (sent && file(sent).exists()) ? "--sentieon-yaml ${sent}" : ''
    def avail_mem = (task.memory ? (task.memory.giga * 0.85).intValue() : 6)
    """
    # --- 1. region depths over the HBA cn_regions (hg38) -> depths.tsv ----------
    # %b (not %s) so printf interprets the \\t escapes into real tab characters
    printf '%b\\n' \\
      'chr16\\t172871\\t176674\\ta3.7' \\
      'chr16\\t169818\\t174075\\ta4.2' \\
      'chr16\\t174075\\t176674\\ta3.7_spec' \\
      'chr16\\t169818\\t172871\\ta4.2_spec' \\
      'chr16\\t169454\\t177522\\tlocus' \\
      'chr16\\t165199\\t169367\\tflank5' \\
      'chr16\\t177425\\t185000\\tflank3' > hba_regions.bed
    printf 'chr16\\t3000000\\t3050000\\tbase\\n' > hba_base.bed

    b=\$(samtools bedcov -Q 0 hba_base.bed ${cram} --reference ${fasta} | awk '{print \$5/(\$3-\$2)}')
    { printf 'base\\t%s\\n' "\$b"
      samtools bedcov -Q 0 hba_regions.bed ${cram} --reference ${fasta} \\
        | awk '{printf "%s\\t%.3f\\n",\$4,\$5/(\$3-\$2)}' ; } > ${prefix}.depths.tsv

    # --- 2. junction-anchored breakpoint evidence -> junc.tsv ------------------
    hba_breakpoints.sh ${cram} ${fasta} ${prefix}.junc.tsv

    # --- 2b. PSV HBA2/HBA1 paralog-specific dosage pileup -> psv.pileup --------
    # advisory HBA2-dosage axis (paralog-immune); -q 0 keeps cross-mapped reads.
    samtools mpileup -f ${fasta} -q 0 -Q 20 -r chr16:173300-177560 ${cram} \\
        > ${prefix}.psv.pileup 2>/dev/null || true

    # --- 3. allelic counts over the HBA locus (common SNPs, MQ30) -> ac.tsv ----
    gatk --java-options "-Xmx${avail_mem}g" SelectVariants \\
        -V ${snp_vcf} -L chr16:169000-178000 -O hba_snps.vcf.gz
    gatk --java-options "-Xmx${avail_mem}g" CollectAllelicCounts \\
        -I ${cram} -R ${fasta} -L hba_snps.vcf.gz \\
        --minimum-mapping-quality 30 -O ${prefix}.ac.tsv

    # --- 4. classify + cross-check -> YAML + crosscheck ------------------------
    hba_classifier.py --sample ${prefix} \\
        --depths ${prefix}.depths.tsv \\
        --allelic-counts ${prefix}.ac.tsv \\
        --breakpoints ${prefix}.junc.tsv ${sy} \\
        --psv-pileup ${prefix}.psv.pileup \\
        --variants-vcf ${prefix}.ac.tsv \\
        --out-yaml ${prefix}.HBA.yaml \\
        --out-crosscheck ${prefix}.HBA.crosscheck.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'HBA:\\n  Copy numbers:\\n    a3.7: 2\\n    a4.2: 2\\n    non-duplication: 2\\n  evidence:\\n    psv_dosage: null\\n' > ${prefix}.HBA.yaml
    printf '# HBA cross-check sample=%s\\n' ${prefix} > ${prefix}.HBA.crosscheck.tsv
    """
}
