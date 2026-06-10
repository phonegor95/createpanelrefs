process ALLELIC_LOH {
    tag "${meta.id}"
    label 'process_medium'

    // GATK (CollectAllelicCounts/SelectVariants/VcfToIntervalList) + bgzip/tabix + python3 in one image.
    container "${params.gatk_container ?: 'https://depot.galaxyproject.org/singularity/gatk4:4.6.1.0--py310hdfd78af_0'}"

    input:
    tuple val(meta), path(pass_vcf), path(cram), path(crai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)
    tuple path(snp_vcf), path(snp_tbi)
    path segdup_bed

    output:
    tuple val(meta), path("*.loh.tsv"),                emit: tsv
    tuple val(meta), path("*.lohpass.vcf.gz"),         emit: lohpass
    tuple val(meta), path("*.lohpass.vcf.gz.tbi"),     emit: lohpass_tbi
    tuple val(meta), path("*.lohflagged.vcf.gz"),      emit: flagged
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def hard = (params.allelic_loh_hard_filter == null || params.allelic_loh_hard_filter) ? '--hard-filter' : ''
    def sd = segdup_bed ? "--segdup-bed ${segdup_bed}" : ''
    def avail_mem = (task.memory ? (task.memory.giga * 0.85).intValue() : 6)
    """
    # 1. Deletion regions from the PASS VCF -> BED (0-based start).
    zcat ${pass_vcf} | awk 'BEGIN{OFS="\t"} !/^#/ && \$5=="<DEL>"{
        end=""; n=split(\$8,info,";"); for(i=1;i<=n;i++){ if(info[i] ~ /^END=/){ e=info[i]; sub("END=","",e); end=e } }
        if(end=="" || end+0 <= \$2){ print "[ALLELIC_LOH] skipping DEL without valid END: "\$1":"\$2 > "/dev/stderr"; next }
        nf=split(\$9,fmt,":"); split(\$10,val,":"); cn=".";
        for(i=1;i<=nf;i++){ if(fmt[i]=="CN") cn=val[i] }
        print \$1, \$2-1, end, cn
    }' > dels.bed

    : > counts.tsv
    if [ -s dels.bed ]; then
        gatk --java-options "-Xmx${avail_mem}g" SelectVariants \\
            -V ${snp_vcf} -L dels.bed --interval-padding 0 -O del_snps.vcf.gz
        if [ \$(zcat del_snps.vcf.gz | grep -vc '^#') -gt 0 ]; then
            gatk --java-options "-Xmx${avail_mem}g" VcfToIntervalList -I del_snps.vcf.gz -O del_snps.interval_list
            gatk --java-options "-Xmx${avail_mem}g" CollectAllelicCounts \\
                -I ${cram} -R ${fasta} -L del_snps.interval_list -O counts.tsv
        fi
    fi

    # 2. Annotate + hard-filter (segdup-overlapping DELs are exempt from filtering).
    allelic_loh_annotate.py \\
        --del-bed dels.bed --counts counts.tsv --sample ${prefix} \\
        --in-vcf ${pass_vcf} ${sd} ${hard} \\
        --out-tsv ${prefix}.loh.tsv \\
        --out-flagged ${prefix}.lohflagged.vcf \\
        --out-pass ${prefix}.lohpass.vcf

    bgzip -f ${prefix}.lohflagged.vcf && tabix -p vcf -f ${prefix}.lohflagged.vcf.gz
    bgzip -f ${prefix}.lohpass.vcf    && tabix -p vcf -f ${prefix}.lohpass.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "sample\tchrom\tstart\tend\tCN\tcovered_SNPs\thet_sites\thet_rate\tLOH_status" > ${prefix}.loh.tsv
    echo | gzip -c > ${prefix}.lohpass.vcf.gz;    touch ${prefix}.lohpass.vcf.gz.tbi
    echo | gzip -c > ${prefix}.lohflagged.vcf.gz
    """
}
