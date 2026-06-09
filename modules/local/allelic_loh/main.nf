process ALLELIC_LOH {
    tag "${meta.id}"
    label 'process_medium'

    // GATK (CollectAllelicCounts/SelectVariants/VcfToIntervalList) + python3 in one image.
    container "${params.gatk_container ?: 'https://depot.galaxyproject.org/singularity/gatk4:4.6.1.0--py310hdfd78af_0'}"

    input:
    tuple val(meta), path(pass_vcf), path(cram), path(crai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)
    tuple path(snp_vcf), path(snp_tbi)

    output:
    tuple val(meta), path("*.loh.tsv"), emit: tsv
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = (task.memory ? (task.memory.giga * 0.85).intValue() : 6)
    """
    # 1. Deletion regions from the PASS VCF -> BED (0-based start). END from INFO, CN from FORMAT.
    zcat ${pass_vcf} | awk 'BEGIN{OFS="\t"} !/^#/ && \$5=="<DEL>"{
        end=\$2; n=split(\$8,info,";"); for(i=1;i<=n;i++){ if(info[i] ~ /^END=/){ e=info[i]; sub("END=","",e); end=e } }
        nf=split(\$9,fmt,":"); split(\$10,val,":"); cn=".";
        for(i=1;i<=nf;i++){ if(fmt[i]=="CN") cn=val[i] }
        print \$1, \$2-1, end, cn
    }' > dels.bed

    if [ ! -s dels.bed ]; then
        echo -e "sample\tchrom\tstart\tend\tCN\tcovered_SNPs\thet_sites\thet_rate\tLOH_status" > ${prefix}.loh.tsv
    else
        # 2. Common SNPs inside the deletion regions only (cheap, indexed).
        gatk --java-options "-Xmx${avail_mem}g" SelectVariants \\
            -V ${snp_vcf} -L dels.bed --interval-padding 0 -O del_snps.vcf.gz

        if [ \$(zcat del_snps.vcf.gz | grep -vc '^#') -eq 0 ]; then
            echo -e "sample\tchrom\tstart\tend\tCN\tcovered_SNPs\thet_sites\thet_rate\tLOH_status" > ${prefix}.loh.tsv
        else
            gatk --java-options "-Xmx${avail_mem}g" VcfToIntervalList -I del_snps.vcf.gz -O del_snps.interval_list
            # 3. Allele counts at those SNP sites from the sample CRAM.
            gatk --java-options "-Xmx${avail_mem}g" CollectAllelicCounts \\
                -I ${cram} -R ${fasta} -L del_snps.interval_list -O counts.tsv
            # 4. Annotate each deletion with LOH evidence (bin/ script on PATH).
            allelic_loh_annotate.py --del-bed dels.bed --counts counts.tsv \\
                --sample ${prefix} --out ${prefix}.loh.tsv
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n 's/^The Genome Analysis Toolkit (GATK) v//p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "sample\tchrom\tstart\tend\tCN\tcovered_SNPs\thet_sites\thet_rate\tLOH_status" > ${prefix}.loh.tsv
    echo '"${task.process}":' > versions.yml
    echo '    gatk4: 4.6.1.0' >> versions.yml
    """
}
