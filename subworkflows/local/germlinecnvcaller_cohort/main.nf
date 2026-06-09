include { GATK4_ANNOTATEINTERVALS                                      } from '../../../modules/nf-core/gatk4/annotateintervals'
include { GATK4_BEDTOINTERVALLIST as GATK4_BEDTOINTERVALLIST_TARGETS   } from '../../../modules/nf-core/gatk4/bedtointervallist'
include { GATK4_BEDTOINTERVALLIST as GATK4_BEDTOINTERVALLIST_EXCLUDE   } from '../../../modules/nf-core/gatk4/bedtointervallist'
include { GATK4_COLLECTREADCOUNTS                                      } from '../../../modules/nf-core/gatk4/collectreadcounts'
include { GATK4_DETERMINEGERMLINECONTIGPLOIDY                          } from '../../../modules/nf-core/gatk4/determinegermlinecontigploidy'
include { GATK4_FILTERINTERVALS                                        } from '../../../modules/nf-core/gatk4/filterintervals'
include { GATK4_GERMLINECNVCALLER                                      } from '../../../modules/nf-core/gatk4/germlinecnvcaller'
include { GATK4_INDEXFEATUREFILE as GATK4_INDEXFEATUREFILE_MAPPABILITY } from '../../../modules/nf-core/gatk4/indexfeaturefile'
include { GATK4_INDEXFEATUREFILE as GATK4_INDEXFEATUREFILE_SEGDUP      } from '../../../modules/nf-core/gatk4/indexfeaturefile'
include { GATK4_INTERVALLISTTOOLS                                      } from '../../../modules/nf-core/gatk4/intervallisttools'
include { GATK4_POSTPROCESSGERMLINECNVCALLS                            } from '../../../modules/nf-core/gatk4/postprocessgermlinecnvcalls'
include { GATK4_PREPROCESSINTERVALS                                    } from '../../../modules/nf-core/gatk4/preprocessintervals'
include { SAMTOOLS_INDEX                                               } from '../../../modules/nf-core/samtools/index'
include { BCFTOOLS_FILTER_GCNV                                         } from '../../../modules/local/bcftools_filter_gcnv'
include { COHORT_RECURRENCE_FILTER                                     } from '../../../modules/local/cohort_recurrence_filter'
include { CNV_QC_OUTLIER                                               } from '../../../modules/local/cnv_qc_outlier'
include { ALLELIC_LOH                                                  } from '../../../modules/local/allelic_loh'
include { ANNOTSV                                                      } from '../../../modules/local/annotsv'
include { KNOTANNOTSV as KNOTANNOTSV_HTML                              } from '../../../modules/nf-core/knotannotsv'
include { KNOTANNOTSV as KNOTANNOTSV_XLSM                              } from '../../../modules/nf-core/knotannotsv'

workflow GERMLINECNVCALLER_COHORT {
    take:
    ch_reads_index // channel: [mandatory] [ val(meta), path(bam/cram), path(bai/crai) ]
    val_pon_name //  string: [optional] name for panel of normals
    val_analysis_type // string: [mandatory] type of analysis ('wes' or 'wgs')
    ch_dict // channel: [optional] [ val(meta), path(dict) ]
    ch_fai // channel: [optional] [ val(meta), path(fai) ]
    ch_fasta // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_exclude_bed // channel: [optional] [ val(meta), path(bed) ]
    ch_user_exclude_interval_list // channel: [optional] [ val(meta), path(intervals) ]
    ch_mappable_regions // channel: [optional] [ val(meta), path(bed) ]
    ch_ploidy_priors // channel: [mandatory] [ path(tsv) ]
    ch_segmental_duplications // channel: [optional] [ val(meta), path(bed) ]
    ch_target_bed // channel: [optional] [ val(meta), path(bed) ]
    ch_user_target_interval_list // channel: [optional] [ val(meta), path(intervals) ]

    main:
    //  Index feature files — only when a real file is provided
    GATK4_INDEXFEATUREFILE_MAPPABILITY(ch_mappable_regions.filter { _meta, regions -> !(regions instanceof List) })
    GATK4_INDEXFEATUREFILE_SEGDUP(ch_segmental_duplications.filter { _meta, segdup -> !(segdup instanceof List) })

    // Bed to interval list conversion — only for WES when bed is provided and no interval list given
    ch_target_bed_interval_list = channel.empty()
    ch_exclude_bed_interval_list = channel.empty()

    if (val_analysis_type == "wes") {
        GATK4_BEDTOINTERVALLIST_TARGETS(
            ch_target_bed.filter { _meta, bed -> !(bed instanceof List) },
            ch_dict,
        )

        GATK4_BEDTOINTERVALLIST_EXCLUDE(
            ch_exclude_bed.filter { _meta, bed -> !(bed instanceof List) },
            ch_dict,
        )

        ch_target_bed_interval_list = GATK4_BEDTOINTERVALLIST_TARGETS.out.interval_list
        ch_exclude_bed_interval_list = GATK4_BEDTOINTERVALLIST_EXCLUDE.out.interval_list
    }

    ch_targets_for_mix = ch_user_target_interval_list
        .combine(ch_target_bed_interval_list.ifEmpty(null))
        .branch { it ->
            intervallistfrompath: it[2].equals(null)
            return [it[0], it[1]]
            intervallistfrombed: !it[2].equals(null)
            return [it[2], it[3]]
        }

    ch_target_interval_list = ch_targets_for_mix.intervallistfrompath
        .mix(ch_targets_for_mix.intervallistfrombed)
        .collect()

    ch_exclude_for_mix = ch_user_exclude_interval_list
        .combine(ch_exclude_bed_interval_list.ifEmpty(null))
        .branch { it ->
            intervallistfrompath: it[2].equals(null)
            return [it[0], it[1]]
            intervallistfrombed: !it[2].equals(null)
            return [it[2], it[3]]
        }

    ch_exclude_interval_list = ch_exclude_for_mix.intervallistfrompath
        .mix(ch_exclude_for_mix.intervallistfrombed)
        .collect()

    GATK4_PREPROCESSINTERVALS(
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_target_interval_list,
        ch_exclude_interval_list,
    )

    GATK4_ANNOTATEINTERVALS(
        GATK4_PREPROCESSINTERVALS.out.interval_list,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_mappable_regions,
        GATK4_INDEXFEATUREFILE_MAPPABILITY.out.index.ifEmpty([[:], []]),
        ch_segmental_duplications,
        GATK4_INDEXFEATUREFILE_SEGDUP.out.index.ifEmpty([[:], []]),
    )

    ch_readcounts_in = ch_reads_index.combine(GATK4_PREPROCESSINTERVALS.out.interval_list.map { it -> it[1] })

    // Collect read counts, and generate models
    GATK4_COLLECTREADCOUNTS(
        ch_readcounts_in,
        ch_fasta,
        ch_fai,
        ch_dict,
    )

    // Sort the read-count files by sample id so every scattered GermlineCNVCaller
    // receives them in an identical order. The per-shard sample order must match
    // across all shards, otherwise PostprocessGermlineCNVCalls rejects them with
    // "The sample name is not the same for all of the shards". Upstream
    // createpanelrefs stops at the cohort model and never needs this; our
    // postprocess fan-out does.
    GATK4_COLLECTREADCOUNTS.out.tsv
        .mix(GATK4_COLLECTREADCOUNTS.out.hdf5)
        .map { _meta, count_file -> count_file }
        .toSortedList { a, b -> a.name <=> b.name }
        .map { sorted_counts -> [[id: val_pon_name], sorted_counts] }
        .set { ch_readcounts_out }


    GATK4_FILTERINTERVALS(
        GATK4_PREPROCESSINTERVALS.out.interval_list,
        ch_readcounts_out,
        GATK4_ANNOTATEINTERVALS.out.annotated_intervals,
    )

    ch_intervallist_out = GATK4_INTERVALLISTTOOLS(GATK4_FILTERINTERVALS.out.interval_list).interval_list.map { _meta, it -> it }.flatten()

    ch_contigploidy_in = ch_readcounts_out
        .combine(GATK4_FILTERINTERVALS.out.interval_list)
        .map { meta, counts, _meta2, il -> [meta, counts, il, []] }

    GATK4_DETERMINEGERMLINECONTIGPLOIDY(
        ch_contigploidy_in,
        [[:], []],
        ch_ploidy_priors,
    )

    ch_cnvcaller_in = ch_readcounts_out
        .combine(ch_intervallist_out)
        .combine(GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.calls)
        .combine(GATK4_ANNOTATEINTERVALS.out.annotated_intervals)
        .map { meta, counts, il, _meta2, calls, _meta3, annotated_intervals -> [meta + [id: il.baseName], counts, il, calls, [], annotated_intervals] }

    GATK4_GERMLINECNVCALLER(ch_cnvcaller_in)

    GATK4_GERMLINECNVCALLER.out.cohortmodel
        .map { _meta, model_dir -> model_dir }
        .collect()
        .map { model_dirs -> [model_shards: model_dirs] }
        .set { ch_gcnv_model_shards }

    GATK4_GERMLINECNVCALLER.out.cohortcalls
        .map { _meta, calls_dir -> calls_dir }
        .collect()
        .map { call_dirs -> [call_shards: call_dirs] }
        .set { ch_gcnv_call_shards }

    GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.calls
        .flatMap { _meta, ploidy_calls ->
            def sample_dirs = ploidy_calls.toFile()
                .listFiles()
                .findAll { it.isDirectory() && it.name.startsWith('SAMPLE_') }
                .sort { a, b -> (a.name - 'SAMPLE_') as Integer <=> (b.name - 'SAMPLE_') as Integer }
            sample_dirs.collect { sample_dir ->
                def sample_index = (sample_dir.name - 'SAMPLE_') as Integer
                def sample_name = new File(sample_dir, 'sample_name.txt').text.trim()
                // Infer sex from the inferred chrX ploidy so downstream segment
                // filtering can apply sex-aware thresholds without a samplesheet
                // column: chrX ploidy >= 2 => female, otherwise male.
                def ploidy_row = new File(sample_dir, 'contig_ploidy.tsv').readLines()
                    .find { line -> def c = line.split('\t')[0]; c == 'chrX' || c == 'X' }
                def chrx_ploidy = ploidy_row ? (ploidy_row.split('\t')[1] as Integer) : null
                def sex = (chrx_ploidy != null && chrx_ploidy >= 2) ? 'F' : 'M'
                [[id: sample_name, sample_index: sample_index, sex: sex], ploidy_calls, sample_index]
            }
        }
        .combine(ch_gcnv_model_shards)
        .combine(ch_gcnv_call_shards)
        .combine(ch_dict)
        .map { row ->
            def meta = row[0]
            def ploidy_calls = row[1]
            def sample_index = row[2]
            def model_shards = row[3].model_shards
            def call_shards = row[4].call_shards
            def dict_file = row[6]
            [meta, model_shards, call_shards, ploidy_calls, dict_file, sample_index]
        }
        .set { ch_postprocess_in }

    GATK4_POSTPROCESSGERMLINECNVCALLS(ch_postprocess_in)

    // Sex-aware filtering of the per-sample genotyped segments (meta.sex set above).
    BCFTOOLS_FILTER_GCNV(GATK4_POSTPROCESSGERMLINECNVCALLS.out.genotyped_segments)

    // ---- Cohort-recurrence (two-arm) soft filter ----------------------------
    // Pool every sample's sex-filtered segments into one cohort process that
    // flags recurrent calls (common CNVs / paralog artifacts). Needs all
    // samples at once, so collect the VCFs and build a sample->sex map. The
    // sex map drives the per-sex denominators for chrX/chrY (autosomes pooled).
    ch_filtered_for_recur = BCFTOOLS_FILTER_GCNV.out.vcf

    ch_recur_sexmap = ch_filtered_for_recur
        .map { meta, _vcf -> "${meta.id}\t${meta.sex}\n" }
        .collectFile(name: 'cohort_sex_map.tsv', sort: true)

    ch_recur_vcfs = ch_filtered_for_recur
        .map { _meta, vcf -> vcf }
        .collect()

    // Nest the collected VCF list inside the tuple BEFORE joining, otherwise
    // combine/merge spread the list into separate tuple elements. Appending the
    // scalar sex_map via combine then keeps the list intact.
    ch_recur_in = ch_recur_vcfs
        .map { vcfs -> [[id: val_pon_name], vcfs] }
        .combine(ch_recur_sexmap)
        .map { meta, vcfs, sex_map -> [meta, sex_map, vcfs] }

    COHORT_RECURRENCE_FILTER(ch_recur_in)

    // ---- Per-sample CNV QC outlier flag (cohort report) ---------------------
    ch_qc_in = ch_recur_vcfs.map { vcfs -> [[id: val_pon_name], vcfs] }
    ch_qc_segdup = ch_segmental_duplications.map { _meta, bed -> bed }.ifEmpty([])
    ch_qc_blacklist = params.qc_blacklist_bed ? file(params.qc_blacklist_bed, checkIfExists: true) : []
    CNV_QC_OUTLIER(ch_qc_in, ch_qc_segdup, ch_qc_blacklist)

    // ---- Allelic LOH confirmation of PASS deletions (orthogonal to depth) ----
    // Gated on params.allelic_snp_vcf. Join each PASS VCF (doubled sample id)
    // back to its alignment (clean meta.id) for CollectAllelicCounts.
    ch_loh_pass = COHORT_RECURRENCE_FILTER.out.pass
        .map { _meta, pass_files -> (pass_files instanceof List) ? pass_files : [pass_files] }
        .flatten()
        .map { f -> [f.name.replaceAll(/\.recurfilt\.pass\.vcf\.gz$/, ''), f] }

    ch_loh_in = ch_loh_pass
        .combine(ch_reads_index.map { m, aln, idx -> [m.id, aln, idx] })
        .filter { doubled, _vcf, sid, _aln, _idx -> doubled == sid || doubled == "${sid}_${sid}" }
        .map { doubled, vcf, _sid, aln, idx -> [[id: doubled], vcf, aln, idx] }

    ch_snp_loh = params.allelic_snp_vcf
        ? channel.value([file(params.allelic_snp_vcf, checkIfExists: true),
                         file(params.allelic_snp_vcf + '.tbi', checkIfExists: true)])
        : channel.value([[], []])

    ch_loh_segdup = ch_segmental_duplications.map { _meta, bed -> bed }.ifEmpty([]).first()

    ALLELIC_LOH(ch_loh_in, ch_fasta, ch_fai, ch_dict, ch_snp_loh, ch_loh_segdup)

    // ---- Per-sample AnnotSV annotation of the survivors ---------------------
    // When allelic LOH ran, annotate its hard-filtered output (no_LOH deletions
    // removed); otherwise annotate the recurrence-filter PASS calls directly.
    ch_annotsv_in = params.allelic_snp_vcf
        ? ALLELIC_LOH.out.lohpass
        : COHORT_RECURRENCE_FILTER.out.pass
            .map { _meta, pass_files -> (pass_files instanceof List) ? pass_files : [pass_files] }
            .flatten()
            .map { f -> [[id: f.name.replaceAll(/\.recurfilt\.pass\.vcf\.gz$/, '')], f] }

    ch_annot_dir = params.annotsv_annotations
        ? channel.value(file(params.annotsv_annotations, checkIfExists: true))
        : channel.value([])

    ANNOTSV(ch_annotsv_in, ch_annot_dir)

    // ---- knotAnnotSV: HTML report + XLSM workbook per sample ----------------
    // 3rd tuple value selects the script (false = HTML, true = XLSM).
    KNOTANNOTSV_HTML(ANNOTSV.out.tsv.map { meta, tsv -> [meta, tsv, false] })
    KNOTANNOTSV_XLSM(ANNOTSV.out.tsv.map { meta, tsv -> [meta, tsv, true] })

    emit:
    // Upstream 1.0.0 emit names consumed by workflows/createpanelrefs.nf
    cnv_model            = GATK4_GERMLINECNVCALLER.out.cohortmodel
    cnv_calls            = GATK4_GERMLINECNVCALLER.out.cohortcalls
    ploidy_model         = GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.model
    read_counts          = ch_readcounts_out
    // Local postprocess fan-out outputs (published via publishDir in conf/modules)
    ploidycalls          = GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.calls
    genotyped_intervals  = GATK4_POSTPROCESSGERMLINECNVCALLS.out.genotyped_intervals
    genotyped_segments   = GATK4_POSTPROCESSGERMLINECNVCALLS.out.genotyped_segments
    denoised_copy_ratios = GATK4_POSTPROCESSGERMLINECNVCALLS.out.denoised_copy_ratios
    filtered_segments    = BCFTOOLS_FILTER_GCNV.out.vcf
    recurfilt_flagged    = COHORT_RECURRENCE_FILTER.out.flagged
    recurfilt_pass       = COHORT_RECURRENCE_FILTER.out.pass
    annotsv              = ANNOTSV.out.tsv
    knotannotsv_html     = KNOTANNOTSV_HTML.out.html
    knotannotsv_xlsm     = KNOTANNOTSV_XLSM.out.xl
    cnv_qc               = CNV_QC_OUTLIER.out.report
    allelic_loh          = ALLELIC_LOH.out.tsv
    allelic_loh_pass     = ALLELIC_LOH.out.lohpass
}
