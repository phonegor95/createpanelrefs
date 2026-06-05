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

workflow GERMLINECNVCALLER_COHORT {
    take:
    ch_input // channel: [mandatory] [ val(meta), path(bam/cram), path(bai/crai) ]
    val_pon_name //  string: [optional] name for panel of normals
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
    //  Prepare references
    GATK4_INDEXFEATUREFILE_MAPPABILITY(ch_mappable_regions)
    GATK4_INDEXFEATUREFILE_SEGDUP(ch_segmental_duplications)

    //Runs for wes analysis, when target_bed file is provided instead of target_interval_list
    GATK4_BEDTOINTERVALLIST_TARGETS(ch_target_bed, ch_dict)

    //Runs for wes analysis, when exclude_bed file is provided instead of target_interval_list
    GATK4_BEDTOINTERVALLIST_EXCLUDE(ch_exclude_bed, ch_dict)

    ch_user_target_interval_list
        .combine(GATK4_BEDTOINTERVALLIST_TARGETS.out.interval_list.ifEmpty(null))
        .branch { it ->
            intervallistfrompath: it[2].equals(null)
            return [it[0], it[1]]
            intervallistfrombed: !it[2].equals(null)
            return [it[2], it[3]]
        }
        .set { ch_targets_for_mix }

    ch_targets_for_mix.intervallistfrompath
        .mix(ch_targets_for_mix.intervallistfrombed)
        .collect()
        .set { ch_target_interval_list }

    ch_user_exclude_interval_list
        .combine(GATK4_BEDTOINTERVALLIST_EXCLUDE.out.interval_list.ifEmpty(null))
        .branch { it ->
            intervallistfrompath: it[2].equals(null)
            return [it[0], it[1]]
            intervallistfrombed: !it[2].equals(null)
            return [it[2], it[3]]
        }
        .set { ch_exclude_for_mix }

    ch_exclude_for_mix.intervallistfrompath
        .mix(ch_exclude_for_mix.intervallistfrombed)
        .collect()
        .set { ch_exclude_interval_list }

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

    // Filter out files that lack indices, and generate them
    ch_input
        .branch { meta, alignment, index ->
            alignment_with_index: index.size() > 0
            return [meta, alignment, index]
            alignment_without_index: index.size() == 0
            return [meta, alignment]
        }
        .set { ch_for_mix }

    SAMTOOLS_INDEX(ch_for_mix.alignment_without_index)

    ch_index = SAMTOOLS_INDEX.out.index

    // Collect alignment files and their indices
    ch_for_mix.alignment_without_index
        .join(ch_index)
        .mix(ch_for_mix.alignment_with_index)
        .combine(GATK4_PREPROCESSINTERVALS.out.interval_list.map { it -> it[1] })
        .set { ch_readcounts_in }

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
        .collect { row -> [sample_id: row[0].id, count_file: row[1]] }
        .map { rows ->
            def sorted_counts = rows
                .sort { a, b -> a.sample_id <=> b.sample_id }
                .collect { it.count_file }
            [[id: val_pon_name], sorted_counts]
        }
        .set { ch_readcounts_out }


    GATK4_FILTERINTERVALS(
        GATK4_PREPROCESSINTERVALS.out.interval_list,
        ch_readcounts_out,
        GATK4_ANNOTATEINTERVALS.out.annotated_intervals,
    )

    GATK4_INTERVALLISTTOOLS(GATK4_FILTERINTERVALS.out.interval_list).interval_list.map { _meta, it -> it }.flatten().set { ch_intervallist_out }

    ch_readcounts_out
        .combine(GATK4_FILTERINTERVALS.out.interval_list)
        .map { meta, counts, _meta2, il -> [meta, counts, il, []] }
        .set { ch_contigploidy_in }

    GATK4_DETERMINEGERMLINECONTIGPLOIDY(
        ch_contigploidy_in,
        [[:], []],
        ch_ploidy_priors,
    )

    ch_readcounts_out
        .combine(ch_intervallist_out)
        .combine(GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.calls)
        .combine(GATK4_ANNOTATEINTERVALS.out.annotated_intervals)
        .map { meta, counts, il, _meta2, calls, _meta3, annotated_intervals -> [meta + [id: il.baseName], counts, il, calls, [], annotated_intervals] }
        .set { ch_cnvcaller_in }

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
                [[id: sample_name, sample_index: sample_index], ploidy_calls, sample_index]
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

    emit:
    cnvmodel             = GATK4_GERMLINECNVCALLER.out.cohortmodel
    ploidymodel          = GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.model
    ploidycalls          = GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.calls
    readcounts           = ch_readcounts_out
    genotyped_intervals  = GATK4_POSTPROCESSGERMLINECNVCALLS.out.genotyped_intervals
    genotyped_segments   = GATK4_POSTPROCESSGERMLINECNVCALLS.out.genotyped_segments
    denoised_copy_ratios = GATK4_POSTPROCESSGERMLINECNVCALLS.out.denoised_copy_ratios
}
