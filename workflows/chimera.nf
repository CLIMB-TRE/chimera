/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap                 } from 'plugin/nf-schema'
include { softwareVersionsToYAML           } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText           } from '../subworkflows/local/utils_nfcore_chimera_pipeline'

include { BWAMEM2_MEM                      } from '../modules/nf-core/bwamem2/mem/main'
include { MINIMAP2_ALIGN                   } from '../modules/nf-core/minimap2/align/main'
include { SYLPH_PROFILE                    } from '../modules/nf-core/sylph/profile/main'
include { SAMTOOLS_SORT                     } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX                    } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_DEPTH                    } from '../modules/nf-core/samtools/depth/main'

include { SYLPH_TAXONOMY                   } from '../modules/local/sylph_taxonomy/sylph_taxonomy'
include { ALIGNMENT_REPORT                 } from '../modules/local/alignment_report/alignment_report'
include { FILTER_BAM                       } from '../modules/local/filter_bam/filter_bam'
include { FILL_SECONDARY_SEQ               } from '../modules/local/fill_secondary_seq/fill_secondary_seq'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CHIMERA {
    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()

    mm2_index = file(params.mm2_index, checkIfExists: true)
    bwa_index = file("${params.bwa_index_prefix}*")
    database_metadata = file(params.database_metadata, checkIfExists: true)

    //
    // Run sylph profiling and taxonomy reporting, unless skipped
    //

    if (!params.skip_sylph) {
        sylph_db = file(params.sylph_db, checkIfExists: true)
        sylph_taxdb = file(params.sylph_taxdb, checkIfExists: true)

        SYLPH_PROFILE(
            ch_samplesheet,
            sylph_db,
        )
        ch_versions = ch_versions.mix(SYLPH_PROFILE.out.versions.first())

        SYLPH_TAXONOMY(
            SYLPH_PROFILE.out.profile_out,
            sylph_taxdb,
        )
    }

    // Run the appropriate aligner based on platform
    ch_samplesheet_branched = ch_samplesheet.branch { meta, _fastq ->
        ont: meta.platform == "ont"
        illumina: meta.platform == "illumina" || meta.platform == "illumina.se"
    }

    MINIMAP2_ALIGN(ch_samplesheet_branched.ont, [[:], mm2_index], true, "bai", false, false)
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())

    BWAMEM2_MEM(ch_samplesheet_branched.illumina, [[:], bwa_index], [[:], []], true)
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())

    // Expand bwa-mem2 XA tags into real secondary alignment records and fill SEQ/QUAL
    // on any secondary records missing sequence, so Illumina secondaries carry sequence
    // the same way minimap2's --secondary-seq already does for ONT.
    FILL_SECONDARY_SEQ(BWAMEM2_MEM.out.bam)

    ch_aligned = MINIMAP2_ALIGN.out.bam.mix(FILL_SECONDARY_SEQ.out.bam)

    // MINIMAP2_ALIGN already emits coordinate-sorted+indexed BAM; FILL_SECONDARY_SEQ emits sorted BAM
    if (params.min_alignment_proportion_filter) {
        FILTER_BAM(ch_aligned)
        ch_for_downstream = FILTER_BAM.out.filtered_bam
    } else {
        ch_for_downstream = ch_aligned
    }

    // ONT: unfiltered/filtered single-end BAM stays coordinate-sorted — only indexing needed
    // Illumina: pair-grouping in FILTER_BAM (when enabled) disrupts coordinate order — always re-sort
    ch_branched = ch_for_downstream.branch { meta, _bam ->
        ont:      meta.platform == "ont"
        illumina: meta.platform == "illumina" || meta.platform == "illumina.se"
    }

    SAMTOOLS_INDEX(ch_branched.ont)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    SAMTOOLS_SORT(ch_branched.illumina, [[:], []], "bai")
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    ch_ont_bam_bai = ch_branched.ont
        .join(SAMTOOLS_INDEX.out.bai, failOnDuplicate: true, failOnMismatch: true)
    ch_illumina_bam_bai = SAMTOOLS_SORT.out.bam
        .join(SAMTOOLS_SORT.out.bai, failOnDuplicate: true, failOnMismatch: true)

    SAMTOOLS_DEPTH(
        ch_ont_bam_bai.mix(ch_illumina_bam_bai),
        [[:], []]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions.first())

    ch_alignment_report_input = SAMTOOLS_DEPTH.out.tsv
        .join(
            ch_branched.ont.mix(SAMTOOLS_SORT.out.bam),
            failOnDuplicate: true,
            failOnMismatch: true,
        )

    scoring_matrix = file(params.alignment_scoring_matrix, checkIfExists: true)
    json_schema = file(params.alignment_scoring_json_schema, checkIfExists: true)

    ALIGNMENT_REPORT(
        ch_alignment_report_input,
        database_metadata,
        scoring_matrix,
        json_schema,
    )


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'chimera_software_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    emit:
    versions = ch_versions // channel: [ path(versions.yml) ]
}
