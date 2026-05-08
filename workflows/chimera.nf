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
include { SAMTOOLS_SORT  as SAMTOOLS_SORT_2 } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX                    } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_DEPTH                    } from '../modules/nf-core/samtools/depth/main'

include { SYLPH_TAXONOMY                   } from '../modules/local/sylph_taxonomy/sylph_taxonomy'
include { ALIGNMENT_REPORT                 } from '../modules/local/alignment_report/alignment_report'
include { FILTER_BAM                       } from '../modules/local/filter_bam/filter_bam'


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
    sylph_db = file(params.sylph_db, checkIfExists: true)
    sylph_taxdb = file(params.sylph_taxdb, checkIfExists: true)

    //
    // Run slyph and alignments to reference db
    //

    SYLPH_PROFILE(
        ch_samplesheet,
        sylph_db,
    )
    ch_versions = ch_versions.mix(SYLPH_PROFILE.out.versions.first())

    // Run the appropriate aligner based on platform
    ch_samplesheet_branched = ch_samplesheet.branch { meta, _fastq ->
        ont: meta.platform == "ont"
        illumina: meta.platform == "illumina" || meta.platform == "illumina.se"
    }

    MINIMAP2_ALIGN(ch_samplesheet_branched.ont, [[:], mm2_index], true, "bai", false, false)
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())

    BWAMEM2_MEM(ch_samplesheet_branched.illumina, [[:], bwa_index], [[:], []], true)
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())

    // MINIMAP2_ALIGN already emits coordinate-sorted+indexed BAM; BWAMEM2_MEM emits sorted BAM
    FILTER_BAM(MINIMAP2_ALIGN.out.bam.mix(BWAMEM2_MEM.out.bam))

    // ONT: filter preserves coordinate order for single-end reads — only indexing needed
    // Illumina: pair-grouping in FILTER_BAM disrupts coordinate order — must re-sort
    ch_filtered_branched = FILTER_BAM.out.filtered_bam.branch { meta, _bam ->
        ont:      meta.platform == "ont"
        illumina: meta.platform == "illumina" || meta.platform == "illumina.se"
    }

    SAMTOOLS_INDEX(ch_filtered_branched.ont)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    SAMTOOLS_SORT_2(ch_filtered_branched.illumina, [[:], []], "bai")
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_2.out.versions.first())

    ch_ont_bam_bai = ch_filtered_branched.ont
        .join(SAMTOOLS_INDEX.out.bai, failOnDuplicate: true, failOnMismatch: true)
    ch_illumina_bam_bai = SAMTOOLS_SORT_2.out.bam
        .join(SAMTOOLS_SORT_2.out.bai, failOnDuplicate: true, failOnMismatch: true)

    SAMTOOLS_DEPTH(
        ch_ont_bam_bai.mix(ch_illumina_bam_bai),
        [[:], []]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions.first())

    SYLPH_TAXONOMY(
        SYLPH_PROFILE.out.profile_out,
        sylph_taxdb,
    )

    ch_alignment_report_input = SAMTOOLS_DEPTH.out.tsv
        .join(
            ch_filtered_branched.ont.mix(SAMTOOLS_SORT_2.out.bam),
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
