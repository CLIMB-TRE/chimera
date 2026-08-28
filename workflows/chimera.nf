/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_chimera_pipeline'

include { RAMMAP_ALIGN           } from '../modules/local/rammap_align/rammap_align'
include { SYLPH_PROFILE          } from '../modules/nf-core/sylph/profile/main'
include { SAMTOOLS_SORT          } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX         } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_DEPTH         } from '../modules/nf-core/samtools/depth/main'

include { SYLPH_TAXONOMY         } from '../modules/local/sylph_taxonomy/sylph_taxonomy'
include { ALIGNMENT_REPORT       } from '../modules/local/alignment_report/alignment_report'
include { FILTER_BAM             } from '../modules/local/filter_bam/filter_bam'


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

    rammap_index = file(params.rammap_index, checkIfExists: true)
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

    RAMMAP_ALIGN(ch_samplesheet, [[:], rammap_index], "bai")
    ch_versions = ch_versions.mix(RAMMAP_ALIGN.out.versions.first())

    ch_aligned = RAMMAP_ALIGN.out.bam

    // RAMMAP_ALIGN emits coordinate-sorted+indexed BAM for every platform
    if (params.min_alignment_proportion_filter) {
        FILTER_BAM(ch_aligned)
        ch_for_downstream = FILTER_BAM.out.filtered_bam
    }
    else {
        ch_for_downstream = ch_aligned
    }

    // ONT: unfiltered/filtered single-end BAM stays coordinate-sorted — only indexing needed
    // Illumina: pair-grouping in FILTER_BAM (when enabled) disrupts coordinate order — always re-sort
    ch_branched = ch_for_downstream.branch { meta, _bam ->
        ont: meta.platform == "ont"
        illumina: meta.platform == "illumina" || meta.platform == "illumina.se"
    }

    SAMTOOLS_INDEX(ch_branched.ont)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    SAMTOOLS_SORT(ch_branched.illumina, [[:], []], "bai")
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    ch_ont_bam_bai = ch_branched.ont.join(SAMTOOLS_INDEX.out.bai, failOnDuplicate: true, failOnMismatch: true)
    ch_illumina_bam_bai = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_SORT.out.bai, failOnDuplicate: true, failOnMismatch: true)

    SAMTOOLS_DEPTH(
        ch_ont_bam_bai.mix(ch_illumina_bam_bai),
        [[:], []],
    )
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions.first())

    ch_alignment_report_input = SAMTOOLS_DEPTH.out.tsv.join(
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
