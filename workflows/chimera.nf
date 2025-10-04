/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_chimera_pipeline'

include { BWAMEM2_MEM            } from '../modules/nf-core/bwamem2/mem/main'
include { MINIMAP2_ALIGN         } from '../modules/nf-core/minimap2/align/main'
include { SYLPH_PROFILE          } from '../modules/nf-core/sylph/profile/main'
include { SAMTOOLS_SORT          } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX         } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_DEPTH         } from '../modules/nf-core/samtools/depth/main'
include { SAMTOOLS_COVERAGE      } from '../modules/nf-core/samtools/coverage/main'

include { SYLPH_TAXONOMY         } from '../modules/local/sylph_taxonomy'
include { ALIGNMENT_REPORT       } from '../modules/local/alignment_stats'


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

    SAMTOOLS_SORT(MINIMAP2_ALIGN.out.bam, [[:], []], "bai")
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    BWAMEM2_MEM(ch_samplesheet_branched.illumina, [[:], bwa_index], [[:], []], true)
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())

    ch_bams = SAMTOOLS_SORT.out.bam.mix(BWAMEM2_MEM.out.bam)

    SAMTOOLS_DEPTH(ch_bams, [[:], []])
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions.first())

    SAMTOOLS_INDEX(BWAMEM2_MEM.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    ch_bams_sorted_indexed = BWAMEM2_MEM.out.bam
        .join(SAMTOOLS_INDEX.out.bai, failOnDuplicate: true, failOnMismatch: true)
        .mix(
            SAMTOOLS_SORT.out.bam.join(
                SAMTOOLS_SORT.out.bai,
                failOnDuplicate: true,
                failOnMismatch: true,
            )
        )

    SAMTOOLS_COVERAGE(ch_bams_sorted_indexed, [[:], []], [[:], []])
    ch_versions = ch_versions.mix(SAMTOOLS_COVERAGE.out.versions.first())

    SYLPH_TAXONOMY(
        SYLPH_PROFILE.out.profile_out,
        sylph_taxdb,
    )

    ch_alignment_report_input = SAMTOOLS_DEPTH.out.tsv.join(SAMTOOLS_COVERAGE.out.coverage, failOnDuplicate: true, failOnMismatch: true)

    ALIGNMENT_REPORT(
        ch_alignment_report_input,
        database_metadata,
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
