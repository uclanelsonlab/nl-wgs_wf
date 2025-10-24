nextflow.enable.dsl = 2

include { FASTQ_PROCESSING } from './subworkflows/local/fastq/main'
include { BAM_PROCESSING } from './subworkflows/local/bam/main'
include { CRAM_PROCESSING } from './subworkflows/local/cram/main'
include { COMMON_ANALYSIS } from './subworkflows/local/common/main'

log.info """\
    SHORTREAD - WGS _ W F   P I P E L I N E
    ===================================
    sample_name  : ${params.sample_name}
    input_type   : ${params.input_type}
    input_1      : ${params.input_1}
    input_2      : ${params.input_2}
    fasta        : ${params.fasta}
    fai          : ${params.fai}
    mt_bed       : ${params.mt_bed}
    """
    .stripIndent(true)

workflow {
    // Prepare reference channels
    ch_fasta = Channel
        .fromPath([params.fasta, params.fai])
        .collect()
        .map { files ->
            def meta = [:]
            meta.id = files[0].baseName
            [ meta, files ]
        }

    ch_fasta_dict = Channel
        .fromPath([params.fasta, params.fai, params.dict])
        .collect()
        .map { files ->
            def meta = [:]
            meta.id = files[0].baseName
            [ meta, files[0], files[1], files[2] ]
        }

    ch_index = Channel
        .fromPath("${params.fasta}.*")
        .filter { it.name.endsWith('.amb') || it.name.endsWith('.ann') || it.name.endsWith('.bwt.2bit.64') || it.name.endsWith('.pac') || it.name.endsWith('.0123') }
        .collect()
        .map { files ->
            def meta = [:]
            meta.id = files[0].baseName
            [ meta, files ]
        }

    // Initialize channels
    ch_sorted_bam = Channel.empty()
    ch_fastp_html = Channel.empty()
    ch_fastp_json = Channel.empty()

    // Route to appropriate subworkflow based on input type
    if (params.input_type == "fastq") {
        // Create FASTQ input channel
        ch_reads = Channel
            .fromPath(params.input_1)
            .map { r1_file ->
                def meta = [:]
                meta.id = params.sample_name ?: r1_file.baseName
                def r2_file = file(params.input_2)
                [ meta, [r1_file, r2_file] ]
            }

        // Debug input channels
        ch_reads.view { "DEBUG: ch_reads = $it" }
        ch_index.view { "DEBUG: ch_index = $it" }
        ch_fasta.view { "DEBUG: ch_fasta = $it" }

        FASTQ_PROCESSING(
            ch_reads,
            ch_index,
            ch_fasta
        )
        
        // Debug outputs from FASTQ_PROCESSING
        FASTQ_PROCESSING.out.sorted_bam.view { "DEBUG: FASTQ_PROCESSING.out.sorted_bam = $it" }
        FASTQ_PROCESSING.out.fastp_html.view { "DEBUG: FASTQ_PROCESSING.out.fastp_html = $it" }
        FASTQ_PROCESSING.out.fastp_json.view { "DEBUG: FASTQ_PROCESSING.out.fastp_json = $it" }
        
        ch_sorted_bam = FASTQ_PROCESSING.out.sorted_bam
        ch_fastp_html = FASTQ_PROCESSING.out.fastp_html
        ch_fastp_json = FASTQ_PROCESSING.out.fastp_json

    } else if (params.input_type == "bam") {
        // Create BAM input channel
        if (params.input_2) {
            ch_input_bam = Channel
                .fromPath([params.input_1, params.input_2])
                .collect()
                .map { files ->
                    def meta = [:]
                    meta.id = params.sample_name ?: files[0].baseName
                    [ meta, files[0], files[1] ]
                }
        } else {
            ch_input_bam = Channel
                .fromPath(params.input_1)
                .map { bam_file ->
                    def meta = [:]
                    meta.id = params.sample_name ?: bam_file.baseName
                    [ meta, bam_file ]
                }
        }

        BAM_PROCESSING(
            ch_input_bam,
            ch_fasta
        )
        
        ch_sorted_bam = BAM_PROCESSING.out.sorted_bam

    } else if (params.input_type == "cram") {
       // Create CRAM input channel
        if (params.input_2) {
            ch_input_cram = Channel
                .fromPath([params.input_1, params.input_2])
                .collect()
                .map { files ->
                    def meta = [:]
                    meta.id = params.sample_name ?: files[0].baseName
                    [ meta, files[0], files[1] ]
                }
        } else {
            ch_input_cram = Channel
                .fromPath(params.input_1)
                .map { cram_file ->
                    def meta = [:]
                    meta.id = params.sample_name ?: cram_file.baseName
                    [ meta, cram_file ]
                }
        }

        CRAM_PROCESSING(
            ch_input_cram,
            ch_fasta
        )
        
        ch_sorted_bam = CRAM_PROCESSING.out.sorted_bam

    } else {
        error "Invalid input_type: ${params.input_type}. Must be 'fastq', 'bam', or 'cram'"
    }

    // Debug final channels before COMMON_ANALYSIS
    ch_sorted_bam.view { "DEBUG: Final ch_sorted_bam = $it" }
    ch_fasta.view { "DEBUG: Final ch_fasta = $it" }
    ch_fasta_dict.view { "DEBUG: Final ch_fasta_dict = $it" }

    // Count items in channels
    ch_sorted_bam.count().view { "DEBUG: ch_sorted_bam count = $it" }
    ch_fasta.count().view { "DEBUG: ch_fasta count = $it" }

    // Run common analysis steps
    COMMON_ANALYSIS(
        ch_sorted_bam,
        ch_fasta,
        ch_fasta_dict,
        ch_fastp_html,
        ch_fastp_json,
        params.input_type
    )
}