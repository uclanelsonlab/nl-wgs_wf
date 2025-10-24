include { SAMTOOLS_CRAM2BAM } from '../../../modules/samtools/main'

workflow CRAM_PROCESSING {
    take:
    ch_input_cram   // channel: [meta, cram] or [meta, cram, crai]
    ch_fasta        // channel: [meta, [fasta, fai]]

    main:
    ch_versions = Channel.empty()

    // Extract CRAM file (ignore index if provided)
    ch_cram_only = ch_input_cram.map { tuple ->
        def meta = tuple[0]
        def cram = tuple[1]
        [meta, cram]
    }

    // Prepare fasta channel with proper structure
    ch_fasta_proper = ch_fasta.map { meta, files -> 
        [meta, files[0], files[1]]  // [meta, fasta, fai]
    }

    // Convert CRAM to BAM
    SAMTOOLS_CRAM2BAM(
        ch_cram_only,
        ch_fasta_proper
    )
    ch_versions = ch_versions.mix(SAMTOOLS_CRAM2BAM.out.versions)

    emit:
    sorted_bam = SAMTOOLS_CRAM2BAM.out.bam // channel: [meta, bam]
    versions   = ch_versions               // channel: [versions.yml]
}