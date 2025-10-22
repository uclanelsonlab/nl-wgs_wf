include { SAMTOOLS_SORT } from '../../../modules/samtools/main'

workflow BAM_PROCESSING {
    take:
    ch_input_bam    // channel: [meta, bam] or [meta, bam, bai]
    ch_fasta        // channel: [meta, [fasta, fai]]

    main:
    ch_versions = Channel.empty()

    // Extract BAM file (ignore index if provided)
    ch_bam_only = ch_input_bam.map { tuple ->
        def meta = tuple[0]
        def bam = tuple[1]
        [meta, bam]
    }

    // Sort BAM to ensure proper coordinate order
    SAMTOOLS_SORT(ch_bam_only)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    emit:
    sorted_bam = SAMTOOLS_SORT.out.bam     // channel: [meta, bam]
    versions   = ch_versions               // channel: [versions.yml]
}