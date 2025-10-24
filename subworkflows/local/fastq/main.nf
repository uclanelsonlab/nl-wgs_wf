include { FASTP } from '../../../modules/fastp/main'
include { BWAMEM2_MEM } from '../../../modules/bwamem2/main'
include { SAMTOOLS_SORT } from '../../../modules/samtools/main'

workflow FASTQ_PROCESSING {
    take:
    ch_reads        // channel: [meta, [fastq_r1, fastq_r2]]
    ch_index        // channel: [meta, bwa_index_files]
    ch_fasta        // channel: [meta, [fasta, fai]]

    main:
    ch_versions = Channel.empty()

    // Debug inputs
    ch_reads.view { "FASTQ_PROCESSING: ch_reads input = $it" }
    ch_index.view { "FASTQ_PROCESSING: ch_index input = $it" }
    ch_fasta.view { "FASTQ_PROCESSING: ch_fasta input = $it" }

    // Quality control and adapter trimming
    FASTP(ch_reads)
    ch_versions = ch_versions.mix(FASTP.out.versions)

    // Alignment
    BWAMEM2_MEM(
        ch_reads,
        ch_index,
        ch_fasta
    )
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)
    BWAMEM2_MEM.out.bam.view { "FASTQ_PROCESSING: BWAMEM2_MEM.out.bam = $it" }

    // Sort BAM
    SAMTOOLS_SORT(BWAMEM2_MEM.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)
    SAMTOOLS_SORT.out.bam.view { "FASTQ_PROCESSING: SAMTOOLS_SORT.out.bam = $it" }

    emit:
    sorted_bam = SAMTOOLS_SORT.out.bam     // channel: [meta, bam]
    fastp_html = FASTP.out.html            // channel: [meta, html]
    fastp_json = FASTP.out.json            // channel: [meta, json]
    versions   = ch_versions               // channel: [versions.yml]
}