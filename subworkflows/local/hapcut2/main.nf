include { BCFTOOLS_VIEW_DIPLOID } from '../../../modules/bcftools/main'
include { HAPCUT2_EXTRACTHAIRS; HAPCUT2_HAPCUT2} from '../../../modules/hapcut2/main'

workflow HAPCUT2_PHASING {
    take:
    ch_bam_bai      // channel: [meta, bam, bai]
    ch_vcf          // channel: [meta, vcf]
    ch_fasta        // channel: [meta, [fasta, fai]]

    main:
    ch_versions = Channel.empty()

    // Filter VCF for diploid genotypes
    BCFTOOLS_VIEW_DIPLOID(ch_vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW_DIPLOID.out.versions)

    // Prepare fasta channel for extractHAIRS
    ch_fasta_files = ch_fasta.map { meta, files ->
        [meta, files[0], files[1]]  // [meta, fasta, fai]
    }

    // Extract haplotype-informative reads
    HAPCUT2_EXTRACTHAIRS(
        ch_bam_bai,
        BCFTOOLS_VIEW_DIPLOID.out.vcf,
        ch_fasta_files.map { meta, fasta, fai -> [meta, fasta] },
        ch_fasta_files.map { meta, fasta, fai -> [meta, fai] }
    )
    ch_versions = ch_versions.mix(HAPCUT2_EXTRACTHAIRS.out.versions)

    // Perform haplotype phasing
    HAPCUT2_HAPCUT2(
        HAPCUT2_EXTRACTHAIRS.out.fragments,
        BCFTOOLS_VIEW_DIPLOID.out.vcf
    )
    ch_versions = ch_versions.mix(HAPCUT2_HAPCUT2.out.versions)

    emit:
    diploid_vcf = BCFTOOLS_VIEW_DIPLOID.out.vcf        // channel: [meta, vcf]
    fragments   = HAPCUT2_EXTRACTHAIRS.out.fragments   // channel: [meta, fragments]
    haplotypes  = HAPCUT2_HAPCUT2.out.haplotypes       // channel: [meta, haplotypes]
    versions    = ch_versions                          // channel: [versions.yml]
}