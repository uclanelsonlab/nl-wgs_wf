include { PICARD_MARKDUPLICATES; PICARD_COLLECT_MULTIPLE_METRICS; PICARD_COLLECT_WGS_METRICS } from '../../../modules/picard/main'
include { SAMTOOLS_INDEX; SAMTOOLS_BAM2CRAM } from '../../../modules/samtools/main'
include { QUALIMAP_BAMQC } from '../../../modules/qualimap/main'
include { MOSDEPTH_BED } from '../../../modules/mosdepth/main'
include { DEEPVARIANT_RUNDEEPVARIANT } from '../../../modules/deepvariant/main'
include { MULTIQC } from '../../../modules/multiqc/main'
include { AUTOMAP } from '../../../modules/automap/main'
include { EXPANSIONHUNTER; EXPANSIONHUNTERDENOVO_PROFILE } from '../../../modules/expansionhunter/main'
include { MANTA_GERMLINE } from '../../../modules/manta/main'
include { CNVPYTOR } from '../../../modules/cnvpytor/main'

workflow COMMON_ANALYSIS {
    take:
    ch_sorted_bam   // channel: [meta, bam]
    ch_fasta        // channel: [meta, [fasta, fai]]
    ch_fasta_dict   // channel: [meta, fasta, fai, dict]
    ch_fastp_html   // channel: [meta, html] (may be empty)
    ch_fastp_json   // channel: [meta, json] (may be empty)
    input_type      // val: string ("fastq", "bam", or "cram")

    main:
    ch_versions = Channel.empty()

    // Mark duplicates
    PICARD_MARKDUPLICATES(
        ch_sorted_bam,
        ch_fasta
    )
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)

    // Index BAM
    SAMTOOLS_INDEX(PICARD_MARKDUPLICATES.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    // QC metrics
    PICARD_COLLECT_MULTIPLE_METRICS(
        SAMTOOLS_INDEX.out.bam,
        ch_fasta
    )
    ch_versions = ch_versions.mix(PICARD_COLLECT_MULTIPLE_METRICS.out.versions)

    PICARD_COLLECT_WGS_METRICS(
        SAMTOOLS_INDEX.out.bam,
        ch_fasta
    )
    ch_versions = ch_versions.mix(PICARD_COLLECT_WGS_METRICS.out.versions)

    QUALIMAP_BAMQC(SAMTOOLS_INDEX.out.bam)
    ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions)

    // Convert to CRAM for mosdepth
    SAMTOOLS_BAM2CRAM(
        SAMTOOLS_INDEX.out.bam,
        ch_fasta
    )
    ch_versions = ch_versions.mix(SAMTOOLS_BAM2CRAM.out.versions)

    // Coverage analysis
    MOSDEPTH_BED(
        ch_fasta_dict,
        SAMTOOLS_BAM2CRAM.out.cram,
        params.mt_bed
    )
    ch_versions = ch_versions.mix(MOSDEPTH_BED.out.versions)

    // Variant calling
    DEEPVARIANT_RUNDEEPVARIANT(
        SAMTOOLS_INDEX.out.bam,
        ch_fasta
    )
    ch_versions = ch_versions.mix(DEEPVARIANT_RUNDEEPVARIANT.out.versions)

    // MultiQC - conditional on input type
    if (input_type == "fastq") {
        MULTIQC(
            params.ref_qc,
            ch_fastp_html,
            ch_fastp_json,
            PICARD_COLLECT_MULTIPLE_METRICS.out.metrics_files,
            PICARD_COLLECT_WGS_METRICS.out.wgs_metrics,
            QUALIMAP_BAMQC.out.results,
            DEEPVARIANT_RUNDEEPVARIANT.out.report,
            MOSDEPTH_BED.out.summary
        )
    } else {
        MULTIQC(
            params.ref_qc,
            Channel.empty(),
            Channel.empty(),
            PICARD_COLLECT_MULTIPLE_METRICS.out.metrics_files,
            PICARD_COLLECT_WGS_METRICS.out.wgs_metrics,
            QUALIMAP_BAMQC.out.results,
            DEEPVARIANT_RUNDEEPVARIANT.out.report,
            MOSDEPTH_BED.out.summary
        )
    }
    ch_versions = ch_versions.mix(MULTIQC.out.versions)

    // Additional variant calling
    AUTOMAP(
        DEEPVARIANT_RUNDEEPVARIANT.out.vcf,
        params.genome
    )
    ch_versions = ch_versions.mix(AUTOMAP.out.versions)

    // Repeat analysis
    EXPANSIONHUNTER(
        SAMTOOLS_INDEX.out.bam,
        ch_fasta,
        params.variant_catalog
    )
    ch_versions = ch_versions.mix(EXPANSIONHUNTER.out.versions)

    EXPANSIONHUNTERDENOVO_PROFILE(
        SAMTOOLS_INDEX.out.bam,
        ch_fasta,
        params.min_anchor_mapq,
        params.max_irr_mapq
    )
    ch_versions = ch_versions.mix(EXPANSIONHUNTERDENOVO_PROFILE.out.versions)

    // Structural variant calling
    MANTA_GERMLINE(
        SAMTOOLS_INDEX.out.bam,
        ch_fasta
    )
    ch_versions = ch_versions.mix(MANTA_GERMLINE.out.versions)

    CNVPYTOR(
        SAMTOOLS_INDEX.out.bam,
        ch_fasta
    )
    ch_versions = ch_versions.mix(CNVPYTOR.out.versions)

    emit:
    multiqc_report = MULTIQC.out.report     // channel: [meta, html]
    vcf           = DEEPVARIANT_RUNDEEPVARIANT.out.vcf  // channel: [meta, vcf]
    versions      = ch_versions             // channel: [versions.yml]
}