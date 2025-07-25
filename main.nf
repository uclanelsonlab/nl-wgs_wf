nextflow.enable.dsl = 2

include { FASTP } from './modules/fastp/main'
include { CNVPYTOR } from './modules/cnvpytor/main'
include { BWAMEM2_MEM } from './modules/bwamem2/main'
include { MANTA_GERMLINE } from './modules/manta/main'
include { PICARD_MARKDUPLICATES } from './modules/picard/main'
include { DEEPVARIANT_RUNDEEPVARIANT } from './modules/deepvariant/main'
include { SAMTOOLS_SORT; SAMTOOLS_BAM2CRAM; SAMTOOLS_INDEX } from './modules/samtools/main'
include { EXPANSIONHUNTER; EXPANSIONHUNTERDENOVO_PROFILE } from './modules/expansionhunter/main'

log.info """\
    SHORTREAD - WGS _ W F   P I P E L I N E
    ===================================
    sample_name  : ${params.sample_name}
    fastq_r1     : ${params.fastq_r1}
    fastq_r2     : ${params.fastq_r2}
    fasta        : ${params.fasta}
    fai          : ${params.fai}
    """
    .stripIndent(true)

// Single channel for paired-end reads
Channel
    .fromPath(params.fastq_r1)
    .map { r1_file ->
        def meta = [:]
        meta.id = params.sample_name ?: r1_file.baseName
        
        // Use the R2 parameter directly
        def r2_file = file(params.fastq_r2)
        
        [ meta, [r1_file, r2_file] ]
    }
    .set { ch_reads }

// Reference genome channel
Channel
    .fromPath([params.fasta, params.fai])
    .collect()
    .map { files ->
        def meta = [:]
        meta.id = files[0].baseName
        [ meta, files ]
    }
    .set { ch_fasta }

// BWA index channel - collect all index files together
Channel
    .fromPath("${params.fasta}.*")
    .filter { it.name.endsWith('.amb') || it.name.endsWith('.ann') || it.name.endsWith('.bwt.2bit.64') || it.name.endsWith('.pac') || it.name.endsWith('.0123') }
    .collect()
    .map { files ->
        def meta = [:]
        meta.id = files[0].baseName
        
        // Check if we have the required BWA-MEM2 files
        def has_0123 = files.any { it.name.endsWith('.0123') }
        def has_amb = files.any { it.name.endsWith('.amb') }
        
        if (!has_amb) {
            error "Missing required BWA index file: ${params.fasta}.amb"
        }
        
        if (!has_0123) {
            log.warn "Missing BWA-MEM2 .0123 file. BWA-MEM2 may fail. Consider running: bwa-mem2 index ${params.fasta}"
        }
        [ meta, files ]
    }
    .set { ch_index }

workflow {
    // Run Fastp for quality control and adapter trimming
    FASTP(ch_reads)
    
    // Run BWA-MEM2 alignment using cleaned reads from FASTP
    BWAMEM2_MEM(
        FASTP.out.reads,
        ch_index,
        ch_fasta
    )
    // Sort BAM file
    SAMTOOLS_SORT(
        BWAMEM2_MEM.out.bam
    )
    // Mark duplicates
    PICARD_MARKDUPLICATES(
        SAMTOOLS_SORT.out.bam,
        ch_fasta
    )
    // Index BAM file
    SAMTOOLS_INDEX(
        PICARD_MARKDUPLICATES.out.bam
    )
    // Run DeepVariant
    DEEPVARIANT_RUNDEEPVARIANT(
        SAMTOOLS_INDEX.out.bam,
        ch_fasta
    )
    // Run ExpansionHunter
    EXPANSIONHUNTER(
        SAMTOOLS_INDEX.out.bam,
        ch_fasta,
        params.variant_catalog
    )
    // Run Manta
    MANTA_GERMLINE(
        SAMTOOLS_INDEX.out.bam,
        ch_fasta
    )
    // Run CNVpytor
    CNVPYTOR(
        SAMTOOLS_INDEX.out.bam,
        ch_fasta
    )
    // Run ExpansionHunterDenovo
    EXPANSIONHUNTERDENOVO_PROFILE(
        SAMTOOLS_INDEX.out.bam,
        ch_fasta,
        params.min_anchor_mapq,
        params.max_irr_mapq
    )
    // Convert BAM to CRAM
    SAMTOOLS_BAM2CRAM(
        SAMTOOLS_INDEX.out.bam,
        ch_fasta
    )
}