nextflow.enable.dsl = 2

include { FASTP } from './modules/fastp/main'
include { AUTOMAP } from './modules/automap/main'
include { MULTIQC } from './modules/multiqc/main'
include { CNVPYTOR } from './modules/cnvpytor/main'
include { BWAMEM2_MEM } from './modules/bwamem2/main'
include { MANTA_GERMLINE } from './modules/manta/main'
include { QUALIMAP_BAMQC } from './modules/qualimap/main'
include { MOSDEPTH_BED } from './modules/mosdepth/main'
include { DEEPVARIANT_RUNDEEPVARIANT } from './modules/deepvariant/main'
include { SAMTOOLS_SORT; SAMTOOLS_BAM2CRAM; SAMTOOLS_INDEX; SAMTOOLS_CRAM2BAM } from './modules/samtools/main'
include { EXPANSIONHUNTER; EXPANSIONHUNTERDENOVO_PROFILE } from './modules/expansionhunter/main'
include { PICARD_MARKDUPLICATES; PICARD_COLLECT_MULTIPLE_METRICS; PICARD_COLLECT_WGS_METRICS } from './modules/picard/main'

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

// Input channels based on input type
if (params.input_type == "fastq") {
    // Single channel for paired-end reads
    Channel
        .fromPath(params.input_1)  // fastq_r1
        .map { r1_file ->
            def meta = [:]
            meta.id = params.sample_name ?: r1_file.baseName
            
            // Use input_2 as fastq_r2
            def r2_file = file(params.input_2)
            
            [ meta, [r1_file, r2_file] ]
        }
        .set { ch_reads }
    
    // Empty channels for BAM/CRAM inputs
    ch_input_bam = Channel.empty()
    ch_input_cram = Channel.empty()
} else if (params.input_type == "bam") {
    // BAM input channel with optional index
    if (params.input_2) {
        // If index is provided (input_2 = bam_index), use it
        Channel
            .fromPath([params.input_1, params.input_2])  // [bam_input, bam_index]
            .collect()
            .map { files ->
                def meta = [:]
                meta.id = params.sample_name ?: files[0].baseName
                [ meta, files[0], files[1] ]  // [meta, bam, bai]
            }
            .set { ch_input_bam }
    } else {
        // If no index provided, just the BAM file (index will be created later)
        Channel
            .fromPath(params.input_1)  // bam_input
            .map { bam_file ->
                def meta = [:]
                meta.id = params.sample_name ?: bam_file.baseName
                [ meta, bam_file ]
            }
            .set { ch_input_bam }
    }
    
    // Empty channels for FASTQ/CRAM inputs
    ch_reads = Channel.empty()
    ch_input_cram = Channel.empty()
} else if (params.input_type == "cram") {
    // CRAM input channel with optional index
    if (params.input_2) {
        // If index is provided (input_2 = cram_index), use it
        Channel
            .fromPath([params.input_1, params.input_2])  // [cram_input, cram_index]
            .collect()
            .map { files ->
                def meta = [:]
                meta.id = params.sample_name ?: files[0].baseName
                [ meta, files[0], files[1] ]  // [meta, cram, crai]
            }
            .set { ch_input_cram }
    } else {
        // If no index provided, just the CRAM file
        Channel
            .fromPath(params.input_1)  // cram_input
            .map { cram_file ->
                def meta = [:]
                meta.id = params.sample_name ?: cram_file.baseName
                [ meta, cram_file ]
            }
            .set { ch_input_cram }
    }
    
    // Empty channels for FASTQ/BAM inputs
    ch_reads = Channel.empty()
    ch_input_bam = Channel.empty()
} else {
    error "Invalid input_type: ${params.input_type}. Must be 'fastq', 'bam', or 'cram'"
}

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

// Reference genome channel with dict for mosdepth
Channel
    .fromPath([params.fasta, params.fai, params.dict])
    .collect()
    .map { files ->
        def meta = [:]
        meta.id = files[0].baseName
        [ meta, files[0], files[1], files[2] ]  // [meta, fasta, fai, dict]
    }
    .set { ch_fasta_dict }

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
    // Initialize channels for different workflow paths
    ch_sorted_bam = Channel.empty()
    ch_fastp_html = Channel.empty()
    ch_fastp_json = Channel.empty()
    
    if (params.input_type == "fastq") {
        // Run Fastp for quality control and adapter trimming
        FASTP(ch_reads)
        
        // Run BWA-MEM2 alignment using cleaned reads from FASTP
        BWAMEM2_MEM(
            ch_reads,
            ch_index,
            ch_fasta
        )
        // Sort BAM file
        SAMTOOLS_SORT(
            BWAMEM2_MEM.out.bam
        )
        
        ch_sorted_bam = SAMTOOLS_SORT.out.bam
        ch_fastp_html = FASTP.out.html
        ch_fastp_json = FASTP.out.json
        
    } else if (params.input_type == "cram") {
        // Convert CRAM to BAM first
        if (params.input_2) {
            // If CRAM has index, we can use it directly for conversion
            SAMTOOLS_CRAM2BAM(
                ch_input_cram.map { meta, cram, crai -> [meta, cram] },
                ch_fasta
            )
        } else {
            // No index provided, just convert
            SAMTOOLS_CRAM2BAM(
                ch_input_cram,
                ch_fasta
            )
        }
        
        ch_sorted_bam = SAMTOOLS_CRAM2BAM.out.bam.map { meta, bam ->
            // Create a tuple with bam and empty index (will be created later)
            [ meta, bam, [] ]
        }
        
    } else if (params.input_type == "bam") {
        if (params.input_2) {
            // If BAM has index, check if it needs sorting or can go directly to mark duplicates
            // For now, we'll still sort to ensure proper coordinate order
            SAMTOOLS_SORT(
                ch_input_bam.map { meta, bam, bai -> [meta, bam] }
            )
            ch_sorted_bam = SAMTOOLS_SORT.out.bam
        } else {
            // No index provided, sort the BAM
            SAMTOOLS_SORT(
                ch_input_bam
            )
            ch_sorted_bam = SAMTOOLS_SORT.out.bam
        }
    }
    
    // Mark duplicates (common step for all input types)
    PICARD_MARKDUPLICATES(
        ch_sorted_bam,
        ch_fasta
    )
    
    // Index BAM file
    SAMTOOLS_INDEX(
        PICARD_MARKDUPLICATES.out.bam
    )
    
    // Run BAM QC to collect metrics
    PICARD_COLLECT_MULTIPLE_METRICS(
        SAMTOOLS_INDEX.out.bam,
        ch_fasta
    )
    PICARD_COLLECT_WGS_METRICS(
        SAMTOOLS_INDEX.out.bam,
        ch_fasta
    )
    QUALIMAP_BAMQC(
        SAMTOOLS_INDEX.out.bam
    )   
    
    // Convert BAM to CRAM for mosdepth (mosdepth works better with CRAM)
    SAMTOOLS_BAM2CRAM(
        SAMTOOLS_INDEX.out.bam,
        ch_fasta
    )
    
    // Run mosdepth for coverage analysis
    MOSDEPTH_BED(
        ch_fasta_dict,
        SAMTOOLS_BAM2CRAM.out.cram,
        params.mt_bed
    )
    
    // Run variant calling first (needed for MultiQC)
    DEEPVARIANT_RUNDEEPVARIANT(
        SAMTOOLS_INDEX.out.bam,
        ch_fasta
    )
    
    // MultiQC - handle optional FASTP outputs
    if (params.input_type == "fastq") {
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
        // For BAM/CRAM input, create empty channels for FASTP outputs
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
    
    // Continue with other variant calling tools
    AUTOMAP(
        DEEPVARIANT_RUNDEEPVARIANT.out.vcf,
        params.genome
    )
    
    // Run repeats calling
    EXPANSIONHUNTER(
        SAMTOOLS_INDEX.out.bam,
        ch_fasta,
        params.variant_catalog
    )
    EXPANSIONHUNTERDENOVO_PROFILE(
        SAMTOOLS_INDEX.out.bam,
        ch_fasta,
        params.min_anchor_mapq,
        params.max_irr_mapq
    )
    
    // Run SV/CNV calling
    MANTA_GERMLINE(
        SAMTOOLS_INDEX.out.bam,
        ch_fasta
    )
    CNVPYTOR(
        SAMTOOLS_INDEX.out.bam,
        ch_fasta
    )  
}