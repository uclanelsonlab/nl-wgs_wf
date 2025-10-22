process PICARD_MARKDUPLICATES {
    tag "$meta.id"
    label 'picard_markduplicates'

    input:
        tuple val(meta), path(bam), path(bam_index)
        tuple val(meta2), path(fasta), path(fai)
    
    output:
        tuple val(meta), path("${prefix}.marked.bam"), emit: bam
        tuple val(meta), path("${prefix}.marked.bam.metrics"), emit: metrics
        path "versions.yml", emit: versions
    
    when:
        task.ext.when == null || task.ext.when
    
    script:
        def args = task.ext.args ?: ''
        prefix = task.ext.prefix ?: "${meta.id}"
        def avail_mem = 3072
        if (!task.memory) {
            log.info '[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this default.'
        } else {
            avail_mem = (task.memory.mega*0.8).intValue()
        }
        
        """
        picard \\
            -Xmx${avail_mem}M \\
            MarkDuplicates \\
            --INPUT ${bam} \\
            --OUTPUT ${prefix}.marked.bam \\
            --METRICS_FILE ${prefix}.marked.bam.metrics \\
            --REFERENCE_SEQUENCE $fasta \\
            --TMP_DIR . \\
            $args
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            picard: \$(picard MarkDuplicates --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
        END_VERSIONS
        """
    
    stub:
        prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.marked.bam
        touch ${prefix}.marked.bam.metrics
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            picard: \$(picard MarkDuplicates --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
        END_VERSIONS
        """
}

process PICARD_COLLECT_MULTIPLE_METRICS {
    tag "$meta.id"
    label "picard_collect_multiple_metrics"
    
    input:
        tuple val(meta), path(bam), path(bam_index)
        tuple val(meta2), path(fasta_files)
    
    output:
        tuple val(meta), path("*_metrics"), emit: metrics_files
        tuple val(meta), path("*.pdf"), emit: pdf_files, optional: true
        path "picard_versions.yml", emit: versions
    
    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        def fasta = fasta_files[0]  // First file is the FASTA
        def avail_mem = 50000
        if (!task.memory) {
            log.info '[Picard CollectMultipleMetrics] Available memory not known - defaulting to 50GB. Specify process memory requirements to change this.'
        } else {
            avail_mem = (task.memory.mega*0.8).intValue()
        }
        """
        java -Xmx${avail_mem}M -jar /usr/picard/picard.jar CollectMultipleMetrics \
            -R ${fasta} \
            -I ${bam} \
            -O ${prefix} \
            --PROGRAM CollectAlignmentSummaryMetrics \
            --PROGRAM CollectInsertSizeMetrics \
            --PROGRAM QualityScoreDistribution \
            --PROGRAM MeanQualityByCycle \
            --PROGRAM CollectBaseDistributionByCycle

        cat <<-END_VERSIONS > picard_versions.yml
        "${task.process}":
            picard: \$(echo \$(picard CollectMultipleMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
        END_VERSIONS
        """
}

process PICARD_COLLECT_WGS_METRICS {
    tag "$meta.id"
    label "picard_collect_wgs_metrics"
    
    input:
        tuple val(meta), path(bam), path(bam_index)
        tuple val(meta2), path(fasta_files)
    
    output:
        tuple val(meta), path("*_metrics"), emit: wgs_metrics
        path "picard_versions.yml", emit: versions
    
    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        def fasta = fasta_files[0]  // First file is the FASTA
        def avail_mem = 50000
        if (!task.memory) {
            log.info '[Picard CollectWgsMetrics] Available memory not known - defaulting to 50GB. Specify process memory requirements to change this.'
        } else {
            avail_mem = (task.memory.mega*0.8).intValue()
        }
        """
        java -Xmx${avail_mem}M -jar /usr/picard/picard.jar CollectWgsMetrics \
            -R ${fasta} \
            -I ${bam} \
            -O ${prefix}.wgs_metrics

        cat <<-END_VERSIONS > picard_versions.yml
        "${task.process}":
            picard: \$(echo \$(picard CollectWgsMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
        END_VERSIONS
        """
}