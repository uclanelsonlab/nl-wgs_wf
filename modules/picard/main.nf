process PICARD_MARKDUPLICATES {
    tag "$meta.id"
    label 'picard_markduplicates'

    input:
        tuple val(meta), path(reads)
        tuple val(meta2), path(fasta_files)

    output:
        tuple val(meta), path("*.MarkDuplicates.bam"), path("*.MarkDuplicates.bai") , emit: bam
        tuple val(meta), path("*.metrics.txt")                       , emit: metrics
        path  "versions.yml"                                         , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        def avail_mem = 3072
        if (!task.memory) {
            log.info '[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
        } else {
            avail_mem = (task.memory.mega*0.8).intValue()
        }
        def fasta = fasta_files[0]  // First file is the FASTA
        
        """
        java -Xmx${avail_mem}M -Djava.io.tmpdir=/temp_dir -jar /usr/picard/picard.jar \\
            MarkDuplicates \\
            --VALIDATION_STRINGENCY LENIENT \\
            --ASSUME_SORT_ORDER coordinate \\
            --PROGRAM_RECORD_ID MarkDuplicates \\
            --CREATE_INDEX true \\
            --REFERENCE_SEQUENCE ${fasta} \\
            --INPUT $reads \\
            --OUTPUT ${prefix}.MarkDuplicates.bam  \\
            --METRICS_FILE ${prefix}.MarkDuplicates.metrics.txt 

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            picard: \$(echo \$(picard MarkDuplicates --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.MarkDuplicates.bam
        touch ${prefix}.MarkDuplicates.metrics.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            picard: \$(echo \$(picard MarkDuplicates --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
        END_VERSIONS
        """
}