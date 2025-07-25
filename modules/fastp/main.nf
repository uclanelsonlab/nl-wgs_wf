process FASTP {
    tag "$meta.id"
    label "fastp"

    input:
        tuple val(meta), path(reads)
    
    output:
        path "*.html", emit: html
        path "*.json", emit: json
        path "versions.yml", emit: versions
    
    when:
        task.ext.when == null || task.ext.when
    
    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        def fastq_r1 = reads[0]
        def fastq_r2 = reads[1]
        """
        fastp \\
            -i ${fastq_r1} \\
            -I ${fastq_r2} \\
            --thread ${task.cpus} \\
            -h ${prefix}_fastp.html \\
            -j ${prefix}_fastp.json \\
            --detect_adapter_for_pe \\
            --compression 4

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed 's/^fastp //')
        END_VERSIONS
        """
    
    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}_fastp.html
        touch ${prefix}_fastp.json
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed 's/^fastp //')
        END_VERSIONS
        """
}