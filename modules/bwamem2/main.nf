process BWAMEM2_MEM {
    tag "$meta.id"
    label 'bwamem2_mem'

    input:
        tuple val(meta), path(reads)
        tuple val(meta2), path(index)
        tuple val(meta3), path(fasta_files)
        
    output:
        tuple val(meta), path("*.bam"), emit: bam
        path  "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        def total_mem = task.memory.toGiga() as Double
        def mem_per_thread = Math.max(1.0, (total_mem * 0.6) / (task.cpus as Double)) as Integer
        def fasta = fasta_files[0]  // First file is the FASTA

        """
        INDEX=\$(find -L ./ -name "*.amb" | sed 's/\\.amb\$//')

        bwa-mem2 mem \\
            -t ${task.cpus} \\
            -M \\
            -R "@RG\\tID:${prefix}\\tCN:CN\\tLB:${prefix}\\tPL:ILLUMINA\\tPU:${prefix}\\tSM:${prefix}" \\
            \$INDEX \\
            ${reads} | \\
        samtools view \\
            -@ ${task.cpus} \\
            -S -b > ${prefix}.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.sorted.bam
        touch ${prefix}.sorted.bam.bai

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
}