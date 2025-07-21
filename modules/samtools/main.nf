process SAMTOOLS_SORT {
    tag "$meta.id"
    label 'samtools_sort'

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: bam
        path "versions.yml", emit: versions
    
    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        samtools sort -@ ${task.cpus} -o ${prefix}.sorted.bam ${bam}
        samtools index ${prefix}.sorted.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
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
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
}