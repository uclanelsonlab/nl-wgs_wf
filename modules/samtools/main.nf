process SAMTOOLS_SORT {
    tag "$meta.id"
    label 'samtools_sort'

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: bam
        path "samtools_versions.yml", emit: versions
    
    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        samtools sort -@ ${task.cpus} -o ${prefix}.sorted.bam ${bam}
        samtools index ${prefix}.sorted.bam

        cat <<-END_VERSIONS > samtools_versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    
    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.sorted.bam
        touch ${prefix}.sorted.bam.bai
        
        cat <<-END_VERSIONS > samtools_versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
}

process SAMTOOLS_BAM2CRAM {
    tag "$meta.id"
    label 'samtools_bam2cram'

    input:
        tuple val(meta), path(bam), path(bai)
        tuple val(meta2), path(fasta_files)

    output:
        tuple val(meta), path("*.cram"), path("*.crai"), emit: cram
        path  "samtools_versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        def fasta = fasta_files[0]  // First file is the FASTA
        """ 
        samtools view -@ ${task.cpus} -T ${fasta} -C --output-fmt-option normal -o ${prefix}.cram ${bam}
        samtools index -@ ${task.cpus} ${prefix}.cram

        cat <<-END_VERSIONS > samtools_versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.cram
        touch ${prefix}.cram.crai

        cat <<-END_VERSIONS > samtools_versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
}

process SAMTOOLS_INDEX {
    tag "$meta.id"
    label 'samtools_index'

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path(bam), path("*.bai")   , emit: bam
        path "samtools_versions.yml"                    , emit: versions
    
    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        samtools index ${bam}

        cat <<-END_VERSIONS > samtools_versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    
    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.bam.bai
        
        cat <<-END_VERSIONS > samtools_versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
}