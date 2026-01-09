process STRLING_EXTRACT {
    tag "$meta.id"
    label 'strling_extract'

    input:
        tuple val(meta), path(aln), path(index)
        tuple val(meta2), path(fasta_files)
        path(str_index)

    output:
        tuple val(meta), path("*.bin"), emit: bin
        path "strling_versions.yml"   , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        def fasta = fasta_files[0]  // First file is the FASTA
        """
        strling extract \\
            -f $fasta \\
            -g $str_index \\
            $args \\
            $aln \\
            ${prefix}.bin

        cat <<-END_VERSIONS > strling_versions.yml
        "${task.process}":
            strling: \$(strling --version 2>&1 | sed 's/^.*strling //; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.bin

        cat <<-END_VERSIONS > strling_versions.yml
        "${task.process}":
            strling: \$(strling --version 2>&1 | sed 's/^.*strling //; s/ .*\$//')
        END_VERSIONS
        """
}

process STRLING_CALL {
    tag "$meta.id"
    label 'strling_call'

    input:
        tuple val(meta), path(aln), path(index)
        tuple val(meta2), path(fasta_files)
        tuple val(meta3), path(bin)

    output:
        tuple val(meta), path("*-results/*-bounds.txt")  , emit: bounds
        tuple val(meta), path("*-results/*-genotype.txt"), emit: genotype
        tuple val(meta), path("*-results/*-unplaced.txt"), emit: unplaced
        path "strling_versions.yml"                      , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        def fasta = fasta_files[0]  // First file is the FASTA
        """
        mkdir -p ${prefix}-results/

        strling call \\
            --output-prefix ${prefix}-results/${prefix} \\
            -f $fasta \\
            $args \\
            $aln \\
            $bin

        cat <<-END_VERSIONS > strling_versions.yml
        "${task.process}":
            strling: \$(strling --version 2>&1 | sed 's/^.*strling //; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        mkdir -p ${prefix}-results/
        touch ${prefix}-results/${prefix}-bounds.txt
        touch ${prefix}-results/${prefix}-genotype.txt
        touch ${prefix}-results/${prefix}-unplaced.txt

        cat <<-END_VERSIONS > strling_versions.yml
        "${task.process}":
            strling: \$(strling --version 2>&1 | sed 's/^.*strling //; s/ .*\$//')
        END_VERSIONS
        """
}