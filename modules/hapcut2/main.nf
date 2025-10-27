process HAPCUT2_EXTRACTHAIRS {
    tag "$meta.id"
    label 'hapcut2'
    
    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(vcf), path(tbi)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fai)
    
    output:
    tuple val(meta), path("${prefix}.fragment_file"), emit: fragments
    path "hapcut2_versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    extractHAIRS \\
        --ref $fasta \\
        --bam $bam \\
        --VCF $vcf \\
        --out ${prefix}.fragment_file \\
        $args
    
    cat <<-END_VERSIONS > hapcut2_versions.yml
    "${task.process}":
        hapcut2: \$(extractHAIRS --help 2>&1 | grep -i version | head -n1 | sed 's/.*version //i')
    END_VERSIONS
    """
    
    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fragment_file
    
    cat <<-END_VERSIONS > hapcut2_versions.yml
    "${task.process}":
        hapcut2: \$(extractHAIRS --help 2>&1 | grep -i version | head -n1 | sed 's/.*version //i')
    END_VERSIONS
    """
}

process HAPCUT2_HAPCUT2 {
    tag "$meta.id"
    label 'hapcut2'
    
    input:
    tuple val(meta), path(fragments)
    tuple val(meta2), path(vcf), path(tbi)
    
    output:
    tuple val(meta), path("${prefix}.haplotype_output_file*"), emit: haplotypes
    path "hapcut2_versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    HAPCUT2 \\
        --fragments $fragments \\
        --VCF $vcf \\
        --output ${prefix}.haplotype_output_file \\
        $args
    
    cat <<-END_VERSIONS > hapcut2_versions.yml
    "${task.process}":
        hapcut2: \$(HAPCUT2 --help 2>&1 | grep -i version | head -n1 | sed 's/.*version //i')
    END_VERSIONS
    """
    
    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.haplotype_output_file
    
    cat <<-END_VERSIONS > hapcut2_versions.yml
    "${task.process}":
        hapcut2: \$(HAPCUT2 --help 2>&1 | grep -i version | head -n1 | sed 's/.*version //i')
    END_VERSIONS
    """
}