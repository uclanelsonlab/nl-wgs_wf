process BCFTOOLS_VIEW_DIPLOID {
    tag "$meta.id"
    label 'bcftools_view_diploid'
    
    input:
    tuple val(meta), path(vcf), path(tbi)
    
    output:
    tuple val(meta), path("${prefix}.diploid.vcf"), emit: vcf
    path "bcftools_versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    bcftools view \\
        -i 'GT="0/0" || GT="0/1" || GT="1/1"' \\
        $args \\
        $vcf > ${prefix}.diploid.vcf
    
    cat <<-END_VERSIONS > bcftools_versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
    
    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.diploid.vcf
    
    cat <<-END_VERSIONS > bcftools_versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}