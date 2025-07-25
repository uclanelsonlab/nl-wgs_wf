process AUTOMAP {
    tag "$meta.id"
    label 'automap'

    input:
        tuple val(meta), path(vcf), path(tbi)
        val genome
    
    output:
        tuple val(meta), path("*.cram"), path("*.crai"), emit: cram
        path  "automap_versions.yml", emit: versions

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        def VERSION = "1.3"
        """
        bgzip ${vcf}
        AutoMap_v1.3.sh --vcf ${prefix}.deepvariant.vcf --out . --genome ${genome} --id ${prefix}

        cat <<-END_VERSIONS > automap_versions.yml
        "${task.process}":
            automap: \$(echo ${VERSION}; echo \$(AutoMap_v1.3.sh) )
        END_VERSIONS
        """
}