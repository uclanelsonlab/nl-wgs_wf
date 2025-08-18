process AUTOMAP {
    tag "$meta.id"
    label 'automap'

    input:
        tuple val(meta), path(vcf), path(tbi)
        val genome
    
    output:
        tuple val(meta), path("*/*.HomRegions.tsv"), path("*/*.HomRegions.pdf") , emit: homregions
        path  "automap_versions.yml"                                            , emit: versions

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        def VERSION = "1.3"
        """
        zcat ${vcf} > ${prefix}.vcf
        AutoMap_v1.3.sh --genome ${genome} --out . --id ${prefix} --vcf ${prefix}.vcf

        cat <<-END_VERSIONS > automap_versions.yml
        "${task.process}":
            automap: \$(echo ${VERSION}; echo \$(AutoMap_v1.3.sh) )
        END_VERSIONS
        """
}