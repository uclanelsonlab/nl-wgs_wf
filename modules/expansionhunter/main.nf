process EXPANSIONHUNTER {
    tag "$meta.id"
    label 'expansionhunter'

    input:
        tuple val(meta), path(bam), path(bai)
        tuple val(meta2), path(fasta_files)
        path(variant_catalog)

    output:
        tuple val(meta), path("*.vcf.gz")        , emit: vcf
        tuple val(meta), path("*.json.gz")       , emit: json
        tuple val(meta), path("*_realigned.bam") , emit: bam
        path "eh_versions.yml"                   , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        def fasta = fasta_files[0]  // First file is the FASTA

        """
        ls ${bai}
        ExpansionHunter \\
            --threads ${task.cpus} \\
            --reads ${bam} \\
            --output-prefix ${prefix} \\
            --reference ${fasta} \\
            --variant-catalog ${variant_catalog}

        bgzip --threads ${task.cpus} ${prefix}.vcf
        bgzip --threads ${task.cpus} ${prefix}.json

        cat <<-END_VERSIONS > eh_versions.yml
        "${task.process}":
            expansionhunter: \$( echo \$(ExpansionHunter --version 2>&1) | head -1 | sed 's/^.*ExpansionHunter v//')
            bgzip: \$(echo \$(bgzip -h 2>&1) | sed 's/^.*Version: //;s/Usage:.*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        echo "" | gzip > ${prefix}.vcf.gz
        echo "" | gzip > ${prefix}.json.gz
        touch ${prefix}_realigned.bam

        cat <<-END_VERSIONS > eh_versions.yml
        "${task.process}":
            expansionhunter: \$( echo \$(ExpansionHunter --version 2>&1) | head -1 | sed 's/^.*ExpansionHunter v//')
            bgzip: \$(echo \$(bgzip -h 2>&1) | sed 's/^.*Version: //;s/Usage:.*//')
        END_VERSIONS
        """
}

process EXPANSIONHUNTERDENOVO_PROFILE {
    tag "$meta.id"
    label 'expansionhunterdenovo_profile'

    input:
        tuple val(meta), path(bam), path(bai)
        tuple val(meta2), path(fasta_files)
        val(min_anchor_mapq)
        val(max_irr_mapq)

    output:
        tuple val(meta), path("*.locus.tsv")        , emit: locus_tsv
        tuple val(meta), path("*.motif.tsv")        , emit: motif_tsv
        tuple val(meta), path("*.str_profile.json") , emit: str_profile
        path "ehd_versions.yml"                     , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        def fasta = fasta_files[0]  // First file is the FASTA

        """
        ExpansionHunterDenovo profile \\
            --reads ${bam} \\
            --reference ${fasta} \\
            --output-prefix ${prefix} \\
            --min-anchor-mapq ${min_anchor_mapq} \\
            --max-irr-mapq ${max_irr_mapq}

        cat <<-END_VERSIONS > ehd_versions.yml
        "${task.process}":
            expansionhunterdenovo: \$(echo \$(ExpansionHunterDenovo --help 2>&1) | sed -e "s/ExpansionHunter Denovo v//;s/ Usage.*//")
        END_VERSIONS
        """
}