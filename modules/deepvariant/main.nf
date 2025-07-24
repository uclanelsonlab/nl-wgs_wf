process DEEPVARIANT_RUNDEEPVARIANT {
    tag "$meta.id"
    label 'deepvariant_run_deepvariant'

    input:
        tuple val(meta), path(bam), path(bai)
        tuple val(meta2), path(fasta_files)

    output:
        tuple val(meta), path("*.deepvariant.vcf.gz"), path("*.deepvariant.vcf.gz.tbi"), emit: vcf
        tuple val(meta), path("*.deepvariant.gvcf.gz"), path("*.deepvariant.gvcf.gz.tbi"), emit: gvcf
        path "deepvariant_versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        def fasta = fasta_files[0]  // First file is the FASTA

        """
        /opt/deepvariant/bin/run_deepvariant \\
            --model_type WGS \\
            --ref ${fasta} \\
            --reads ${bam} \\
            --output_vcf ${prefix}.deepvariant.vcf.gz \\
            --output_gvcf ${prefix}.deepvariant.gvcf.gz \\
            --num_shards ${task.cpus} --logging_dir ./logs

        cat <<-END_VERSIONS > deepvariant_versions.yml
        "${task.process}":
            deepvariant_callvariants: \$(echo \$(/opt/deepvariant/bin/run_deepvariant --version) | sed 's/^.*version //; s/ .*\$//' )
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        echo "" | gzip > ${prefix}.deepvariant.vcf.gz
        echo "" | gzip > ${prefix}.deepvariant.gvcf.gz

        cat <<-END_VERSIONS > deepvariant_versions.yml
        "${task.process}":
            deepvariant_callvariants: \$(echo \$(/opt/deepvariant/bin/run_deepvariant --version) | sed 's/^.*version //; s/ .*\$//' )
        END_VERSIONS
        """
}