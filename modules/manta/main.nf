process MANTA_GERMLINE {
    tag "$meta.id"
    label 'manta_germline'

    input:
        tuple val(meta), path(bam), path(bai)
        tuple val(meta2), path(fasta_files)

    output:
        tuple val(meta), path("*candidate_small_indels.vcf.gz")    , emit: candidate_small_indels_vcf
        tuple val(meta), path("*candidate_small_indels.vcf.gz.tbi"), emit: candidate_small_indels_vcf_tbi
        tuple val(meta), path("*candidate_sv.vcf.gz")              , emit: candidate_sv_vcf
        tuple val(meta), path("*candidate_sv.vcf.gz.tbi")          , emit: candidate_sv_vcf_tbi
        tuple val(meta), path("*diploid_sv.vcf.gz")                , emit: diploid_sv_vcf
        tuple val(meta), path("*diploid_sv.vcf.gz.tbi")            , emit: diploid_sv_vcf_tbi
        path "manta_versions.yml"                                  , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        def fasta = fasta_files[0]  // First file is the FASTA
        """
        configManta.py \\
            --bam ${bam} \\
            --reference $fasta \\
            --runDir manta \\
            --referenceFasta ${fasta}

        python manta/runWorkflow.py -m local -j $task.cpus

        mv manta/results/variants/candidateSmallIndels.vcf.gz \\
            ${prefix}.candidate_small_indels.vcf.gz
        mv manta/results/variants/candidateSmallIndels.vcf.gz.tbi \\
            ${prefix}.candidate_small_indels.vcf.gz.tbi
        mv manta/results/variants/candidateSV.vcf.gz \\
            ${prefix}.candidate_sv.vcf.gz
        mv manta/results/variants/candidateSV.vcf.gz.tbi \\
            ${prefix}.candidate_sv.vcf.gz.tbi
        mv manta/results/variants/diploidSV.vcf.gz \\
            ${prefix}.diploid_sv.vcf.gz
        mv manta/results/variants/diploidSV.vcf.gz.tbi \\
            ${prefix}.diploid_sv.vcf.gz.tbi

        cat <<-END_VERSIONS > manta_versions.yml
        "${task.process}":
            manta: \$( configManta.py --version )
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        echo "" | gzip > ${prefix}.candidate_small_indels.vcf.gz
        touch ${prefix}.candidate_small_indels.vcf.gz.tbi
        echo "" | gzip > ${prefix}.candidate_sv.vcf.gz
        touch ${prefix}.candidate_sv.vcf.gz.tbi
        echo "" | gzip > ${prefix}.diploid_sv.vcf.gz
        touch ${prefix}.diploid_sv.vcf.gz.tbi

        cat <<-END_VERSIONS > manta_versions.yml
        "${task.process}":
            manta: \$( configManta.py --version )
        END_VERSIONS
        """
}