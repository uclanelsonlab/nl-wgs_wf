process MOSDEPTH_BED {
    label "mosdepth_bed"

    input:
        tuple val(meta), path(fasta), path(fai), path(dict)  
        tuple val(meta), path(cram), path(crai)
        path mt_bed


    output:
        tuple val(meta), path("*.mosdepth.global.dist.txt"), emit: global_dist
        tuple val(meta), path("*.mosdepth.region.dist.txt"), emit: region_dist
        tuple val(meta), path("*.mosdepth.summary.txt"),     emit: summary
        tuple val(meta), path("*.per-base.bed.gz"),          emit: perbase
        tuple val(meta), path("*.per-base.bed.gz.csi"),      emit: perbase_index
        tuple val(meta), path("*.regions.bed.gz"),           emit: regions_bed
        tuple val(meta), path("*.regions.bed.gz.csi"),       emit: regions_bed_index
        path "mosdepth_versions.yml",      emit: versions
    
    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        ls ${crai}
        mosdepth -t ${task.cpus} -n --by ${mt_bed} -f ${fasta} ${prefix}_MT_genes_hg38 ${cram}
        mosdepth -t ${task.cpus} -n --fast-mode --by 500 -f ${fasta} ${prefix}_genome_hg38 ${cram}

        cat <<-END_VERSIONS > mosdepth_versions.yml
        "${task.process}":
            mosdepth: \$(echo \$(mosdepth --version 2>&1) | awk '{print \$2}' )
        END_VERSIONS
        """
}