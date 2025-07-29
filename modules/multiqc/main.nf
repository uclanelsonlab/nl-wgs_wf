process MULTIQC {
    tag "$meta.id"
    label "multiqc"
        
    input:
        path fastp_html 
        path fastp_json 
        tuple val(meta), path(metrics_files)
        tuple val(meta), path(wgs_metrics)
        tuple val(meta), path(qualimap_results)
    
    output:
        tuple val(meta), path("*multiqc_report.html"), emit: report
        tuple val(meta), path("*_data")              , emit: data
        tuple val(meta), path("*_plots")             , emit: plots
        path "multiqc_versions.yml"                  , emit: versions
    
    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        multiqc . \
            --filename ${prefix}_multiqc_report.html \
            --outdir .

        cat <<-END_VERSIONS > multiqc_versions.yml
        "${task.process}":
            multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        mkdir multiqc_data
        mkdir multiqc_plots
        touch multiqc_report.html

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
        END_VERSIONS
        """
}