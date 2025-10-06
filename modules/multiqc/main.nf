process MULTIQC {
    tag "$meta.id"
    label "multiqc"
        
    input:
        path ref_qc
        path fastp_html 
        path fastp_json 
        tuple val(meta), path(metrics_files)
        tuple val(meta2), path(wgs_metrics)
        tuple val(meta3), path(qualimap_results)
        tuple val(meta4), path(deepvariant_report)
    
    output:
        tuple val(meta), path("*multiqc_report.html"), emit: report
        tuple val(meta), path("*_data")              , emit: data
        path "multiqc_versions.yml"                  , emit: versions
    
    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        # Debug: List all files in current directory
        echo "Files in current directory:"
        ls -la
        # Unzip the reference QC files
        tar -xvzf ${ref_qc}
        
        echo "Input files:"
        echo "FASTP HTML: ${fastp_html}"
        echo "FASTP JSON: ${fastp_json}"
        echo "Metrics files: ${metrics_files}"
        echo "WGS metrics: ${wgs_metrics}"
        echo "Qualimap results: ${qualimap_results}"
        echo "DeepVariant report: ${deepvariant_report}"
        
        multiqc . \\
            --filename ${prefix}_multiqc_report.html \\
            --outdir . \\
            --verbose

        cat <<-END_VERSIONS > multiqc_versions.yml
        "${task.process}":
            multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        mkdir multiqc_data
        touch ${prefix}_multiqc_report.html

        cat <<-END_VERSIONS > multiqc_versions.yml
        "${task.process}":
            multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
        END_VERSIONS
        """
}