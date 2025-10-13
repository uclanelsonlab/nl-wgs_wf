process AUTOMAP {
    tag "$meta.id"
    label 'automap'

    input:
        tuple val(meta), path(vcf), path(tbi)
        val genome
    
    output:
        tuple val(meta), path("*.HomRegions.tsv"), emit: tsv
        tuple val(meta), path("*.HomRegions.pdf"), emit: pdf, optional: true
        path  "automap_versions.yml", emit: versions

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        def VERSION = "1.3"
        
        // Validate genome parameter
        if (!(genome in ['hg19', 'hg38'])) {
            error "Invalid genome parameter: ${genome}. Must be 'hg19' or 'hg38'"
        }
        
        """
        # Debug: Show genome parameter
        echo "DEBUG: Using genome parameter: '${genome}'"
        
        # Extract VCF if compressed
        zcat ${vcf} > ${prefix}.vcf
        
        # Run AutoMap
        AutoMap_v1.3.sh --genome ${genome} --out . --id ${prefix} --vcf ${prefix}.vcf
        
        # Move files to expected location if they're in a subdirectory
        if [ -d ${prefix} ]; then
            echo "Moving files from ${prefix}/ directory"
            mv ${prefix}/*.HomRegions.tsv . 2>/dev/null || true
            mv ${prefix}/*.HomRegions.pdf . 2>/dev/null || true
        fi
        
        # List all files for debugging
        echo "Files created:"
        find . -name "*.HomRegions.*" -type f
        
        # Ensure at least the TSV file exists
        if [ ! -f ${prefix}.HomRegions.tsv ]; then
            echo "ERROR: TSV output file not created"
            exit 1
        fi
        
        # Create a dummy PDF if it doesn't exist (to satisfy the output requirement)
        if [ ! -f ${prefix}.HomRegions.pdf ]; then
            echo "WARNING: PDF output file not created, creating dummy file"
            touch ${prefix}.HomRegions.pdf
        fi

        cat <<-END_VERSIONS > automap_versions.yml
        "${task.process}":
            automap: \$(echo ${VERSION})
        END_VERSIONS
        """
}