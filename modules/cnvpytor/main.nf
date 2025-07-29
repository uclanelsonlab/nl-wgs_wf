process CNVPYTOR {
    tag "$meta.id"
    label 'sv_cnvpytor'

    input:
        tuple val(meta), path(bam), path(bai)
        tuple val(meta2), path(fasta_files)

    output:
        tuple val(meta), path("*.pytor")            , emit: cnvpytor_pytor
        tuple val(meta), path("*1000.vcf")          , emit: cnvpytor_vcf
        tuple val(meta), path("*filtered.vcf")      , emit: cnvpytor_filtered_vcf
        tuple val(meta), path("*manhattan_plot.png"), emit: cnvpytor_manhattan_plot
        path "cnvpytor_versions.yml"                , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        def fasta = fasta_files[0]  // First file is the FASTA
        """
        cnvpytor -root ${prefix}.pytor -rd ${bam} -T ${fasta}
        cnvpytor -root ${prefix}.pytor -his 1000
        cnvpytor -root ${prefix}.pytor -partition 1000
        cnvpytor -root ${prefix}.pytor -view 1000 <<EOF
        set Q0_range -1 0.5
        set p_N 0 0.5
        set size_range 300 inf
        set dG_range 100000 inf
        set p_range 0 0.05
        set print_filename ${prefix}_cnvpytor_1000_filtered.vcf
        print calls
        EOF

        cnvpytor -root ${prefix}.pytor -view 1000 <<EOF 
        set style bmh
        set rd_use_mask
        set file_titles ${prefix}
        manhattan
        save ${prefix}_manhattan_plot.png
        EOF

        python3 <<CODE
        import cnvpytor,os
        binsizes = "1000".split(" ")
        for binsize in binsizes:
            file_list = "${prefix}.pytor".split(" ")
            app = cnvpytor.Viewer(file_list, params={} )
            outputfile = "{}_{}.{}".format("${prefix}_cnvpytor",binsize.strip(),"vcf")
            app.print_filename = outputfile
            app.bin_size = int(binsize)
            app.print_calls_file()
        CODE

        cat <<-END_VERSIONS > cnvpytor_versions.yml
        "${task.process}":
            cnvpytor: \$(cnvpytor --version | sed -n 's/.*CNVpytor \\(.*\\)/\\1/p')
        END_VERSIONS
        """
}