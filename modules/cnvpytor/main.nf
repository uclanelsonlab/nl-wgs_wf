process CNVPYTOR {
    tag "$meta.id"
    label 'sv_cnvpytor'

    input:
        tuple val(meta), path(bam), path(bai)
        tuple val(meta2), path(fasta_files)

    output:
        tuple val(meta), path("*.pytor") , emit: cnvpytor_pytor
        tuple val(meta), path("*.tsv")   , emit: cnvpytor_tab
        tuple val(meta), path("*.vcf")   , emit: cnvpytor_vcf
        path "cnvpytor_versions.yml"     , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        def fasta = fasta_files[0]  // First file is the FASTA
        """
        cnvpytor -root ${prefix}.pytor -rd ${bam} -T ${fasta}
        cnvpytor -root ${prefix}.pytor -his 1000
        cnvpytor -root ${prefix}.pytor -partition 1000
        cnvpytor -root ${prefix}.pytor -call 1000 > ${prefix}.tsv

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