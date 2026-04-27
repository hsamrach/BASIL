process speciation_gtdbtk {
    label 'gtdbtk'
    publishDir "${params.outdir}/speciation/", mode: 'copy'

    input:
    path collected_dir
    path gtdbtk_db

    output:
    path "gtdbtk_output", emit: gtdbtk_reports

    script:
    """
        export GTDBTK_DATA_PATH=${params.gtdbtk_db}

        gtdbtk classify_wf \\
            --genome_dir "${collected_dir}" \\
            --out_dir gtdbtk_output \\
            --extension fasta \\
            --cpus ${params.cpus} \\
            --pplacer_cpus ${params.cpus} \\
            --skip_ani_screen \\
            --force
    """
}
