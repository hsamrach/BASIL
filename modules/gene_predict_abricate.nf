process gene_predict_abricate {
    label 'abricate'

    publishDir "${params.outdir}/gene_prediction", mode: 'copy', pattern: "abricate_report/*.tsv"
    publishDir "${params.outdir}/Report", mode: 'copy', pattern: "abricate_report.html"

    input:
    path collected_dir

    output:
    path "abricate_report/*.tsv"
    path "abricate_report.html"

    script:

    """
    mkdir -p abricate_report

    ls ${collected_dir}/*.fasta > abricate_list.txt

    ABRICATE_DBS=\$(abricate --list | awk 'NR>1 {print \$1}')

    for DB in \$ABRICATE_DBS; do
        echo "Running abricate with database: \$DB"
        abricate \\
            --threads ${params.cpus} \\
            --db "\$DB" \\
            --fofn abricate_list.txt \\
            --nopath \\
            --minid ${params.minid_genes} \\
            --mincov ${params.mincov_genes} \\
            > "abricate_report/\${DB}.tsv"
    done

    # generate html report from tsv outputs
    python3 ${projectDir}/python/abricate_report.py abricate_report abricate_report.html ${params.minid_genes} ${params.mincov_genes}
    """
}
