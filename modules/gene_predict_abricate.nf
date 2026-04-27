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
    def extra_opts = params.abricate_opt?.trim() ? params.abricate_opt : ''

    """
    mkdir -p abricate_report

    ls ${collected_dir}/*.fasta > abricate_list.txt

    ABRICATE_DBS=\$(abricate --list | awk 'NR>1 {print \$1}')

    for DB in \$ABRICATE_DBS; do
        echo "Running abricate with database: \$DB"
        abricate ${extra_opts} \\
            --db "\$DB" \\
            --fofn abricate_list.txt \\
            --nopath \\
            > "abricate_report/\${DB}.tsv"
    done

    # generate html report from tsv outputs
    python3 ${projectDir}/bin/abricate_report.py abricate_report abricate_report.html
    """
}