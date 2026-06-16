process gene_predict_amr_consensus {

    label 'amr_consensus'

    publishDir "${params.outdir}/Report", mode: 'copy', pattern: "amr_consensus_summary.tsv"
    publishDir "${params.outdir}/Report", mode: 'copy', pattern: "amr_consensus_report.html"

    input:
    path(abricate_tsv_dir)
    path(abritamr_summary)

    output:
    path("amr_consensus_report.html"),   emit: html_report
    path("amr_consensus_summary.tsv"),   emit: tsv_summary

    script:
    """

    if [ -d "${abricate_tsv_dir}" ]; then
        ABRICATE_DIR="${abricate_tsv_dir}"
    else
        mkdir -p abricate_tsvs
        cp ${abricate_tsv_dir} abricate_tsvs/ 2>/dev/null || true
        ABRICATE_DIR="abricate_tsvs"
    fi
    python3 ${projectDir}/python/amr_consensus_report.py \\
        --abricate_dir "\${ABRICATE_DIR}" \\
        --abritamr     "${abritamr_summary}" \\
        --output        amr_consensus_report.html \\
        --output_tsv    amr_consensus_summary.tsv \\
        --minid_genes  ${params.minid_genes} \\
        --mincov_genes ${params.mincov_genes}
    """
}
