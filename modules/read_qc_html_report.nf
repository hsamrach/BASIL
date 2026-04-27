process read_qc_html_report {
    label 'read_qc_report'
    
    publishDir "${params.outdir}/Report", mode: 'copy'

    input:
    path(qc_files)
    path(status_files)

    output:
    path "read_qc_summary_report.html", emit: qc_html

    script:
    """
    python3 ${projectDir}/bin/read_qc_summary_report.py . read_qc_summary_report.html
    """
}