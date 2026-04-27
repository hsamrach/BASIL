process main_report {
    label 'main_report'
    
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path "*"

    output:
    path "Report/*.html"

    script:
    """
    multiqc . \\
        -o Report \\
        --title "Main_Report" \\
        --filename "main_report.html" \\
        --force \\
        --cl-config '{"max_table_rows": 100000}' \\
        --exclude samtools
    """
}
