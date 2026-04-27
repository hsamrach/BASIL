process assembly_html_report {
    label 'assembly_report'

    publishDir "${params.outdir}/Report", mode: 'copy'

    input:
    path(assembly_files)
    path(polishing_files)

    output:
    path "assembly_report.html", emit: html

    script:
    """
    python3 ${projectDir}/bin/assembly_report.py . assembly_report.html
    """
}