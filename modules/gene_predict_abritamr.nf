process gene_predict_abritamr {
    label 'abritamr'
    
    publishDir "${params.outdir}/gene_prediction/abritamr", mode: 'copy', pattern: "*.txt"
    publishDir "${params.outdir}/Report", mode: 'copy', pattern: "abritamr_report.html"

    input:
    path collected_dir

    output:
    path "summary_partials.txt"
    path "summary_virulence.txt"
    path "summary_matches.txt"
    path "abritamr_report.html"


    script:
    // only use options if it's a non-empty string
    def extra_opts = ''
    if (params.abritamr_opt && params.abritamr_opt instanceof String && params.abritamr_opt.trim()) {
        extra_opts = params.abritamr_opt
    }
    
    """
    # create tab-delimited file for abritamr input
    for fasta in ${collected_dir}/*.fasta; do
        if [[ -f "\$fasta" && -s "\$fasta" ]]; then
            sample_name=\$(basename "\$fasta" .fasta)
            abs_path=\$(realpath "\$fasta")
            echo -e "\${sample_name}\t\${abs_path}" >> abritamr_list.tab
        fi
    done
    
    # run abritamr
    abritamr run \\
        --contigs abritamr_list.tab \\
        --prefix abritamr_output \\
        ${extra_opts}
    
    # generate html report from summary_matches.txt
    python3 ${projectDir}/bin/abritamr_report.py summary_matches.txt abritamr_report.html
    """
}
