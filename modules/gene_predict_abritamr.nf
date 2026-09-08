process gene_predict_abritamr {
    label 'abritamr'
    
    publishDir "${params.outdir}/gene_prediction/abritamr", mode: 'copy', pattern: "*.txt"
    publishDir "${params.outdir}/gene_prediction/abritamr", mode: 'copy', pattern: "merged_abritamr_amr.tsv"
    publishDir "${params.outdir}/Report", mode: 'copy', pattern: "abritamr_report.html"

    input:
    path collected_dir

    output:
    path "summary_partials.txt"
    path "summary_virulence.txt"
    path "summary_matches.txt"
    path "abritamr_report.html"
    path "merged_abritamr.tsv"
    path "merged_abritamr_amr.tsv"

    script:
    def species_opt = params.mutation?.toString()?.trim() ? "--species ${params.mutation.toString().trim()}" : ""

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
    mkdir abritamr_results
    cd abritamr_results

    abritamr run \\
        --contigs ../abritamr_list.tab \\
        -j ${params.cpus} \\
        --identity ${params.minid_genes.toDouble() / 100.0} \\
        ${species_opt}

    mv summary_partials.txt ../
    mv summary_virulence.txt ../
    mv summary_matches.txt ../

    cd ..

    # generate html report
    python3 ${projectDir}/python/abritamr_report.py summary_matches.txt abritamr_report.html ${params.minid_genes}

    # merge per-sample amrfinder results
    python3 ${projectDir}/python/merge_abritamr.py \\
        --input_dir abritamr_results \\
        --output merged_abritamr.tsv \\
        --mincov_genes ${params.mincov_genes} \\
        --minid_genes ${params.minid_genes}
    """
}