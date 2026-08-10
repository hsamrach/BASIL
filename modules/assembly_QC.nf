process assembly_QC {
    label 'assembly_QC'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(filtered_fasta), path(polished_assembly)
    path checkm2_db

    output:
    path "${sample_id}_quast_filt", emit: quast_reports
    path "${sample_id}_quast_polished", emit: quast_polish_reports
    path "${sample_id}_checkm2", emit: checkm2_reports

    script:
    """
    # QUAST run
    if [[ -f "${polished_assembly}" && -s "${polished_assembly}" ]]; then
        quast ${polished_assembly} \\
            -o ${sample_id}_quast_polished \\
            --threads ${params.cpus} \\
            --min-contig ${params.constants.quast_min_contig}
    else
        mkdir -p ${sample_id}_quast_polished
    fi
    
    quast ${filtered_fasta} \\
        -o ${sample_id}_quast_filt \\
        --threads ${params.cpus} \\
        --min-contig ${params.constants.quast_min_contig}

    # CheckM2 run
    checkm2 predict \\
        --threads ${params.cpus} \\
        --input . \\
        --output-directory ${sample_id}_checkm2 \\
        --database_path ${checkm2_db} \\
        --extension fasta \\
        --force
    """
}
