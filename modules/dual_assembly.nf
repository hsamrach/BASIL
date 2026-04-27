process dual_assembly {
    label 'assembly'

    tag "${sample_id}"
    
    input:
    tuple val(sample_id), path(r1), path(r2), path(merged), path(unmerged_r1), path(unmerged_r2), path(qc_summary)

    output: 
    path "${sample_id}_assembly_summary.txt", emit: summary 
    tuple val(sample_id), path("${sample_id}_best.fasta"), emit: best_assembly

    script:
    """
    ##########################################################
    # Step 1: parse QC summary to get read length and coverage
    ##########################################################
    python3 ${projectDir}/bin/parse_qc.py ${qc_summary}

    #####################################
    # Step 2: load assembly configuration
    #####################################
    source assembly_config.txt
    
    echo "Assembly configuration:"
    echo "  Mean read length: \$MEAN_READ_LENGTH bp"
    echo "  Coverage: \$COVERAGE X"
    echo "  K-mers (paired-only): \$KMERS_SPADES"
    echo "  Assembly mode: \$ASSEMBLY_MODE"

    ############################################    
    # Step 3: assembly 1 - Paired-end reads only
    ############################################
    
    spades.py \\
        --pe1-1 ${r1} \\
        --pe1-2 ${r2} \\
        -k \$KMERS_SPADES \\
        -o ${sample_id}_paired_only \\
        --threads ${params.cpus} \\
        --cov-cutoff auto \\
        \$ASSEMBLY_MODE
    
    ################################################
    # Step 4: assembly 2 - Paired-end + merged reads
    ################################################
    
    spades.py \\
        --pe1-1 ${unmerged_r1} \\
        --pe1-2 ${unmerged_r2} \\
        --pe1-s ${merged} \\
        -k \$KMERS_SPADES \\
        -o ${sample_id}_paired_plus_merged \\
        --threads ${params.cpus} \\
        --cov-cutoff auto \\
        \$ASSEMBLY_MODE
    
    ##############################################
    # Step 5: run briefly QUAST on both assemblies
    ##############################################
    
    # assembly 1
    quast ${sample_id}_paired_only/${params.asm_used}.fasta \\
        -o ${sample_id}_paired_only_quast \\
        --threads ${params.cpus} \\
        --min-contig ${params.quast_min_contig} \\
        --fast
    
    # assembly 2
    quast ${sample_id}_paired_plus_merged/${params.asm_used}.fasta \\
        -o ${sample_id}_paired_plus_merged_quast \\
        --threads ${params.cpus} \\
        --min-contig ${params.quast_min_contig} \\
        --fast

    ############################################    
    # Step 6: compare and select best assembly
    ############################################

    python3 ${projectDir}/bin/compare_assemblies.py ${sample_id} ${params.min_contig_length} ${params.min_contig_cov} ${params.asm_used}
    
    #####################################################
    # Step 7: export the selected assembly for polishing
    #####################################################

    BEST_CONTIGS=\$(cat best_assembly.txt)
    cp "\$BEST_CONTIGS" "${sample_id}_best.fasta"
    echo "Selected assembly exported to ${sample_id}_best.fasta"
"""
}
