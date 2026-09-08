process filter_contig {
    label 'filter'

    tag "${sample_id}"

    publishDir "${params.outdir}/assembled_genomes", mode: 'copy', saveAs: { filename ->
        filename.endsWith("_filt.fasta") ? filename : null
    }

    input:
    tuple val(sample_id), path(assembly_fasta)

    output:
    tuple val(sample_id), path("${sample_id}_filt.fasta"), emit: filtered_assembly
    path "${sample_id}_filtering_report.txt", emit: filter_report

    script:
    """
    INPUT_FILE="${assembly_fasta}"
    MIN_LENGTH=${params.min_contig_length}
    MIN_COV=${params.min_contig_cov}
    OUTPUT_FILE="${sample_id}_filt.fasta"
    REPORT_FILE="${sample_id}_filtering_report.txt"

    echo "Input: \$INPUT_FILE"
    echo "Filtering criteria:"
    echo "  - Minimum length: \$MIN_LENGTH bp"
    echo "  - Minimum coverage: \$MIN_COV X"

    TOTAL_CONTIGS=\$(grep -c "^>" "\$INPUT_FILE")
    TOTAL_BP=\$(seqkit stats -T "\$INPUT_FILE" | awk 'NR==2 {print \$5}')

    echo "Total input contigs: \$TOTAL_CONTIGS"
    echo "Total input bases:   \$TOTAL_BP bp"

    seqkit fx2tab "\$INPUT_FILE" | \\
    awk -F "\\t" -v min_cov="\$MIN_COV" '{
        header = \$1
        seq = \$2

        cov_val = ""
        depth_val = ""

        if (header ~ /_cov_[0-9]+\\.?[0-9]*/) {
            cov_val = header
            sub(/.*_cov_/, "", cov_val)
            sub(/_.*/, "", cov_val)
            sub(/[^0-9.].*/, "", cov_val)
        }

        if (header ~ /_polish_mean_depth_[0-9]+\\.?[0-9]*/) {
            depth_val = header
            sub(/.*_polish_mean_depth_/, "", depth_val)
            sub(/X.*/, "", depth_val)
            sub(/[^0-9.].*/, "", depth_val)
        }

        pass = 0
        if (cov_val != "" && depth_val != "") {
            if (cov_val+0 >= min_cov+0 && depth_val+0 >= min_cov+0) pass = 1
        } else if (cov_val != "") {
            if (cov_val+0 >= min_cov+0) pass = 1
        } else if (depth_val != "") {
            if (depth_val+0 >= min_cov+0) pass = 1
        }

        if (pass) print ">" header "\\n" seq
    }' | \\
    seqkit seq -m "\$MIN_LENGTH" > "\$OUTPUT_FILE"

    if [[ -s "\$OUTPUT_FILE" ]]; then
        PASSED_CONTIGS=\$(grep -c "^>" "\$OUTPUT_FILE")
        PASSED_BP=\$(seqkit stats -T "\$OUTPUT_FILE" | awk 'NR==2 {print \$5}')
    else
        PASSED_CONTIGS=0
        PASSED_BP=0
    fi

    echo ""
    echo "=========================================="
    echo "FILTERING SUMMARY"
    echo "=========================================="
    echo "Total input contigs: \$TOTAL_CONTIGS"
    echo "Total output contigs: \$PASSED_CONTIGS"
    echo ""
    echo "Total input bases:  \$TOTAL_BP bp"
    echo "Total output bases: \$PASSED_BP bp"

    if [[ "\$TOTAL_BP" -gt 0 ]]; then
        RETENTION=\$(awk -v total_in="\$TOTAL_BP" -v total_out="\$PASSED_BP" 'BEGIN {printf "%.2f", (total_out/total_in)*100}')
        echo "Retention rate: \$RETENTION %"
    else
        RETENTION=0
    fi

    echo "Output written to: \$OUTPUT_FILE"

    {
    echo "Sample ID: ${sample_id}"
    echo ""
    echo "======================================================================"
    echo "CONTIG FILTERING REPORT"
    echo "======================================================================"
    echo ""
    echo "Input assembly: \$INPUT_FILE"
    echo ""
    echo "Filtering Criteria:"
    echo "  Minimum length: \$MIN_LENGTH bp"
    echo "  Minimum coverage: \$MIN_COV X"
    echo ""
    echo "Results:"
    echo "  Total input contigs: \$TOTAL_CONTIGS"
    echo "  Total output contigs: \$PASSED_CONTIGS"
    echo ""
    echo "Assembly Size:"
    echo "  Input:  \$TOTAL_BP bp"
    echo "  Output: \$PASSED_BP bp"
    echo "  Retention: \$RETENTION %"
    } > "\$REPORT_FILE"
    """
}
