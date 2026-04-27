process paired_end_QC {
    label 'pe_qc'

    tag "${sample_id}"
    publishDir "${params.outdir}/filtered_reads", mode: 'copy', saveAs: { filename ->
    if (filename.endsWith("_trimmed_R1.fastq.gz")) filename
    else if (filename.endsWith("_trimmed_R2.fastq.gz")) filename
    else if (filename.endsWith("_merged.fastq.gz")) filename
    else if (filename.endsWith("_unmerged_R1.fastq.gz")) filename
    else if (filename.endsWith("_unmerged_R2.fastq.gz")) filename
    else null
}

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id),
          path("${sample_id}_trimmed_R1.fastq.gz"),
          path("${sample_id}_trimmed_R2.fastq.gz"),
          path("${sample_id}_merged.fastq.gz"),
          path("${sample_id}_unmerged_R1.fastq.gz"),
          path("${sample_id}_unmerged_R2.fastq.gz"),
          path("${sample_id}_QC_summary.txt"),
          path("${sample_id}_qc_status.txt"), emit: reads
    path "${sample_id}.json", emit: json

    script:
    def genome_size_value = params.genome_size ?: "None"
    def read_correction_flag = params.pe_no_correction ? "--correction" : ""
    def kmc_cmd = params.genome_size ? "" : """
        # run KMC on R1 only
        mkdir -p kmc_tmp

        echo "${sample_id}_trimmed_R1.fastq.gz" > kmc_input.txt

        kmc \\
            -k${params.kmer_size} \\
            -t${params.cpus} \\
            -ci${params.kmer_min_count} \\
            @kmc_input.txt \\
            ${sample_id}_kmc \\
            kmc_tmp 2>&1 | tee kmc_output.txt

        # parse unique counted k-mers directly from KMC stdout
        grep -i "unique counted k" kmc_output.txt | \\
            awk -F: '{gsub(/ /, "", \$2); print \$2}' > ${sample_id}.kmersize

        rm -rf kmc_tmp kmc_input.txt kmc_output.txt
    """

    """
    ##########
    # 1. FASTP
    ##########

    fastp \\
        -i ${r1} \\
        -I ${r2} \\
        -o ${sample_id}_trimmed_R1.fastq.gz \\
        -O ${sample_id}_trimmed_R2.fastq.gz \\
        -u ${params.pe_quality_fail_rate} \\
        --length_required ${params.pe_min_length} \\
        --thread ${params.cpus} \\
        --json ${sample_id}.json \\
        ${read_correction_flag} \\
        -q ${params.pe_base_depth} \\
        ${params.pe_extra_opt}

    #########################################
    # 2. run KMC if genome size not provided
    #########################################
    ${kmc_cmd}

    ######################################################################
    # 3. conditional coverages and genome size estimation for downsampling
    ######################################################################
    python3 ${projectDir}/bin/qc_genome_size.py ${sample_id} ${genome_size_value} ${params.max_genome_cov} ${params.min_genome_cov} ${params.genome_size_max} ${params.genome_size_min} ${params.genome_size_fallback} ${params.warn_assembly_cov}

    ################################################################
    # 4. if downsampling is needed, perform downsampling with seqkit
    ################################################################
    DOWNSAMPLE_FACTOR=\$(cat downsample.flag)

    if [ "\$DOWNSAMPLE_FACTOR" != "SKIP" ]; then
        echo "Downsampling reads with factor: \$DOWNSAMPLE_FACTOR"

        seqkit sample2 -p \$DOWNSAMPLE_FACTOR -s 11 \\
            -o ${sample_id}_downsampled_R1.fastq.gz ${sample_id}_trimmed_R1.fastq.gz
        seqkit sample2 -p \$DOWNSAMPLE_FACTOR -s 11 \\
            -o ${sample_id}_downsampled_R2.fastq.gz ${sample_id}_trimmed_R2.fastq.gz

        mv ${sample_id}_downsampled_R1.fastq.gz ${sample_id}_trimmed_R1.fastq.gz
        mv ${sample_id}_downsampled_R2.fastq.gz ${sample_id}_trimmed_R2.fastq.gz
    fi

    #######################################
    # 5. merge trimmed reads for assembly 2
    #######################################

    fastp \\
        -i ${sample_id}_trimmed_R1.fastq.gz \\
        -I ${sample_id}_trimmed_R2.fastq.gz \\
        --out1 ${sample_id}_unmerged_R1.fastq.gz \\
        --out2 ${sample_id}_unmerged_R2.fastq.gz \\
        --merge \\
        --merged_out ${sample_id}_merged.fastq.gz \\
        -Q \\
        --disable_length_filtering \\
        --disable_adapter_trimming \\
        --dont_eval_duplication \\
        --thread ${params.cpus} \\
        2>&1 | grep -v "^Detecting" || true

    """
}
