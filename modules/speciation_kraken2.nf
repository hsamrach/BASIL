process speciation_kraken2_filtered {
    label 'kraken2_filtered'
    publishDir "${params.outdir}/speciation/kraken2", mode: 'copy'

    input:
    tuple val(sample_ids), path(filtered_fastas)
    path kraken2_db

    output:
    path "*.report", emit: kraken2_filt_reports
    path "*.output", emit: kraken2_filt_outputs

    script:
    def sample_list   = sample_ids instanceof List ? sample_ids : [sample_ids]
    def filtered_list = filtered_fastas instanceof List ? filtered_fastas : [filtered_fastas]
    """
    SAMPLE_IDS=(${sample_list.join(' ')})
    FILTERED_FASTAS=(${filtered_list.join(' ')})

    for i in "\${!SAMPLE_IDS[@]}"; do
        SAMPLE_ID="\${SAMPLE_IDS[\$i]}"

        kraken2 \\
            --db ${kraken2_db} \\
            --threads ${params.cpus} \\
            --use-names \\
            --output "\${SAMPLE_ID}_filt.output" \\
            --report "\${SAMPLE_ID}_filt.report" \\
            "\${FILTERED_FASTAS[\$i]}"
    done
    """
}

process speciation_kraken2_polished {
    label 'kraken2_polished'
    publishDir "${params.outdir}/speciation/kraken2", mode: 'copy'

    input:
    tuple val(sample_ids), path(polished_assemblies)
    path kraken2_db

    output:
    path "*.report", emit: kraken2_polished_reports
    path "*.output", emit: kraken2_polished_outputs

    script:
    def sample_list   = sample_ids instanceof List ? sample_ids : [sample_ids]
    def polished_list = polished_assemblies instanceof List ? polished_assemblies : [polished_assemblies]
    """
    SAMPLE_IDS=(${sample_list.join(' ')})
    POLISHED_FASTAS=(${polished_list.join(' ')})

    for i in "\${!SAMPLE_IDS[@]}"; do
        SAMPLE_ID="\${SAMPLE_IDS[\$i]}"

        if [[ -f "\${POLISHED_FASTAS[\$i]}" && -s "\${POLISHED_FASTAS[\$i]}" ]]; then
            kraken2 \\
                --db ${kraken2_db} \\
                --threads ${params.cpus} \\
                --use-names \\
                --output "\${SAMPLE_ID}_polished.output" \\
                --report "\${SAMPLE_ID}_polished.report" \\
                "\${POLISHED_FASTAS[\$i]}"
        fi
    done
    """
}

process speciation_kraken2_reads {
    label 'kraken2_reads'
    publishDir "${params.outdir}/speciation/kraken2", mode: 'copy'

    input:
    tuple val(sample_ids), path(r1s), path(r2s)
    path kraken2_db

    output:
    path "*.report", emit: kraken2_reads_reports
    path "*.output", emit: kraken2_reads_outputs

    script:
    def sample_list = sample_ids instanceof List ? sample_ids : [sample_ids]
    def r1_list     = r1s instanceof List ? r1s : [r1s]
    def r2_list     = r2s instanceof List ? r2s : [r2s]
    """
    SAMPLE_IDS=(${sample_list.join(' ')})
    R1S=(${r1_list.join(' ')})
    R2S=(${r2_list.join(' ')})

    for i in "\${!SAMPLE_IDS[@]}"; do
        SAMPLE_ID="\${SAMPLE_IDS[\$i]}"

        kraken2 \\
            --db ${kraken2_db} \\
            --threads ${params.cpus} \\
            --use-names \\
            --paired \\
            --output "\${SAMPLE_ID}_reads.output" \\
            --report "\${SAMPLE_ID}_reads.report" \\
            "\${R1S[\$i]}" "\${R2S[\$i]}"
    done
    """
}
