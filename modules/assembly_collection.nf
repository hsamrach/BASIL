process assembly_collection {

    input:
    tuple val(sample_ids), path(filtered_fastas), path(polished_assemblies)

    output:
    path "collected_assemblies", emit: collected_dir

    script:
    def sample_list   = sample_ids instanceof List ? sample_ids : [sample_ids]
    def filtered_list = filtered_fastas instanceof List ? filtered_fastas : [filtered_fastas]
    def polished_list = polished_assemblies instanceof List ? polished_assemblies : [polished_assemblies]

    """
    mkdir -p collected_assemblies

    SAMPLE_IDS=(${sample_list.join(' ')})
    FILTERED_FASTAS=(${filtered_list.join(' ')})
    POLISHED_FASTAS=(${polished_list.join(' ')})

    N=\${#SAMPLE_IDS[*]}
    FILT_LEN=\${#FILTERED_FASTAS[*]}
    POL_LEN=\${#POLISHED_FASTAS[*]}

    for i in \$(seq 0 \$(( N - 1 ))); do
        SAMPLE_ID="\${SAMPLE_IDS[\$i]}"

        FILT_IDX=\$i
        if [[ \$FILT_LEN -eq 1 ]]; then FILT_IDX=0; fi
        if [[ -f "\${FILTERED_FASTAS[\$FILT_IDX]}" && -s "\${FILTERED_FASTAS[\$FILT_IDX]}" ]]; then
            cp "\${FILTERED_FASTAS[\$FILT_IDX]}" "collected_assemblies/\$(basename "\${FILTERED_FASTAS[\$FILT_IDX]}")"
        fi

        POL_IDX=\$i
        if [[ \$POL_LEN -eq 1 ]]; then POL_IDX=0; fi
        if [[ -f "\${POLISHED_FASTAS[\$POL_IDX]}" && -s "\${POLISHED_FASTAS[\$POL_IDX]}" ]]; then
            cp "\${POLISHED_FASTAS[\$POL_IDX]}" "collected_assemblies/\$(basename "\${POLISHED_FASTAS[\$POL_IDX]}")"
        fi
    done

    echo "total assemblies collected: \$(ls collected_assemblies/*.fasta | wc -l)"
    """
}
