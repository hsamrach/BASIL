process gene_predict_amr_consensus {

    label 'amr_consensus'

    publishDir "${params.outdir}/gene_prediction", mode: 'copy', pattern: "amr_consensus.tsv"
    publishDir "${params.outdir}/gene_prediction", mode: 'copy', pattern: "amr_crosscheck.tsv"
    publishDir "${params.outdir}/Report", mode: 'copy', pattern: "amr_consensus_report.html"
    publishDir "${params.outdir}/gene_prediction", mode: 'copy', pattern: "amr_sequences.fasta"

    input:
    path(abricate_tsv_dir)
    path(abritamr_amr_tsv)
    path(reference_dir)

    output:
    path("amr_consensus.tsv"), emit: consensus_tsv
    path("amr_consensus_report.html"), emit: html_report
    path("amr_crosscheck.tsv")
    path("amr_sequences.fasta")

    script:
    """
    set -euo pipefail

    if [ -d "${abricate_tsv_dir}" ]; then
        ABRICATE_DIR="${abricate_tsv_dir}"
    else
        mkdir -p abricate_tsvs
        cp ${abricate_tsv_dir} abricate_tsvs/ 2>/dev/null || true
        ABRICATE_DIR="abricate_tsvs"
    fi

    if [ -d "${reference_dir}" ]; then
        REF_DIR="${reference_dir}"
    else
        mkdir -p reference_genomes
        for item in ${reference_dir}; do
            if [ -d "\$item" ]; then
                find "\$item" -maxdepth 1 -type f \\( -name '*.fasta' -o -name '*.fa' -o -name '*.fna' -name '*.faa' \\) -exec cp {} reference_genomes/ \\;
            elif [ -f "\$item" ]; then
                cp "\$item" reference_genomes/
            fi
        done
        REF_DIR="reference_genomes"
    fi

    NCBI_TSV=""
    CARD_TSV=""
    RESFINDER_TSV=""
    ARGANNOT_TSV=""

    shopt -s nullglob
    for tsv in "\${ABRICATE_DIR}"/*.tsv; do
        base=\$(basename "\$tsv")
        case "\${base,,}" in
            ncbi.tsv) NCBI_TSV="\$tsv" ;;
            card.tsv) CARD_TSV="\$tsv" ;;
            resfinder.tsv) RESFINDER_TSV="\$tsv" ;;
            argannot.tsv) ARGANNOT_TSV="\$tsv" ;;
        esac
    done

    CROSSCHECK_ARGS=(
        --reference-dir "\${REF_DIR}"
        --outdir .
        --abritamr "${abritamr_amr_tsv}"
        --threads "${params.cpus}"
        --coverage "${params.constants.cross_cov}"
        --identity "${params.constants.cross_id}"
    )
    [ -n "\$NCBI_TSV" ] && CROSSCHECK_ARGS+=( --ncbi "\$NCBI_TSV" )
    [ -n "\$CARD_TSV" ] && CROSSCHECK_ARGS+=( --card "\$CARD_TSV" )
    [ -n "\$RESFINDER_TSV" ] && CROSSCHECK_ARGS+=( --resfinder "\$RESFINDER_TSV" )
    [ -n "\$ARGANNOT_TSV" ] && CROSSCHECK_ARGS+=( --argannot "\$ARGANNOT_TSV" )

    python3 ${projectDir}/python/crosscheck_amr.py "\${CROSSCHECK_ARGS[@]}"

    python3 ${projectDir}/python/amr_consensus_report.py \\
        --consensus_tsv amr_consensus.tsv \\
        --output amr_consensus_report.html \\
        --minid_genes "${params.minid_genes}" \\
        --mincov_genes "${params.mincov_genes}"
    """
}
