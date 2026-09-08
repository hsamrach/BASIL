process meta_merge {
    label 'meta_merge'
    tag  "meta_merge"

    publishDir "${params.outdir}/", mode: 'copy', pattern: '*.csv'

    input:
    path input_files

    output:
    path "meta_merged.csv", emit: merged_table

    script:

    def file_args = (input_files instanceof List
        ? input_files.collect { it.name }
        : [input_files.name]
    ).join(' ')

    """
    python3 ${projectDir}/python/meta_merge.py ${file_args} -o meta_merged.csv
    """
}
