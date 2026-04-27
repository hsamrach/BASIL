#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// arguments typo validation
@groovy.transform.Field
static final Set VALID_PARAMS = [
    'r1', 'r2', 'reads_dir', 'reads_tsv', 'dir_depth',
    'outdir', 'pe_quality_fail_rate', 'pe_min_length', 'pe_base_depth',
    'pe_no_correction', 'pe_extra_opt', 'min_genome_cov', 'max_genome_cov',
    'genome_size', 'min_contig_length', 'min_contig_cov', 'asm_used', 'skip_polish',
    'only_filtered', 'only_polished', 'checkm2_db', 'kraken2_db', 'gtdbtk_db',
    'abritamr_opt', 'abricate_opt', 'parallel_run', 'cpus', 'ram', 
    'kmer_min_count', 'kmer_size', 'genome_size_max', 'genome_size_min', 
    'genome_size_fallback', 'warn_assembly_cov', 'quast_min_contig',
    'pol_mindepth', 'pol_minmq', 'pol_minqual', 'pol_mingap', 
    'pol_bwa_B', 'pol_bwa_K', 'pol_bwa_L', 'help', 'version'
] as Set

def help() {
    log.info"""
    Usage: basil --reads_dir <reads_directory> [options]
           basil --reads_tsv <samples.tsv> [options]
           basil --r1 <read1.fastq(.gz)> --r2 <read2.fastq(.gz)> [options]

    Input/Output:
    --r1 FILE                           Input read 1 file (default: null)
    --r2 FILE                           Input read 2 file (default: null)
    --reads_dir DIR                     Directory containing paired-end reads (default: null)
    --reads_tsv FILE                    TSV file with columns: sample, read1, read2 (default: null)
    --dir_depth N                       Depth of directory for searching reads in --reads_dir (default: 1)
    --outdir DIR                        Output directory (default: BASIL_out)

    Paired-end QC:
    --pe_quality_fail_rate N            Maximum percentage of low-quality bases allowed (default: 40)
    --pe_min_length N                   Minimum read length required (default: 50)
    --pe_base_depth N                   Phred score threshold for qualified bases (default: 15)
    --pe_no_correction                  Disable read correction for paired-end reads (default: enabled)
    --pe_extra_opt "STRING"             Extra fastp options — note: 5 options are already used by default (start with pe_*)
    --min_genome_cov N                  Minimum genome coverage for pass assembly (default: < 20)
    --max_genome_cov N                  Maximum genome coverage for downsampling (default: > 150)
    --genome_size N                     Expected genome size. e.g. --genome_size 5000000 (Automatically calculate if not provided)

    Assembly analysis:
    --asm_used "STRING"                 SPAdes output to use for downstream analysis (choices: contigs/scaffolds, default: contigs)
    --min_contig_length N               Minimum contig length for filtered assembly (default: ≤ 300)
    --min_contig_cov N                  Minimum contig coverage for filtered assembly (default: ≤ 2)
    --skip_polish                       Skip polishing step (default: disabled)
    --only_filtered                     Only filtered contigs will be submitted to process speciation, and gene prediction (default: disabled)
    --only_polished                     Only polished contigs will be submitted to process speciation, and gene prediction (default: disabled)
    --checkm2_db FILE                   Path to CheckM2 database (mandatory) "/path/checkm2_database/uniref100.KO.1.dmnd/" (default: null)
    --kraken2_db DIR                    Path to Kraken2 database to trigger kraken2 step (default: skipped)
    --gtdbtk_db DIR                     Path to GTDB-Tk database to trigger GTDB-Tk step (default: skipped)
    --abritamr_opt "STRING"             Use at least an option of abritamr to trigger abritamr step (default: "skipped")
    --abricate_opt "STRING"             Use at least an option of abricate to trigger abricate step (default: "skipped")

    Resources control:
    --parallel_run N                    Number of sample runs in parallel (default: 1)
    --cpus N                            CPUs in GB per sample (default: 8)
    --ram N                             RAM in GB per sample (default: 16)
    -resume                             Resume work (built-in nextflow function)
    -profile "STRING"                   Alternative use of profile platform (choices: apptainer/singularity/docker/mamba, default: apptainer)
    --version                           Show version and exit
    --help                              Show this help message and exit
    """
}

def version() {
    log.info """
    BASIL Version 1.0
    """.stripIndent()
}

def validateParams(params) {
    def invalid = params.keySet() - VALID_PARAMS
    if (invalid) {
        def formatted = invalid.sort().collect { "--${it}" }.join(', ')
        error """\
        ERROR: Unrecognised option(s): ${formatted}

        This may be a typo or an unsupported option.
        Run with --help for a full list of supported parameters.
        """.stripIndent()
    }
}

// call modules
include {paired_end_QC} from './modules/paired_end_QC.nf'
include {read_qc_html_report} from './modules/read_qc_html_report.nf'
include {dual_assembly} from './modules/dual_assembly.nf'
include {filter_contig} from './modules/filter_contig.nf'
include {polish_contig} from './modules/polish_contig.nf'
include {assembly_html_report} from './modules/assembly_html_report.nf'
include {assembly_QC} from './modules/assembly_QC.nf'
include {speciation_kraken2_reads} from './modules/speciation_kraken2.nf'
include {speciation_kraken2_polished} from './modules/speciation_kraken2.nf'
include {speciation_kraken2_filtered} from './modules/speciation_kraken2.nf'
include {assembly_collection} from './modules/assembly_collection.nf'
include {speciation_gtdbtk} from './modules/speciation_gtdbtk.nf'
include {gene_predict_abritamr} from './modules/gene_predict_abritamr.nf'
include {gene_predict_abricate} from './modules/gene_predict_abricate.nf'
include {main_report} from './modules/main_report.nf'

// function to create input channel from various sources
def createInputChannel() {
    if (params.reads_tsv ?: false) {
        log.info "Using TSV input: ${params.reads_tsv}"
        return channel
            .fromPath(params.reads_tsv, checkIfExists: true)
            .splitCsv(header: false, sep: '\t')
            .map { row ->
                def sample = row[0]
                def r1 = file(row[1], checkIfExists: true)
                def r2 = file(row[2], checkIfExists: true)
                tuple(sample, r1, r2)
            }

    } else if (params.reads_dir) {
        log.info "Searching for reads in: ${params.reads_dir}"

        def clean_path = params.reads_dir.replaceAll('/$', '')

        def pattern_root = "${clean_path}/*_{R,}{1,2}*.{fastq,fq}*"
        def pattern_sub  = "${clean_path}/**/*_{R,}{1,2}*.{fastq,fq}*"

        def ch_root = Channel.fromFilePairs(pattern_root, flat: true, checkIfExists: false)

        def ch_sub = (params.dir_depth > 1)
            ? Channel.fromFilePairs(pattern_sub, flat: true, checkIfExists: false, maxDepth: params.dir_depth)
            : Channel.empty()

        def ch = ch_root.mix(ch_sub)

        return ch
            .ifEmpty {
                error """
                No paired-end reads found in: ${params.reads_dir}

                Tried patterns:
                - ROOT : ${pattern_root}
                - SUB  : ${pattern_sub ?: 'N/A'}

                Search depth  : ${params.dir_depth == -1 ? 'unlimited' : params.dir_depth}

                Expected naming conventions:
                - sample_R1.fastq(.gz*) / sample_R2.fastq(.gz*)
                - sample_1.fastq(.gz*)  / sample_2.fastq(.gz*)
                - sample_R1.fq(.gz*)    / sample_R2.fq(.gz*)
                - sample_1.fq(.gz*)     / sample_2.fq(.gz*)

                Please check:
                1. Path is correct : ${clean_path}
                2. Files exist     : ls ${clean_path}/*_*{1,2}*
                3. Files are readable
                """
            }
            .map { sample_id, r1, r2 ->
                if ([r1, r2].any { it.name.endsWith('.zip') }) {
                    error "ZIP files are not supported: ${r1}, ${r2}. Please convert to .fastq.gz"
                }
                tuple(sample_id, r1, r2)
            }

    } else if (params.r1 && params.r2) {
        def sample_id = file(params.r1).name
            .replaceAll(/_R?1.*$/, '')
            .replaceAll(/\.(fastq|fq)(\.gz)?$/, '')
        log.info "Using single sample input: ${sample_id}, R1=${params.r1}, R2=${params.r2}"
        return channel
            .of(tuple(sample_id, file(params.r1, checkIfExists: true), file(params.r2, checkIfExists: true)))

    } else if ((params.help ?: false)) {
        return channel.empty()

    } else {
        error """
        ERROR: No input provided!

        Choose one of the following input methods:
          1. TSV file       : --reads_tsv samples.tsv
          2. Directory path : --reads_dir /path/to/reads/
          3. Single sample  : --r1 sample_R1.fastq(.gz) --r2 sample_R2.fastq(.gz)

        Run with --help for more information.
        """
    }
}

workflow {
    main:

validateParams(params)

if (params.help) {
help()
System.exit(0)
}

if (params.version) {
version()
System.exit(0)
}

    log.info """
    ╔══════════════════════════════════════════╗
    ║               BASIL Workflow             ║
    ╚══════════════════════════════════════════╝

    Output_dir : ${file(params.outdir).toAbsolutePath()}
    """.stripIndent()

// option validation
if (params.only_filtered && params.only_polished) {
    error "ERROR: --only_filtered and --only_polished cannot be used together"
}
if (params.only_polished && params.skip_polish) {
    error "ERROR: --only_polished cannot be used with --skip_polish"
}
if (params.checkm2_db == null) {
    error "ERROR: --checkm2_db is mandatory. Provide the path to the CheckM2 database file."
}

def checkm2_db = file(params.checkm2_db)

if (!checkm2_db.isFile()) {
    error "ERROR: --checkm2_db must be a file (got: '${params.checkm2_db}'). " +
          "Provide the full path to the .dmnd file, not the directory."
}
// create input channel
def input_ch = createInputChannel()

// execute paired-end QC
paired_end_QC(input_ch)

// branch pass / fail based on qc_status.txt
paired_end_QC.out.reads
    .branch {
        pass: it[7].text.trim() == "PASS"
        fail: true
    }
    .set { qc_branch }

// log skipped samples
qc_branch.fail.view { sample_id, r1, r2, merged, unmerged_r1, unmerged_r2, qc, status ->
    "SKIPPED ${sample_id} — ${status.text.trim()}"
}

// strip qc_status.txt before passing downstream
def pass_reads_ch = qc_branch.pass
    .map { sample_id, r1, r2, merged, unmerged_r1, unmerged_r2, qc, status ->
        tuple(sample_id, r1, r2, merged, unmerged_r1, unmerged_r2, qc)
    }

// all samples (pass + fail) — derived from branch, not raw channel
def all_samples_ch = qc_branch.pass.mix(qc_branch.fail)

// all reads for kraken2 — derived from all_samples_ch
def all_reads_ch = all_samples_ch
    .map { sample_id, r1, r2, merged, unmerged_r1, unmerged_r2, qc, status ->
        tuple(sample_id, r1, r2, merged, unmerged_r1, unmerged_r2, qc)
    }

// multiMap all_samples_ch for report inputs — avoids consuming it twice
all_samples_ch
    .multiMap { sample_id, r1, r2, merged, unmerged_r1, unmerged_r2, qc, status ->
        qc_files:     qc
        status_files: status
    }
    .set { report_inputs }

// execute dual assembly
dual_assembly(pass_reads_ch)

// execute polishing and filtering
def final_assembly_ch
def filter_input_ch

if (!params.skip_polish) {

    def polish_input_ch = dual_assembly.out.best_assembly
        .join(
            pass_reads_ch.map { sample_id, r1, r2, merged, unmerged_r1, unmerged_r2, qc ->
                tuple(sample_id, r1, r2)
            }
        )

    polish_contig(polish_input_ch)

    filter_input_ch = polish_contig.out.polished_assembly
        .map { fasta ->
            def sample_id = fasta.baseName.replaceAll('_polished$', '')
            tuple(sample_id, fasta)
        }

} else {

    filter_input_ch = dual_assembly.out.best_assembly
}

filter_contig(filter_input_ch)

if (!params.skip_polish) {

    // join post-polish filtered + polished assemblies
    final_assembly_ch = filter_contig.out.filtered_assembly
        .join(
            polish_contig.out.polished_assembly.map { fasta ->
                def sample_id = fasta.baseName.replaceAll('_polished$', '')
                tuple(sample_id, fasta)
            }
        )

} else {

    // skip_polish mode carries only filtered assemblies downstream
    final_assembly_ch = filter_contig.out.filtered_assembly
        .map { sample_id, filtered_fasta ->
            tuple(sample_id, filtered_fasta, file("NO_FILE"))
        }
}

// assembly_QC runs on filtered assemblies and polished assemblies when available
assembly_QC(final_assembly_ch)

// execute qc_summary html report
read_qc_html_report(
    report_inputs.qc_files.collect(),
    report_inputs.status_files.collect()
)

// execute assembly html report
def polishing_report_ch = params.skip_polish
    ? Channel.of(file("NO_FILE"))
    : polish_contig.out.polishing_report.collect()

assembly_html_report(
    dual_assembly.out.summary
        .mix(filter_contig.out.filter_report)
        .collect(),
    polishing_report_ch
)

// execute kraken2
if (params.kraken2_db) {
    kraken2_db_ch = channel.fromPath(params.kraken2_db, type: 'dir', checkIfExists: true)

    // kraken2_reads runs first, on all samples including QC-failed ones
    kraken2_reads_collected_ch = all_reads_ch
        .map { sample_id, r1, r2, merged, unmerged_r1, unmerged_r2, qc -> tuple(sample_id, r1, r2) }
        .collect()
        .map { list ->
            def items  = list.collate(3)
            def ids    = items.collect { it[0] }
            def r1s    = items.collect { it[1] }
            def r2s    = items.collect { it[2] }
            tuple(ids, r1s, r2s)
        }

    speciation_kraken2_reads(kraken2_reads_collected_ch, kraken2_db_ch)

    // kraken2 on polished assembly — skipped if only_filtered or skip_polish is set
    if (!params.only_filtered && !params.skip_polish) {
        kraken2_pol_collected_ch = final_assembly_ch
            .join(
                pass_reads_ch
                    .map { sample_id, r1, r2, merged, unmerged_r1, unmerged_r2, qc -> tuple(sample_id, r1, r2) }
            )
            .collect()
            .map { list ->
                def items               = list.collate(5)
                def sample_ids          = items.collect { it[0] }
                def polished_assemblies = items.collect { it[2] }
                def r1s                 = items.collect { it[3] }
                def r2s                 = items.collect { it[4] }
                tuple(sample_ids, polished_assemblies, r1s, r2s)
            }

        speciation_kraken2_polished(
            kraken2_pol_collected_ch.map { ids, pol, r1s, r2s -> tuple(ids, pol) },
            kraken2_db_ch
        )
    }

    // kraken2 on filtered assembly — skipped only if only_polished is set
    if (!params.only_polished) {
        kraken2_filt_collected_ch = final_assembly_ch
            .join(
                pass_reads_ch
                    .map { sample_id, r1, r2, merged, unmerged_r1, unmerged_r2, qc -> tuple(sample_id, r1, r2) }
            )
            .collect()
            .map { list ->
                def items               = list.collate(5)
                def sample_ids          = items.collect { it[0] }
                def filtered_fastas     = items.collect { it[1] }
                def r1s                 = items.collect { it[3] }
                def r2s                 = items.collect { it[4] }
                tuple(sample_ids, filtered_fastas, r1s, r2s)
            }

        speciation_kraken2_filtered(
            kraken2_filt_collected_ch.map { ids, filt, r1s, r2s -> tuple(ids, filt) },
            kraken2_db_ch
        )
    }
}

// batch analysis: collect all assembled genomes
def collected_assembly_ch = final_assembly_ch
    .collect()
    .map { list ->
        def items               = list.collate(3)
        def sample_ids          = items.collect { it[0] }
        def filtered_fastas     = items.collect { it[1] }
        def polished_assemblies = items.collect { it[2] }

        if (params.only_polished) {
            tuple(sample_ids, [file("NO_FILE")], polished_assemblies)
        } else if (params.only_filtered || params.skip_polish) {
            tuple(sample_ids, filtered_fastas, [file("NO_FILE")])
        } else {
            tuple(sample_ids, filtered_fastas, polished_assemblies)
        }
    }

// only collect assemblies if downstream analysis needs them
def needs_collection = params.gtdbtk_db           ||
                       params.abricate_opt?.trim() ||
                       params.abritamr_opt?.trim()

if (needs_collection) {
    assembly_collection(collected_assembly_ch)

    if (params.gtdbtk_db) {
        gtdbtk_db_ch = channel.fromPath(params.gtdbtk_db, type: 'dir', checkIfExists: true)
        speciation_gtdbtk(assembly_collection.out.collected_dir, gtdbtk_db_ch)
    }

    if (params.abricate_opt?.trim()) {
        gene_predict_abricate(assembly_collection.out.collected_dir)
    }

    if (params.abritamr_opt?.trim()) {
        gene_predict_abritamr(assembly_collection.out.collected_dir)
    }
}

// aggregate reports
def all_done = paired_end_QC.out.json
    .mix(assembly_QC.out.quast_reports)
    .mix(assembly_QC.out.checkm2_reports)

if (!params.skip_polish) {
    all_done = all_done
        .mix(assembly_QC.out.quast_polish_reports)
        .mix(polish_contig.out.polishing_report)
}

if (params.kraken2_db) {
    // reads always included
    all_done = all_done
        .mix(speciation_kraken2_reads.out.kraken2_reads_reports)

    if (!params.only_polished) {
        all_done = all_done
            .mix(speciation_kraken2_filtered.out.kraken2_filt_reports)
    }

    if (!params.only_filtered && !params.skip_polish) {
        all_done = all_done
            .mix(speciation_kraken2_polished.out.kraken2_polished_reports)
    }
}

if (params.gtdbtk_db) {
    all_done = all_done.mix(speciation_gtdbtk.out.gtdbtk_reports)
}

main_report(all_done.collect())
}
