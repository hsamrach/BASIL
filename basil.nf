#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// All rights reserved. © 2026 BASIL, Samrach Han, and the BASIL development team.
// BASIL is released under GPL lisence 3.0. For terms of use see LICENSE.

def help() {
    log.info"""
    Usage: basil --reads_dir <reads_directory> [options]
           basil --reads_tsv <samples.tsv> [options]
           basil --r1 <read1.fastq(.gz)> --r2 <read2.fastq(.gz)> [options]
           basil --contigs_dir <contigs_directory> [options]

    Input/Output:
    --r1 FILE                           Input read 1 file (default: null)
    --r2 FILE                           Input read 2 file (default: null)
    --reads_dir DIR                     Directory containing paired-end reads (default: null)
    --reads_tsv FILE                    TSV file with columns: sample, read1, read2 (default: null)
    --contigs_dir DIR                   Directory containing pre-assembled contigs (.fasta, .fas, .fa, .fna) (default: null)
    --dir_depth N                       Depth of directory for searching input files in --reads_dir/--contigs_dir (default: 1)
    --outdir DIR                        Output directory (default: BASIL_out)

    Paired-end QC:
    --pe_quality_fail_rate N            Maximum percentage of low-quality bases allowed (default: 40). "-u" in fastp
    --pe_min_length N                   Minimum read length required (default: 50). "--length_required" in fastp
    --pe_base_depth N                   Phred score threshold for qualified bases (default: 15). "-q" in fastp
    --pe_no_correction                  Disable read correction for paired-end reads (default: enabled). "--correction" in fastp
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
    --checkm2_db FILE                   Path to CheckM2 database (mandatory) "/path/checkm2_database/uniref100.KO.1.dmnd/" (saved for future runs)
    --kraken2_db DIR                    Path to Kraken2 database (saved for future runs, auto-runs if available)
    --gtdbtk_db DIR                     Path to GTDB-Tk database (saved for future runs, auto-runs if available)
    --minid_genes N                     Minimum identity for abricate & abritamr (default: 90)
    --mincov_genes N                    Minimum coverage for abricate (default: 90). abritamr is constantly used 90% coverage (unmodifiedable).
    --mutation "STRING"                 Provide species name for point mutation analysis by abritamr (default: null)
                                        Available species: "Acinetobacter_baumannii,Burkholderia_cepacia,Burkholderia_pseudomallei,
                                        Burkholderia_mallei,Campylobacter,Citrobacter_freundii,Clostridioides_difficile,
                                        Corynebacterium_diphtheriae,Enterobacter_asburiae,Enterobacter_cloacae,Enterococcus_faecalis,
                                        Enterococcus_faecium,Escherichia,Klebsiella_oxytoca,Klebsiella_pneumoniae,Neisseria_gonorrhoeae,
                                        Neisseria_meningitidis,Pseudomonas_aeruginosa,Salmonella,Serratia_marcescens,Staphylococcus_aureus,
                                        Staphylococcus_pseudintermedius,Streptococcus_agalactiae,Streptococcus_pneumoniae,Streptococcus_pyogenes,
                                        Vibrio_cholerae,Vibrio_vulfinicus,Vibrio_parahaemolyticus"

    Resources control:
    --meta_merge "f1.csv,f2.xls,..."    Merge metadata and genomic output file (.xlsx, .xls, .tsv, .csv, .tab) (default: null)
                                        Use it independently.
    --parallel_run N                    Number of sample runs in parallel (default: 1)
    --cpus N                            CPUs in GB per sample (default: 8)
    --ram N                             RAM in GB per sample (default: 16)
    -resume                             Resume work (built-in nextflow function)
    -profile "STRING"                   Alternative use of profile platform (choices: apptainer/singularity/docker/mamba, default: apptainer)
    --clear_saved_db "STRING"           Clear saved database paths (choices: checkm2, kraken2, gtdbtk, all)
    --version                           Show version and exit
    --help                              Show this help message and exit
    """
}

def version() {
    log.info """
    BASIL Version 1.1
    """.stripIndent()
}

def validateParams(params) {
    def VALID_PARAMS = [
        'r1', 'r2', 'reads_dir', 'reads_tsv', 'contigs_dir', 'dir_depth',
        'outdir', 'pe_quality_fail_rate', 'pe_min_length', 'pe_base_depth',
        'pe_no_correction', 'pe_extra_opt', 'min_genome_cov', 'max_genome_cov',
        'genome_size', 'min_contig_length', 'min_contig_cov', 'asm_used', 'skip_polish',
        'only_filtered', 'only_polished', 'checkm2_db', 'kraken2_db', 'gtdbtk_db',
        'minid_genes', 'mincov_genes', 'mutation', 'parallel_run', 'cpus', 'ram',
        'help', 'version', 'clear_saved_db', 'constants', 'meta_merge'
    ] as Set

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

def getConfigDir() {
    def configDir = new File(System.getProperty("user.home"), ".basil")
    if (!configDir.exists()) {
        configDir.mkdirs()
    }
    return configDir
}

def loadConfig() {
    def configFile = new File(getConfigDir(), "config.json")
    if (!configFile.exists()) {
        return [:]
    }
    
    try {
        def slurper = new groovy.json.JsonSlurper()
        return slurper.parse(configFile)
    } catch (Exception e) {
        log.warn "Failed to load config: ${e.message}"
        return [:]
    }
}

def saveConfig(config) {
    def configFile = new File(getConfigDir(), "config.json")
    try {
        def json = new groovy.json.JsonOutput().toJson(config)
        configFile.text = groovy.json.JsonOutput.prettyPrint(json)
        log.info "Configuration saved to: ${configFile.absolutePath}"
    } catch (Exception e) {
        log.warn "Failed to save config: ${e.message}"
    }
}

def getOrSaveDatabase(params, configKey, paramKey, dbType = 'dir') {
    def config = loadConfig()
    def savedPath = config[configKey]
    def providedPath = params[paramKey]
    
    // If user provided a path, validate and save it
    if (providedPath) {
        def dbFile = file(providedPath)
        if (dbType == 'dir' && !dbFile.isDirectory()) {
            error "ERROR: --${paramKey} must be a directory (got: '${providedPath}')"
        } else if (dbType == 'file' && !dbFile.isFile()) {
            error "ERROR: --${paramKey} must be a file (got: '${providedPath}')"
        }
        
        config[configKey] = providedPath
        saveConfig(config)
        log.info "Saved ${configKey} to configuration"
        return providedPath
    }
    
    // If no path provided, try to use saved path
    if (savedPath) {
        def dbFile = file(savedPath)
        if (dbType == 'dir' && dbFile.isDirectory()) {
            log.info "Using saved ${configKey}: ${savedPath}"
            return savedPath
        } else if (dbType == 'file' && dbFile.isFile()) {
            log.info "Using saved ${configKey}: ${savedPath}"
            return savedPath
        } else {
            log.warn "Saved ${configKey} no longer exists: ${savedPath}"
            config.remove(configKey)
            saveConfig(config)
            return null
        }
    }
    
    return null
}

def clearSavedDatabase(configKey) {
    def config = loadConfig()
    if (config.remove(configKey)) {
        saveConfig(config)
        log.info "Cleared saved ${configKey}"
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
include {speciation_kraken2_contigs_dir} from './modules/speciation_kraken2.nf'
include {assembly_collection} from './modules/assembly_collection.nf'
include {format_convert_contigs} from './modules/assembly_collection.nf'
include {speciation_gtdbtk} from './modules/speciation_gtdbtk.nf'
include {gene_predict_abritamr} from './modules/gene_predict_abritamr.nf'
include {gene_predict_abricate} from './modules/gene_predict_abricate.nf'
include {gene_predict_amr_consensus} from './modules/amr_consensus.nf'
include {main_report} from './modules/main_report.nf'
include {meta_merge}   from './modules/meta_merge.nf'

// function to create input channel from various sources
def createInputChannel() {
    if (params.reads_tsv ?: false) {
        log.info "Using TSV input: ${file(params.reads_tsv).toAbsolutePath()}"
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
        log.info "Searching for reads in: ${file(params.reads_dir).toAbsolutePath()}"

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
                No paired-end reads found in: ${file(params.reads_dir).toAbsolutePath()}

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
        log.info "Using individual input: ${sample_id}, R1=${params.r1}, R2=${params.r2}"
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
          4. Contigs dir    : --contigs_dir /path/to/contigs/

        Run with --help for more information.
        """
    }
}

def createContigsInputChannel() {
    if (!(params.contigs_dir ?: false)) {
        return channel.empty()
    }

    log.info "Searching for contig files in: ${file(params.contigs_dir).toAbsolutePath()}"
    
    def clean_path = params.contigs_dir.replaceAll('/$', '')
    
    // Pattern for FASTA files with accepted extensions
    def pattern = "${clean_path}/*.{fasta,fas,fa,fna}"
    
    return channel
        .fromPath(pattern, checkIfExists: false)
        .ifEmpty {
            error """
            No FASTA files found in: ${file(params.contigs_dir).toAbsolutePath()}

            Searched pattern: ${pattern}

            Supported file extensions:
            - .fasta
            - .fas
            - .fa
            - .fna

            Please check:
            1. Path is correct       : ${clean_path}
            2. Files exist           : ls ${clean_path}/*
            3. Files have supported extensions
            """
        }
        .map { fasta_file ->
            def sample_id = fasta_file.baseName
            tuple(sample_id, fasta_file)
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

// Handle clear_saved_db parameter
if (params.clear_saved_db) {
    def clearDbOption = params.clear_saved_db?.toString()?.toLowerCase()
    if (clearDbOption == 'checkm2') {
        clearSavedDatabase('checkm2_db')
    } else if (clearDbOption == 'kraken2') {
        clearSavedDatabase('kraken2_db')
    } else if (clearDbOption == 'gtdbtk') {
        clearSavedDatabase('gtdbtk_db')
    } else if (clearDbOption == 'all') {
        clearSavedDatabase('checkm2_db')
        clearSavedDatabase('kraken2_db')
        clearSavedDatabase('gtdbtk_db')
        log.info "All saved database paths have been cleared"
    } else {
        log.warn "Unknown clear_saved_db option: '${clearDbOption}'. Please use one of: checkm2, kraken2, gtdbtk, all"
    }
    System.exit(0)
}

if (params.meta_merge ?: false) {

    def raw_paths = params.meta_merge
        .toString()
        .split(',')
        *.trim()
        .findAll { it }

    if (raw_paths.isEmpty()) {
        error """\
        ERROR: --meta_merge requires a comma-separated list of file paths.
        Usage : --meta_merge 'metadata.xlsx,quast.tsv,checkm2.csv'
        """.stripIndent()
    }

    log.info """
    Files    : ${raw_paths.join(', ')}
    Output   : ${file(params.outdir).toAbsolutePath()}/meta_merged.csv
    """.stripIndent()

    def meta_files_ch = Channel.fromList(
        raw_paths.collect { p -> file(p, checkIfExists: true) }
    )

    meta_merge(meta_files_ch.collect())
    return
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

// load or get databases
def checkm2_db_path = getOrSaveDatabase(params, 'checkm2_db', 'checkm2_db', 'file')
def kraken2_db_path = getOrSaveDatabase(params, 'kraken2_db', 'kraken2_db', 'dir')
def gtdbtk_db_path = getOrSaveDatabase(params, 'gtdbtk_db', 'gtdbtk_db', 'dir')

// CheckM2 is mandatory - validate early
if (!checkm2_db_path) {
    error """
    ERROR: CheckM2 database path not found!
    
    CheckM2 database is mandatory. Please provide it with:
    --checkm2_db /path/to/uniref100.KO.1.dmnd
    
    Once provided, it will be saved and automatically used in future runs.
    """.stripIndent()
}

// Validate CheckM2 database
def checkm2_db = file(checkm2_db_path)
if (!checkm2_db.isFile()) {
    error "ERROR: CheckM2 database file not found at: ${checkm2_db_path}"
}

// checkm2_db channel made
def checkm2_db_ch = channel.value(file(checkm2_db_path))

// Check if using contigs_dir mode
if (params.contigs_dir ?: false) {
    
    // Create input channel from contigs directory
    def contigs_input_ch = createContigsInputChannel()
    
    // Convert formats
    def contig_files_ch = contigs_input_ch.map { sample_id, fasta_file -> fasta_file }
    format_convert_contigs(contig_files_ch.collect())

    
    // Create collected assembly channel
def collected_contig_ch = format_convert_contigs.out.collected_dir
    .map { fasta_file ->
        tuple(fasta_file.baseName, fasta_file)
    }
    .collect()
    .map { list ->
        def sample_ids   = list.collect { it[0] }
        def filtered_fas = list.collect { it[1] }

        tuple(sample_ids, filtered_fas, [file("NO_FILE")])
    }
    
    // Batch analysis: assembly collection
    def needs_collection = gtdbtk_db_path || true

    if (needs_collection) {

        if (kraken2_db_path) {

            kraken2_db_ch = channel.fromPath(
                kraken2_db_path,
                type: 'dir',
                checkIfExists: true
            )

            speciation_kraken2_contigs_dir(
            format_convert_contigs.out.collected_dir.collect(),
            kraken2_db_ch
        )
        }

        if (gtdbtk_db_path) {

            gtdbtk_db_ch = channel.fromPath(
                gtdbtk_db_path,
                type: 'dir',
                checkIfExists: true
            )

            speciation_gtdbtk(
                format_convert_contigs.out.collected_dir.collect(),
                gtdbtk_db_ch
            )
        }

        gene_predict_abricate(
            format_convert_contigs.out.collected_dir.collect()
        )

        gene_predict_abritamr(
            format_convert_contigs.out.collected_dir.collect()
        )

        gene_predict_amr_consensus(
            gene_predict_abricate.out[0],
            gene_predict_abritamr.out[2]
        )
    }
    
    // Aggregate reports from contigs_dir mode
    def all_done = Channel.empty()

    if (kraken2_db_path) {
        all_done = all_done.mix(
            speciation_kraken2_contigs_dir.out.kraken2_contigs_reports
        )
    }

    if (gtdbtk_db_path) {
        all_done = all_done.mix(
            speciation_gtdbtk.out.gtdbtk_reports
        )
    }
    
    if (needs_collection) {
        all_done = all_done
            .mix(gene_predict_abricate.out[0])
            .mix(gene_predict_abricate.out[1])
            .mix(gene_predict_abritamr.out[0])
            .mix(gene_predict_abritamr.out[1])
            .mix(gene_predict_abritamr.out[2])
            .mix(gene_predict_abritamr.out[3])
            .mix(gene_predict_amr_consensus.out.html_report)
            .mix(gene_predict_amr_consensus.out.tsv_summary)
    }
    
    main_report(all_done.collect())
    
    return
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
assembly_QC(final_assembly_ch, checkm2_db_ch)

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
if (kraken2_db_path) {
    kraken2_db_ch = channel.fromPath(kraken2_db_path, type: 'dir', checkIfExists: true)

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
def needs_collection = gtdbtk_db_path ||
                       true  // Abricate and AbritAMR always run

if (needs_collection) {
    assembly_collection(collected_assembly_ch)

    if (gtdbtk_db_path) {
        gtdbtk_db_ch = channel.fromPath(gtdbtk_db_path, type: 'dir', checkIfExists: true)
        speciation_gtdbtk(assembly_collection.out.collected_dir, gtdbtk_db_ch)
    }

    gene_predict_abricate(assembly_collection.out.collected_dir)

    gene_predict_abritamr(assembly_collection.out.collected_dir)

    gene_predict_amr_consensus(
        gene_predict_abricate.out[0],
        gene_predict_abritamr.out[2]
    )
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

if (kraken2_db_path) {
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

if (gtdbtk_db_path) {
    all_done = all_done.mix(speciation_gtdbtk.out.gtdbtk_reports)
}

// Abricate and AbritAMR reports always included
if (needs_collection) {
    all_done = all_done
        .mix(gene_predict_abricate.out[0])
        .mix(gene_predict_abricate.out[1])
        .mix(gene_predict_abritamr.out[0])
        .mix(gene_predict_abritamr.out[1])
        .mix(gene_predict_abritamr.out[2])
        .mix(gene_predict_abritamr.out[3])
        .mix(gene_predict_amr_consensus.out.html_report)
        .mix(gene_predict_amr_consensus.out.tsv_summary)
}

main_report(all_done.collect())
}
