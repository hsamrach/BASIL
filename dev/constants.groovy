def kmer_size = 21
def kmer_min_count = 10
def genome_size_max = 10000000
def genome_size_min = 600000
def genome_size_fallback = 4500000
def warn_assembly_cov = 40.0
def quast_min_contig = 0
def pol_mindepth = 0.25
def pol_minmq = 60
def pol_minqual = 3
def pol_mingap = 10
def pol_bwa_K = 100000000
def pol_bwa_B = 4
def pol_bwa_L = 20

// Return as a map for access throughout the pipeline
params.constants = [
    kmer_size: kmer_size,
    kmer_min_count: kmer_min_count,
    genome_size_max: genome_size_max,
    genome_size_min: genome_size_min,
    genome_size_fallback: genome_size_fallback,
    warn_assembly_cov: warn_assembly_cov,
    quast_min_contig: quast_min_contig,
    pol_mindepth: pol_mindepth,
    pol_minmq: pol_minmq,
    pol_minqual: pol_minqual,
    pol_mingap: pol_mingap,
    pol_bwa_K: pol_bwa_K,
    pol_bwa_B: pol_bwa_B,
    pol_bwa_L: pol_bwa_L
]
