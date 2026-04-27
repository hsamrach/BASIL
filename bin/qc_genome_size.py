#!/usr/bin/env python3

import sys
import json
import glob

sample            = sys.argv[1]
genome_size_param = float(sys.argv[2]) if sys.argv[2] != "None" else None
max_genome_cov    = float(sys.argv[3])
min_genome_cov    = float(sys.argv[4])  # coverage below this = skip assembly

genome_size_max      = float(sys.argv[5])  # genome above this = overestimated
genome_size_min      = float(sys.argv[6])  # genome below this = underestimated
genome_size_fallback = float(sys.argv[7])  # fallback for both over/under estimation
warn_assembly_cov    = float(sys.argv[8])  # coverage below this = warn but continue

# read fastp json file
with open(f"{sample}.json") as f:
    data = json.load(f)

after            = data["summary"]["after_filtering"]
total_reads      = after["total_reads"]
total_bases      = after["total_bases"]
mean_read_length = total_bases / total_reads if total_reads > 0 else 0

# conditional genome size estimation and correction
original_genome_size = None
correction_reason    = None

if genome_size_param is not None:
    genome_size       = genome_size_param
    estimation_method = "User provided"

else:
    # read unique k-mer count from KMC output
    kmersize_files = glob.glob("*.kmersize")
    if not kmersize_files:
        raise FileNotFoundError("No .kmersize file found — KMC output not parsed correctly")

    with open(kmersize_files[0]) as f:
        content = f.read().strip()

    if not content:
        raise ValueError("Empty .kmersize file — KMC may have failed")

    # genome size estimation based on unique k-mer count
    genome_size       = int(float(content))
    estimation_method = f"Genome size estimated automatically"

    # size corrections — size_cap checked FIRST
    if genome_size > genome_size_max:
        original_genome_size = genome_size
        genome_size           = genome_size_fallback
        correction_reason     = "size_cap"
        estimation_method     = f"Fallback {genome_size_fallback:,} bp (Genome size exceeded {genome_size_max:,} bp cap)"

    elif genome_size < genome_size_min:
        original_genome_size = genome_size
        genome_size           = genome_size_fallback
        correction_reason     = "genome_too_small"
        estimation_method     = f"Fallback {genome_size_fallback:,} bp (Genome size underestimated: {original_genome_size:,.0f} bp)"

# calculate coverage and downsampling factor
coverage           = total_bases / genome_size if genome_size > 0 else 0
needs_downsampling = coverage > max_genome_cov
downsample_factor  = max_genome_cov / coverage if needs_downsampling else 1.0

# additional information to be reported 
warnings = []

if genome_size_param is None:

    if correction_reason == "size_cap":
        warnings.append(f"Genome size exceeded cap ({original_genome_size:,.0f} bp > {genome_size_max:,} bp)")
        warnings.append(f"Automatically set to fallback: {genome_size_fallback:,} bp")
        warnings.append(f"If it is not an expected result, consider using --genome_size for accurate coverage")

    elif correction_reason == "genome_too_small":
        warnings.append(f"Genome size underestimated ({original_genome_size:,.0f} bp < {genome_size_min:,} bp)")
        warnings.append(f"Automatically set to fallback: {genome_size_fallback:,} bp")
        warnings.append(f"If it is not an expected result, consider using --genome_size for accurate coverage")

    elif genome_size < 1e6:
        warnings.append("Estimated genome size is smaller than 1 Mbp")
        warnings.append("If it is not an expected result, consider using --genome_size parameter for accurate coverage")

    # min coverage borderline check
    if coverage < min_genome_cov:
        warnings.append(f"Insufficient sequencing depth for reliable assembly")
        warnings.append(f"This sample will be excluded from assembly")
        warnings.append(f"If you are intended to proceed with assembly, consider using --min_genome_cov to set a lower threshold")
    elif coverage < warn_assembly_cov:
        warnings.append(f"Low coverage warning: {coverage:.1f}X (recommended: ≥{warn_assembly_cov}X)")
        warnings.append(f"Assembly may be fragmented")

# summaries
print(f"  Sample         : {sample}")
print(f"  Total reads    : {total_reads:,}")
print(f"  Total bases    : {total_bases:,}")
print(f"  Read length    : {mean_read_length:.2f} bp")
print(f"  Genome size    : {genome_size:,.0f} bp  [{estimation_method}]")
print(f"  Coverage       : {coverage:.2f}X")
if warnings:
    for w in warnings:
        print(f"  {w}")

# write qc summary
with open(f"{sample}_QC_summary.txt", "w") as out:
    out.write(f"Sample ID: {sample}\n\n")
    out.write("==== FASTP AFTER FILTERING ====\n")
    out.write(f"Total reads: {total_reads}\n")
    out.write(f"Total bases: {total_bases}\n")
    out.write(f"Mean read length: {mean_read_length:.2f}\n\n")

    out.write("==== GENOME SIZE ====\n")
    out.write(f"Method: {estimation_method}\n")
    out.write(f"Genome size: {genome_size:,.0f} bp\n")
    out.write(f"Estimated coverage: {coverage:.2f} X\n")

    if warnings:
        out.write("\n")
        for w in warnings:
            out.write(f"\u26A0 {w}\n")

    if needs_downsampling:
        out.write("\n")
        out.write("==== DOWNSAMPLING ====\n")
        out.write(f"Read coverage exceeds maximum threshold ({max_genome_cov}X)\n")
        out.write(f"Downsampling factor: {downsample_factor:.4f}\n")
        out.write(f"Target coverage: {max_genome_cov}X\n")

# write downsampled flag
with open("downsample.flag", "w") as f:
    f.write(f"{downsample_factor}\n" if needs_downsampling else "SKIP\n")

# write qc status
qc_pass = coverage >= min_genome_cov

with open(f"{sample}_qc_status.txt", "w") as f:
    if qc_pass:
        f.write(f"PASS\n")
    else:        
        f.write(f"FAIL\tInsufficient coverage: {coverage:.1f}X (minimum: {min_genome_cov}X)\n")