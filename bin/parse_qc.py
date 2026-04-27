#!/usr/bin/env python3

import sys

qc_file          = sys.argv[1]
mean_read_length = 0
coverage         = 0

# parse QC summary
with open(qc_file) as f:
    for line in f:
        if line.startswith("Mean read length:"):
            mean_read_length = float(line.split(":")[1].strip())
        elif line.startswith("Estimated coverage:"):
            coverage = float(line.split(":")[1].strip().replace("X", "").strip())


# k-mer calculation
def calculate_kmer_sizes(avr_len, cov):
    kspades = "21,33"

    if cov > 50:
        if avr_len >= 74:
            kspades += ",55"
        if avr_len >= 96:
            kspades += ",77"
        if avr_len >= 118:
            kspades += ",99"
        if avr_len >= 145:
            k = int(avr_len - 30)
            if k % 2 == 0:
                k += 1
            if k >= 127:
                k = 127
            kspades += f",{k}"
    else:
        if avr_len >= 82:
            kspades += ",55"
        if avr_len >= 115:
            kspades += ",77"
        if avr_len >= 151:
            kspades += ",99"
        if avr_len >= 251:
            kspades += ",127"

    return kspades


kmers_spades = calculate_kmer_sizes(mean_read_length, coverage)
assembly_mode = "--isolate" if coverage > 50 else "--careful"

# write config
with open("assembly_config.txt", "w") as f:
    f.write(f"MEAN_READ_LENGTH={mean_read_length}\n")
    f.write(f"COVERAGE={coverage}\n")
    f.write(f"KMERS_SPADES={kmers_spades}\n")
    f.write(f"ASSEMBLY_MODE={assembly_mode}\n")

print(f"Parsed QC metrics:")
print(f"  Mean read length: {mean_read_length:.2f} bp")
print(f"  Coverage:         {coverage:.2f}X")
print(f"  K-mers (SPAdes):  {kmers_spades}")
print(f"  Assembly mode:    {assembly_mode}")