#!/usr/bin/env python3

import os
import sys

sample     = sys.argv[1]
min_length = int(sys.argv[2])
min_coverage = float(sys.argv[3])
asm_used   = sys.argv[4]

def parse_quast_report(report_file):
    """parse QUAST report.tsv file"""
    metrics = {}
    if not os.path.exists(report_file):
        return metrics

    with open(report_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("# contigs"):
                metrics['contigs'] = int(line.split("\t")[1])
            elif line.startswith("Total length"):
                metrics['total_length'] = int(line.split("\t")[1])
            elif line.startswith("N50"):
                metrics['n50'] = int(line.split("\t")[1])

    return metrics


def score_assembly(metrics):
    """scoring: prioritize N50, then contig count, then length"""
    n50_score      = metrics.get('n50', 0)
    contig_penalty = metrics.get('contigs', 999999)
    length_score   = metrics.get('total_length', 0)
    return (n50_score * 1000) - (contig_penalty * 100) + (length_score * 0.01)


# parse QUAST reports
assembly1_metrics = parse_quast_report(f"{sample}_paired_only_quast/report.tsv")
assembly2_metrics = parse_quast_report(f"{sample}_paired_plus_merged_quast/report.tsv")

print("\n" + "=" * 70)
print("QUAST METRICS COMPARISON")
print("=" * 70)

print(f"\nAssembly 1 (Paired-only):")
print(f"  Contigs:      {assembly1_metrics.get('contigs', 0):,}")
print(f"  Total length: {assembly1_metrics.get('total_length', 0):,} bp")
print(f"  N50:          {assembly1_metrics.get('n50', 0):,} bp")

print(f"\nAssembly 2 (Paired + Merged):")
print(f"  Contigs:      {assembly2_metrics.get('contigs', 0):,}")
print(f"  Total length: {assembly2_metrics.get('total_length', 0):,} bp")
print(f"  N50:          {assembly2_metrics.get('n50', 0):,} bp")

# best assembly selection
score1 = score_assembly(assembly1_metrics)
score2 = score_assembly(assembly2_metrics)

if score2 > score1:
    best_assembly = 2
    best_contigs  = f"{sample}_paired_plus_merged/{asm_used}.fasta"
else:
    best_assembly = 1
    best_contigs  = f"{sample}_paired_only/{asm_used}.fasta"

print(f"\n{'=' * 70}")
print(f"BEST ASSEMBLY SELECTED: Assembly {best_assembly}")
print(f"{'=' * 70}")

# decide which assembly is better based on metrics
reasons = []
m1 = assembly1_metrics
m2 = assembly2_metrics

if best_assembly == 2:
    if m2.get('n50', 0) > m1.get('n50', 0):
        reasons.append(f"Higher N50 ({m2.get('n50', 0):,} vs {m1.get('n50', 0):,} bp)")
    if m2.get('contigs', 999) < m1.get('contigs', 999):
        reasons.append(f"Fewer contigs ({m2.get('contigs', 0):,} vs {m1.get('contigs', 0):,})")
    if m2.get('total_length', 0) > m1.get('total_length', 0):
        reasons.append(f"Larger genome ({m2.get('total_length', 0):,} vs {m1.get('total_length', 0):,} bp)")
else:
    if m1.get('n50', 0) > m2.get('n50', 0):
        reasons.append(f"Higher N50 ({m1.get('n50', 0):,} vs {m2.get('n50', 0):,} bp)")
    if m1.get('contigs', 999) < m2.get('contigs', 999):
        reasons.append(f"Fewer contigs ({m1.get('contigs', 0):,} vs {m2.get('contigs', 0):,})")
    if m1.get('total_length', 0) > m2.get('total_length', 0):
        reasons.append(f"Larger genome ({m1.get('total_length', 0):,} vs {m2.get('total_length', 0):,} bp)")

if reasons:
    print("\nReasons:")
    for reason in reasons:
        print(f"  {reason}")
else:
    print("  (Similar quality metrics)")

# write assembly summary
with open(f"{sample}_assembly_summary.txt", "w") as out:
    out.write(f"Sample ID: {sample}\n\n")

    # assembly configuration
    with open("assembly_config.txt") as cfg:
        out.write("=" * 70 + "\n")
        out.write("ASSEMBLY CONFIGURATION\n")
        out.write("=" * 70 + "\n")
        for line in cfg:
            out.write(line)
        out.write("\n")

    # assembly 1 metrics
    out.write("=" * 70 + "\n")
    out.write("ASSEMBLY 1: Paired-end reads only\n")
    out.write("=" * 70 + "\n")
    out.write(f"Number of contigs: {m1.get('contigs', 0):,}\n")
    out.write(f"Total assembly length: {m1.get('total_length', 0):,} bp\n")
    out.write(f"N50: {m1.get('n50', 0):,} bp\n\n")

    # assembly 2 metrics
    out.write("=" * 70 + "\n")
    out.write("ASSEMBLY 2: Paired-end + merged reads\n")
    out.write("=" * 70 + "\n")
    out.write(f"Number of contigs: {m2.get('contigs', 0):,}\n")
    out.write(f"Total assembly length: {m2.get('total_length', 0):,} bp\n")
    out.write(f"N50: {m2.get('n50', 0):,} bp\n\n")

    # comparison deltas
    out.write("=" * 70 + "\n")
    out.write("COMPARISON\n")
    out.write("=" * 70 + "\n")

    if m1 and m2:
        contig_diff = m2.get('contigs', 0)      - m1.get('contigs', 0)
        length_diff = m2.get('total_length', 0) - m1.get('total_length', 0)
        n50_diff    = m2.get('n50', 0)          - m1.get('n50', 0)

        pct = lambda diff, base: f" ({diff/base*100:+.1f}%)" if base > 0 else ""

        out.write(f"Δ Contigs: {contig_diff:+,}{pct(contig_diff, m1.get('contigs', 1))}\n")
        out.write(f"Δ Length: {length_diff:+,} bp{pct(length_diff, m1.get('total_length', 1))}\n")
        out.write(f"Δ N50: {n50_diff:+,} bp{pct(n50_diff, m1.get('n50', 1))}\n\n")

    # best assembly selection
    out.write("=" * 70 + "\n")
    out.write(f"BEST ASSEMBLY: Assembly {best_assembly}\n")
    out.write("=" * 70 + "\n")
    if reasons:
        for reason in reasons:
            out.write(f"- {reason}\n")
    else:
        out.write("- Similar quality metrics\n")
    out.write(f"\nSelected assembly for polishing: {sample}_best.fasta\n")
    out.write(f"Filtering will run after polishing (min_length={min_length} bp, min_coverage={min_coverage} X)\n")

# write best assembly path for next step
with open("best_assembly.txt", "w") as f:
    f.write(best_contigs)

print(f"\nSummary written to {sample}_assembly_summary.txt")
print(f"Best assembly path written to best_assembly.txt")
