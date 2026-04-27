#!/usr/bin/env python3

import sys
import re

def read_fasta(path):
    seqs = {}
    name = None
    buf = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if name:
        seqs[name] = "".join(buf)
    return seqs

def extract_ori_id(polished_header):

    m = re.search(r'_ori_([^_]+_[0-9]+)_len_', polished_header)
    if m:
        return m.group(1)

    if "_polish_mean_depth_" in polished_header:
        return polished_header.split("_polish_mean_depth_")[0]

    if polished_header.endswith("_pilon"):
        return polished_header[:-6]

    return polished_header

def build_before_index(before):
    """
    Build a case-insensitive lookup index from the before FASTA.
    Maps simplified keys like 'node_1' to full header names.
    
    SPAdes headers: NODE_1_length_X_cov_Y
    After tolower + truncation in awk: node_1
    """
    index = {}
    for name in before:
        m = re.match(r'(NODE_[0-9]+)', name, re.IGNORECASE)
        if m:
            key = m.group(1).lower().replace('node_', 'node_')
            index[key] = name
        index[name.lower()] = name
    return index

def count_diff(s1, s2):

    if not s1:
        return None, None, None, None

    min_len = min(len(s1), len(s2))
    subs = sum(1 for a, b in zip(s1, s2) if a != b)
    len_diff = len(s2) - len(s1)
    ins  = max(0,  len_diff)
    dels = max(0, -len_diff)
    total = subs + ins + dels
    return total, subs, ins, dels

def main():
    if len(sys.argv) < 3:
        print("Usage: compare_polished.py <before.fasta> <after.fasta>", file=sys.stderr)
        sys.exit(1)

    before_path = sys.argv[1]
    after_path  = sys.argv[2]

    # Optional: sample_id as 3rd arg for output filename
    sample_id = sys.argv[3] if len(sys.argv) > 3 else None
    out_path = f"{sample_id}_polished_status.txt" if sample_id else "polished_status.txt"

    before = read_fasta(before_path)
    after  = read_fasta(after_path)

    # Build case-insensitive index for before sequences
    before_index = build_before_index(before)

    matched   = 0
    unmatched = 0

    with open(out_path, "w") as out:
        out.write("contig\ttotal_changes\tsubs\tins\tdels\tstatus\n")
        for contig, seq in after.items():
            ori_key = extract_ori_id(contig)

            # Try direct match first, then index lookup
            before_seq = before.get(ori_key) or \
                         before.get(ori_key.upper()) or \
                         before.get(before_index.get(ori_key.lower(), ""), "")

            total, subs, ins, dels = count_diff(before_seq, seq)

            if total is None:
                status = "NO_MATCH"
                out.write(f"{contig}\tNA\tNA\tNA\tNA\t{status}\n")
                unmatched += 1
            else:
                status = "CHANGED" if total > 0 else "UNCHANGED"
                out.write(f"{contig}\t{total}\t{subs}\t{ins}\t{dels}\t{status}\n")
                matched += 1

    print(f"Compared {matched} matched contigs, {unmatched} unmatched.", file=sys.stderr)

if __name__ == "__main__":
    main()