#!/usr/bin/env python3
import os
import sys
import glob
import argparse
import subprocess
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from itertools import combinations
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

class QC:
    samples_processed = 0
    hits_per_db = {'ncbi': 0, 'card': 0, 'resfinder': 0, 'abritamr': 0, 'argannot': 0}
    seqs_extracted = 0
    missing_genomes = 0
    missing_contigs = 0
    overlap_candidates = 0
    blast_comparisons = 0
    sequence_matches = 0
    different_sequences = 0
    consensus_genes = 0
    single_db_genes = 0
    multi_db_genes = 0

def read_input_tsv(filepath, db_name):
    if not filepath or not os.path.exists(filepath):
        return pd.DataFrame()

    try:
        df_raw = pd.read_csv(filepath, sep='\t', dtype=str)
    except Exception as e:
        logging.error(f"Failed to read {filepath}: {e}")
        return pd.DataFrame()

    records = []
    for _, row in df_raw.iterrows():
        try:
            sample = str(row.iloc[0]).strip()
            
            for ext in ['.fasta', '.fa', '.fna', '.faa']:
                if sample.lower().endswith(ext):
                    sample = sample[:-len(ext)]
                    break

            if db_name == 'abritamr':
                # Columns: 3=Contig, 4=Start, 5=Stop, 6=Strand, 7=Gene
                contig = row.iloc[2]
                start, stop = int(row.iloc[3]), int(row.iloc[4])
                strand = row.iloc[5]
                gene = row.iloc[6]
            else:
                # Columns: 2=Contig, 3=Start, 4=Stop, 5=Strand, 6=Gene
                contig = row.iloc[1]
                start, stop = int(row.iloc[2]), int(row.iloc[3])
                strand = row.iloc[4]
                gene = row.iloc[5]
            
            records.append({
                'Sample_ID': sample,
                'Database': db_name,
                'Contig': contig,
                'Start': start,
                'Stop': stop,
                'Strand': strand,
                'Gene_Name': gene
            })
            QC.hits_per_db[db_name] += 1
        except (ValueError, IndexError):
            continue

    return pd.DataFrame(records)

def find_genome(sample_id, ref_dir):
    """Finds the corresponding FASTA file for a sample."""
    extensions = ['.fasta', '.fa', '.fna', '.faa']
    for ext in extensions:
        path = os.path.join(ref_dir, f"{sample_id}{ext}")
        if os.path.exists(path):
            return path
    return None

_genome_cache = {}

def load_genome(genome_path):
    if genome_path not in _genome_cache:
        records = {}
        try:
            for rec in SeqIO.parse(genome_path, "fasta"):
                records.setdefault(rec.id, rec)
        except Exception as e:
            logging.warning(f"Failed to parse genome {genome_path}: {e}")
        _genome_cache[genome_path] = records
    return _genome_cache[genome_path]

def extract_sequence(record, genome_path):
    contig = record['Contig']
    start = min(record['Start'], record['Stop'])
    stop = max(record['Start'], record['Stop'])
    strand = record['Strand']

    try:
        genome = load_genome(genome_path)
        rec = genome.get(str(contig))
        if rec is None:
            QC.missing_contigs += 1
            return None

        seq_str = str(rec.seq[start - 1:stop])

        if not seq_str:
            QC.missing_contigs += 1
            return None

        if strand == '-':
            seq_str = str(Seq(seq_str).reverse_complement())

        return seq_str.upper()

    except Exception as e:
        logging.error(f"Sequence extraction failed for {contig}:{start}-{stop} in {genome_path}: {e}")
        QC.missing_contigs += 1
        return None

def is_overlapping(a, b):
    if a['Contig'] != b['Contig']:
        return False
    # Max of starts <= Min of stops
    return max(a['Start'], b['Start']) <= min(a['Stop'], b['Stop'])

def get_connected_components(nodes, edges):
    adj = {n: set() for n in nodes}
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    
    visited = set()
    components = []
    for n in nodes:
        if n not in visited:
            comp = []
            queue = [n]
            while queue:
                curr = queue.pop(0)
                if curr not in visited:
                    visited.add(curr)
                    comp.append(curr)
                    queue.extend(list(adj[curr] - visited))
            components.append(comp)
    return components

def main():
    parser = argparse.ArgumentParser(description="Cross-check and consolidate AMR genes using flexible BLAST thresholds.")
    parser.add_argument('--abritamr', help="abritamr output TSV")
    parser.add_argument('--argannot', help="ARG-ANNOT output TSV")
    parser.add_argument('--card', help="CARD output TSV")
    parser.add_argument('--ncbi', help="NCBI output TSV")
    parser.add_argument('--resfinder', help="ResFinder output TSV")
    parser.add_argument('--reference-dir', required=True, help="Directory containing reference genomes")
    parser.add_argument('--outdir', required=True, help="Output directory")
    
    # New flexible threshold arguments
    parser.add_argument('--coverage', type=float, default=90.0, help="Minimum sequence coverage percentage (default: 90.0)")
    parser.add_argument('--identity', type=float, default=100.0, help="Minimum sequence identity percentage (default: 100.0)")
    parser.add_argument('--threads', type=int, default=8, help="Number of threads for BLAST (default: 8)")
    
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    logging.info(f"Reading input files... (Thresholds: {args.identity}% ID, {args.coverage}% Cov)")
    df_list = []
    for path, db in [(args.ncbi, 'ncbi'), (args.card, 'card'), (args.resfinder, 'resfinder'), 
                     (args.abritamr, 'abritamr'), (args.argannot, 'argannot')]:
        df = read_input_tsv(path, db)
        if not df.empty:
            df_list.append(df)
            
    if not df_list:
        logging.error("No valid inputs parsed. Exiting.")
        sys.exit(1)
        
    all_genes_df = pd.concat(df_list, ignore_index=True)
    all_samples = all_genes_df['Sample_ID'].unique()

    # Open output files
    fasta_out_path = os.path.join(args.outdir, 'amr_sequences.fasta')
    blast_out_path = os.path.join(args.outdir, 'blast_results.tsv')
    crosscheck_out_path = os.path.join(args.outdir, 'amr_crosscheck.tsv')
    consensus_out_path = os.path.join(args.outdir, 'amr_consensus.tsv')
    
    fasta_out = open(fasta_out_path, 'w')
    
    blast_records = []
    crosscheck_records = []
    consensus_records = []
    consensus_id_counter = 1
    
    db_priority = ['ncbi', 'card', 'resfinder', 'abritamr', 'argannot']

    sample_seqs_map = {}
    sample_order = []

    for sample in all_samples:
        QC.samples_processed += 1
        sample_df = all_genes_df[all_genes_df['Sample_ID'] == sample].copy()
        
        genome_path = find_genome(sample, args.reference_dir)
        if not genome_path:
            logging.warning(f"Missing reference genome for sample {sample}. Skipping.")
            QC.missing_genomes += 1
            continue
            
        sample_seqs = {}
        fasta_ids_to_row = {}
        
        for idx, row in sample_df.iterrows():
            seq = extract_sequence(row, genome_path)
            if seq:
                QC.seqs_extracted += 1
                fasta_id = f"{row['Sample_ID']}|{row['Database']}|{row['Gene_Name']}|{row['Contig']}|{row['Start']}|{row['Stop']}|{row['Strand']}"
                if fasta_id in sample_seqs:
                    fasta_id = f"{fasta_id}|{idx}"
                
                sample_seqs[fasta_id] = seq
                fasta_ids_to_row[fasta_id] = row
                
                fasta_out.write(f">{fasta_id}\n{seq}\n")
        
        if not sample_seqs:
            continue
            
        candidates = []
        fasta_id_list = list(sample_seqs.keys())
        for a_id, b_id in combinations(fasta_id_list, 2):
            row_a = fasta_ids_to_row[a_id]
            row_b = fasta_ids_to_row[b_id]
            if is_overlapping(row_a, row_b):
                candidates.append((a_id, b_id))
                QC.overlap_candidates += 1
                
        sample_order.append(sample)
        sample_seqs_map[sample] = (sample_seqs, fasta_ids_to_row, candidates, fasta_id_list)

    fasta_out.close()

    # ---- Per-sample consensus using per-sample BLAST (keeps E-values db-local) ----
    for sample in sample_order:
        sample_seqs, fasta_ids_to_row, candidates, fasta_id_list = sample_seqs_map[sample]

        match_edges = []
        blast_dict = {}
        if candidates:
            tmp_fasta = os.path.join(args.outdir, f"tmp_{sample}.fasta")
            tmp_db = os.path.join(args.outdir, f"tmp_{sample}_db")
            tmp_blast = os.path.join(args.outdir, f"tmp_{sample}_blast.tsv")

            with open(tmp_fasta, 'w') as tf:
                for fid, seq in sample_seqs.items():
                    tf.write(f">{fid}\n{seq}\n")

            subprocess.run(["makeblastdb", "-in", tmp_fasta, "-dbtype", "nucl", "-out", tmp_db], capture_output=True)

            blast_cmd = [
                "blastn", "-query", tmp_fasta, "-db", tmp_db, "-num_threads", str(args.threads),
                "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen",
                "-out", tmp_blast
            ]
            subprocess.run(blast_cmd, capture_output=True)

            if os.path.exists(tmp_blast):
                with open(tmp_blast, 'r') as bf:
                    for line in bf:
                        parts = line.strip().split('\t')
                        if len(parts) >= 14:
                            qid, sid = parts[0], parts[1]
                            if qid != sid:
                                QC.blast_comparisons += 1
                                pident, length = float(parts[2]), int(parts[3])
                                evalue, bitscore = float(parts[10]), float(parts[11])
                                qlen, slen = int(parts[12]), int(parts[13])

                                qcov = (length / qlen) * 100
                                scov = (length / slen) * 100

                                blast_records.append({
                                    'Sample_ID': sample, 'query_ID': qid, 'subject_ID': sid,
                                    'identity': pident, 'alignment_length': length,
                                    'query_coverage': qcov, 'subject_coverage': scov,
                                    'evalue': evalue, 'bit_score': bitscore
                                })

                                blast_dict[(qid, sid)] = {
                                    'pident': pident, 'length': length,
                                    'qcov': qcov, 'scov': scov,
                                    'evalue': evalue, 'bitscore': bitscore
                                }

        for a_id, b_id in candidates:
            a_row = fasta_ids_to_row[a_id]
            b_row = fasta_ids_to_row[b_id]
            seq_a = sample_seqs[a_id]
            seq_b = sample_seqs[b_id]
            
            b_res = blast_dict.get((a_id, b_id)) or blast_dict.get((b_id, a_id))
            
            is_match = False
            match_decision = "Different_gene"
            
            # Rule check: BLAST metrics >= user thresholds
            if b_res and b_res['pident'] >= args.identity and b_res['scov'] >= args.coverage:
                is_match = True
                match_decision = "Same_gene"
                match_edges.append((a_id, b_id))
                QC.sequence_matches += 1
            else:
                QC.different_sequences += 1
                
            crosscheck_records.append({
                'Sample_ID': sample,
                'Query_Database': a_row['Database'], 'Query_Gene': a_row['Gene_Name'],
                'Subject_Database': b_row['Database'], 'Subject_Gene': b_row['Gene_Name'],
                'Query_Contig': a_row['Contig'], 'Subject_Contig': b_row['Contig'],
                'Query_Start': a_row['Start'], 'Query_Stop': a_row['Stop'],
                'Subject_Start': b_row['Start'], 'Subject_Stop': b_row['Stop'],
                'Query_Length': len(seq_a), 'Subject_Length': len(seq_b),
                'Alignment_Length': b_res['length'] if b_res else 0,
                'Percent_Identity': b_res['pident'] if b_res else 0,
                'Query_Coverage': b_res['qcov'] if b_res else 0,
                'Subject_Coverage': b_res['scov'] if b_res else 0,
                'Evalue': b_res['evalue'] if b_res else "NA",
                'Bit_Score': b_res['bitscore'] if b_res else "NA",
                'Sequence_Match': "TRUE" if is_match else "FALSE",
                'Match_Decision': match_decision
            })

        if candidates:
            for tmp_file in glob.glob(os.path.join(args.outdir, f"tmp_{sample}*")):
                os.remove(tmp_file)

        components = get_connected_components(fasta_id_list, match_edges)
        
        for comp in components:
            QC.consensus_genes += 1
            if len(comp) > 1:
                QC.multi_db_genes += 1
                match_status = "Consensus_match"
            else:
                QC.single_db_genes += 1
                match_status = "Single_database"
                
            comp_rows = [fasta_ids_to_row[fid] for fid in comp]
            comp_rows.sort(key=lambda x: db_priority.index(x['Database']))
            rep = comp_rows[0]
            
            db_genes = {db: '' for db in db_priority}
            for row in comp_rows:
                db_genes[row['Database']] = row['Gene_Name']
                
            consensus_records.append({
                'Sample_ID': sample,
                'Consensus_ID': f"AMR_{consensus_id_counter:06d}",
                'Contig': rep['Contig'],
                'Start': rep['Start'],
                'Stop': rep['Stop'],
                'Strand': rep['Strand'],
                'Consensus_Gene': rep['Gene_Name'],
                'Supporting_Databases': ",".join([r['Database'] for r in comp_rows]),
                'Number_of_Databases': len(comp_rows),
                'Sequence_Length': len(sample_seqs[comp[0]]),
                'Sequence': sample_seqs[comp[0]],
                'NCBI_Gene': db_genes['ncbi'],
                'CARD_Gene': db_genes['card'],
                'ResFinder_Gene': db_genes['resfinder'],
                'AbritAMR_Gene': db_genes['abritamr'],
                'ARGANNOT_Gene': db_genes['argannot'],
                'Match_Status': match_status
            })
            consensus_id_counter += 1

    logging.info("Writing output files...")
    if consensus_records:
        pd.DataFrame(consensus_records).to_csv(consensus_out_path, sep='\t', index=False)
    if crosscheck_records:
        pd.DataFrame(crosscheck_records).to_csv(crosscheck_out_path, sep='\t', index=False)
    if blast_records:
        pd.DataFrame(blast_records).to_csv(blast_out_path, sep='\t', index=False)

    print("\n" + "="*50)
    print(" AMR CROSSCHECK QUALITY CONTROL REPORT")
    print("="*50)
    print(f"Samples processed:                  {QC.samples_processed}")
    print(f"Total sequences extracted:          {QC.seqs_extracted}")
    print("\n--- AMR hits per database ---")
    for db, count in QC.hits_per_db.items():
        print(f"  {db.upper().ljust(15)}: {count}")
    print("\n--- Errors / Missing ---")
    print(f"Missing reference genomes:          {QC.missing_genomes}")
    print(f"Missing contigs (in fasta):         {QC.missing_contigs}")
    print("\n--- Candidate Comparisons ---")
    print(f"Coordinate-overlap candidates:      {QC.overlap_candidates}")
    print(f"BLAST comparisons generated:        {QC.blast_comparisons}")
    print(f"Sequence matches (>= {args.identity}% ID, >= {args.coverage}% Cov): {QC.sequence_matches}")
    print(f"Different sequences identified:     {QC.different_sequences}")
    print("\n--- Final Consensus ---")
    print(f"Total consensus genes generated:    {QC.consensus_genes}")
    print(f"Single-database genes:              {QC.single_db_genes}")
    print(f"Multi-database consensus genes:     {QC.multi_db_genes}")
    print("="*50 + "\n")
    logging.info(f"Pipeline completed successfully. Results saved in '{args.outdir}'.")

if __name__ == '__main__':
    main()