#!/usr/bin/env python3

import argparse
import glob
import json
import os
import re
import sys

import pandas as pd

DB_META = {
    1: {"name": "NCBI AMRFinderPlus",          "key": "ncbi",          "sup": "\u00b9", "color": "#1976d2"},
    2: {"name": "CARD",                        "key": "card",          "sup": "\u00b2", "color": "#e53935"},
    3: {"name": "ResFinder",                   "key": "resfinder",     "sup": "\u00b3", "color": "#0a9e74"},
    4: {"name": "ARG-ANNOT",                   "key": "argannot",      "sup": "\u2074", "color": "#7b1fa2"},
    5: {"name": "AMRFinderPlus (abritAMR)",    "key": "amrfinderplus", "sup": "\u2075", "color": "#00838f"},
}
# Unicode superscripts:  ¹ ² ³ ⁴ ⁵

def clean_abritamr_gene(gene: str) -> str:
    gene = gene.strip()
    gene = re.sub(r"\*+\s*$", "", gene)
    gene = re.sub(r"\(\*+\)\s*$", "", gene)
    gene = re.sub(r"\([A-Z][a-z]{2,}\w*\)\s*$", "", gene)
    return gene.strip()


_ENZYME_NUM_PAREN_RE = re.compile(r"\((\d+'{0,2})\)")


def strip_enzyme_number_parens(gene: str) -> str:
    """Remove parentheses that wrap a bare enzyme/regiochemistry number.

    Different databases report the same gene with or without parentheses
    around its numeric designation, e.g. 'aac(3)-IId' (parenthesised) vs
    'aac3-IId' (bare). Both forms should collapse onto the same gene during
    comparison, so this is applied only when building the matching key in
    normalize_gene() -- deliberately NOT in prepare_gene_for_display(). That
    way the original spelling (often the parenthesised, standard-nomenclature
    form) is still available as a display candidate, and best_canonical()
    can bring the parentheses back for the report even though matching
    ignored them. Parentheses are stripped only when their content is purely
    numeric, optionally followed by prime marks (e.g. "(6')", "(2'')").
    Parentheses wrapping letters -- single-letter alleles like 'erm(A)', or
    the species/leading-class annotations handled elsewhere for ARG-ANNOT
    and abritAMR -- are left untouched here.
    """
    return _ENZYME_NUM_PAREN_RE.sub(r"\1", gene)


def prepare_gene_for_display(gene: str, db_id=None) -> str:
    g = gene.strip()

    if db_id == 4:
        g = re.sub(r"^\([^)]*\)", "", g).strip()
    elif db_id == 3:
        g = re.sub(r"_(\d+)$", "", g).strip()
    elif db_id == 5:
        g = clean_abritamr_gene(g)

    return g.strip()

def normalize_qnr(g: str) -> str:
    return re.sub(r'^(qnr)-([a-z]\d+)$', r'\1\2', g)

def normalize_gene(gene: str, db_id=None) -> str:
    g = prepare_gene_for_display(gene, db_id)
    if g.lower().startswith("bla"):
        g = g[3:]
    g = strip_enzyme_number_parens(g)
    g = re.sub(r"\(([A-Za-z])\)", r"\1", g)
    g = g.lower()
    g = g.replace("_", "-")
    g = normalize_qnr(g) 
    g = re.sub(r"-+", "-", g)
    return g.strip("-")


_QNR_REPORT_RE = re.compile(r"^qnr[-_]?([a-z])[-_]?(\d+)$", re.IGNORECASE)


def canonicalize_qnr_gene(gene: str) -> str:
    m = _QNR_REPORT_RE.match(gene.strip())
    if m:
        return f"qnr{m.group(1).upper()}{m.group(2)}"
    return gene


def format_gene_for_report(gene: str, had_bla: bool = False) -> str:
    gene = canonicalize_qnr_gene(gene)
    if had_bla and not gene.lower().startswith("bla"):
        gene = f"bla{gene}"
    return gene


def best_canonical(raw_names) -> str:
    def rank(n):
        has_bla     = int(n.lower().startswith("bla"))
        no_paren    = int("(" not in n)
        upper_start = int(bool(n) and n[0].isupper() and not n.lower().startswith("bla"))
        return (has_bla, no_paren, upper_start, len(n), n.lower())

    return min(raw_names, key=rank)


def make_entry_payload(gene: str, db_ids: list) -> dict:
    db_ids = sorted(db_ids)
    sups = "".join(DB_META[i]["sup"] for i in db_ids)
    return {
        "gene": gene,
        "display": f"{gene}{sups}",
        "db_ids": db_ids,
        "db_names": [DB_META[i]["name"] for i in db_ids],
        "db_count": len(db_ids),
    }

def parse_abricate_tsv(tsv_path: str) -> dict:
    sample_genes: dict = {}
    try:
        with open(tsv_path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.rstrip("\n\r")
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 6:
                    continue
                raw = os.path.basename(parts[0])
                for ext in (".fasta", ".fna", ".fa", ".fas"):
                    if raw.lower().endswith(ext):
                        raw = raw[: -len(ext)]
                        break
                sample_id = raw
                gene = parts[5].strip()
                if not gene:
                    continue
                sample_genes.setdefault(sample_id, []).append(gene)
    except OSError as exc:
        print(f"  Warning: cannot open {tsv_path}: {exc}", file=sys.stderr)
    return sample_genes


def parse_abricate_dir(abricate_dir: str) -> dict:
    results: dict = {}
    for db_id, meta in DB_META.items():
        if db_id == 5:
            continue  # AMRFinderPlus handled separately
        key = meta["key"]
        candidates = (
            glob.glob(os.path.join(abricate_dir, f"{key}.tsv"))
            + glob.glob(os.path.join(abricate_dir, f"{key.upper()}.tsv"))
            + glob.glob(os.path.join(abricate_dir, f"*{key}*.tsv"))
        )

        seen: set = set()
        candidates = [c for c in candidates if not (c in seen or seen.add(c))]

        found = next((c for c in candidates if os.path.isfile(c)), None)
        if found:
            genes = parse_abricate_tsv(found)
            results[db_id] = genes
            total = sum(len(v) for v in genes.values())
            print(
                f"  Parsed {len(genes):>4} samples, {total:>6} hits  "
                f"— {meta['name']:12s} ({os.path.basename(found)})"
            )
        else:
            results[db_id] = {}
            print(f"  Warning: no TSV found for '{key}' in {abricate_dir}")
    return results


def parse_abritamr_tsv(tsv_path: str) -> dict:

    sample_genes: dict = {}
    if not os.path.isfile(tsv_path):
        print(f"  Warning: AbritAMR file not found: {tsv_path}", file=sys.stderr)
        return sample_genes
    try:
        df = pd.read_csv(tsv_path, sep="\t", dtype=str).fillna("")
    except Exception as exc:
        print(f"  Warning: cannot parse {tsv_path}: {exc}", file=sys.stderr)
        return sample_genes

    if df.empty:
        return sample_genes

    isolate_col = df.columns[0]
    class_cols  = list(df.columns[1:])

    for _, row in df.iterrows():
        sample_id = str(row[isolate_col]).strip()
        genes: list = []
        for col in class_cols:
            val = str(row[col]).strip()
            if not val or val.lower() == "nan":
                continue
            for raw_g in val.split(","):
                g = prepare_gene_for_display(raw_g.strip(), 5)
                if g:
                    genes.append(g)
        if genes:
            sample_genes[sample_id] = genes

    total = sum(len(v) for v in sample_genes.values())
    print(
        f"  Parsed {len(sample_genes):>4} samples, {total:>6} hits  "
        f"— AMRFinderPlus  ({os.path.basename(tsv_path)})"
    )
    return sample_genes

def build_consensus(all_db_genes: dict) -> dict:

    all_samples: set = set()
    for db_genes in all_db_genes.values():
        all_samples.update(db_genes.keys())

    consensus: dict = {}
    for sample_id in sorted(all_samples):
        gene_map: dict = {}  # norm_key → {raw_names: set, db_ids: set}
        for db_id, db_genes in all_db_genes.items():
            for gene in db_genes.get(sample_id, []):
                display_gene = prepare_gene_for_display(gene, db_id)
                if not display_gene:
                    continue
                key = normalize_gene(gene, db_id)
                if not key:
                    key = display_gene.lower()
                    print(
                        f"  Warning: gene '{gene}' (sample {sample_id}) normalized to an "
                        f"empty matching key; reporting as '{display_gene}' instead of dropping it",
                        file=sys.stderr,
                    )
                if key not in gene_map:
                    gene_map[key] = {"display_names": set(), "db_ids": set()}
                gene_map[key]["display_names"].add(display_gene)
                gene_map[key]["db_ids"].add(db_id)

        entries = []
        for _key, info in sorted(gene_map.items()):
            canonical = best_canonical(info["display_names"])
            had_bla = any(n.lower().startswith("bla") for n in info["display_names"])
            canonical = format_gene_for_report(canonical, had_bla)
            entries.append(make_entry_payload(canonical, info["db_ids"]))

        entries.sort(key=lambda t: t["gene"].lower())
        consensus[sample_id] = {
            "entries": entries,
            "multi_entries": [entry for entry in entries if entry["db_count"] >= 2],
            "single_entries": [entry for entry in entries if entry["db_count"] == 1],
            "total_count": len(entries),
        }
    return consensus


def write_tsv(consensus: dict, tsv_path: str) -> None:
    header = [
        "Sample_ID", "Gene", "DB_Support_Count",
        *[meta["name"] for _, meta in DB_META.items()],
        "Display",
    ]
    rows = []
    for sample_id, sample_data in sorted(consensus.items()):
        for entry in sample_data["entries"]:
            db_ids = entry["db_ids"]
            rows.append({
                "Sample_ID":       sample_id,
                "Gene":            entry["gene"],
                "DB_Support_Count": len(db_ids),
                **{
                    meta["name"]: "1" if db_id in db_ids else "0"
                    for db_id, meta in DB_META.items()
                },
                "Display":         entry["display"],
            })
    pd.DataFrame(rows, columns=header).to_csv(tsv_path, sep="\t", index=False)
    print(f"  TSV summary : {tsv_path}")


def _build_js_data(consensus: dict):
    """Return (json_str, db_gene_counts) for injection into HTML."""
    db_gene_counts = {db_id: 0 for db_id in DB_META}
    rows = []

    for sample_id, sample_data in sorted(consensus.items()):
        entries = sample_data["entries"]
        multi_entries = sample_data["multi_entries"]
        single_entries = sample_data["single_entries"]

        for entry in entries:
            for db_id in entry["db_ids"]:
                db_gene_counts[db_id] += 1

        rows.append({
            "sample_id":          sample_id,
            "total_gene_count":   sample_data["total_count"],
            "multi_gene_count":   len(multi_entries),
            "single_gene_count":  len(single_entries),
            "entries":            entries,
            "multi_entries":      multi_entries,
            "single_entries":     single_entries,
            "total_gene_str":     ", ".join(entry["display"] for entry in entries),
            "multi_gene_str":     ", ".join(entry["display"] for entry in multi_entries),
            "single_gene_str":    ", ".join(entry["display"] for entry in single_entries),
        })

    return json.dumps(rows), db_gene_counts


def generate_html(consensus: dict, output_path: str, minid_genes: float, mincov_genes: float) -> None:

    all_samples = sorted(consensus.keys())
    data_json, db_gene_counts = _build_js_data(consensus)
    db_meta_json = json.dumps({str(k): v for k, v in DB_META.items()})

    # Global statistics
    total_samples = len(all_samples)
    all_unique_genes: set = set()
    for sample_data in consensus.values():
        for entry in sample_data["entries"]:
            all_unique_genes.add(entry["gene"].lower())
    total_unique_genes = len(all_unique_genes)
    avg_genes = (
        round(sum(len(v) for v in consensus.values()) / total_samples, 1)
        if total_samples else 0
    )

    # HTML building blocks
    sample_options = "".join(
        f'<label class="multi-option">'
        f'<input type="checkbox" value="{s}" onchange="onMultiChange()"/>{s}'
        f"</label>"
        for s in all_samples
    )

    db_stat_cards = "".join(
        f'<label class="stat-card db-stat db-toggle" data-db-id="{db_id}">'
        f'<input type="checkbox" class="db-toggle-checkbox" data-db-id="{db_id}" '
        f'checked onchange="onDbToggle({db_id})"/>'
        f'<div class="db-toggle-body">'
        f'<div class="stat-label">{m["sup"]} {m["name"]}</div>'
        f'<div class="stat-value" style="color:{m["color"]}">{db_gene_counts[db_id]}</div>'
        f'</div>'
        f'</label>'
        for db_id, m in DB_META.items()
    )

    # Full HTML 
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1.0"/>
<title>BASIL // AMR Consensus Report</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/xlsx/0.18.5/xlsx.full.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/jspdf/2.5.1/jspdf.umd.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/jspdf-autotable/3.8.2/jspdf.plugin.autotable.min.js"></script>
<style>
  :root {{
    --bg:        #f8f9fa;
    --surface:   #ffffff;
    --surface2:  #f0f2f5;
    --border:    #dee2e6;
    --accent:    #0a9e74;
    --accent2:   #e53935;
    --accent3:   #e67e00;
    --text:      #1a1f2e;
    --muted:     #6c757d;
    --font-main: Calibri, Arial, sans-serif;
  }}
  *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ background: var(--bg); color: var(--text); font-family: var(--font-main); min-height: 100vh; }}

  /* HEADER */
  header {{
    background: var(--surface); border-bottom: 2px solid var(--accent);
    box-shadow: 0 2px 8px rgba(0,0,0,0.06);
    padding: 24px 40px; display: flex; align-items: center;
    justify-content: space-between; gap: 16px; flex-wrap: wrap;
  }}
  .logo {{ font-size: 1.4rem; font-weight: 600; color: var(--accent);
           letter-spacing: 2px; text-transform: uppercase; }}
  .logo span {{ color: inherit; }}

  /* EXPORT */
  .export-wrap {{ position: relative; }}
  .export-btn {{
    background: var(--accent); color: #fff; border: none; padding: 10px 22px;
    border-radius: 6px; font-family: var(--font-main); font-size: 0.75rem;
    font-weight: 700; letter-spacing: 1px; cursor: pointer; transition: opacity 0.2s;
  }}
  .export-btn:hover {{ opacity: 0.85; }}
  .export-menu {{
    display: none; position: absolute; top: calc(100% + 6px); right: 0;
    background: var(--surface2); border: 1px solid var(--border); border-radius: 8px;
    overflow: hidden; min-width: 190px; z-index: 100;
    box-shadow: 0 8px 24px rgba(0,0,0,0.4);
  }}
  .export-wrap:hover .export-menu,
  .export-wrap:focus-within .export-menu {{ display: block; }}
  .export-menu button {{
    display: flex; align-items: center; gap: 10px; width: 100%;
    background: transparent; border: none; color: var(--text); padding: 11px 18px;
    font-family: var(--font-main); font-size: 0.72rem; letter-spacing: 0.5px;
    cursor: pointer; transition: background 0.15s; text-align: left;
  }}
  .export-menu button:hover {{ background: var(--surface); color: var(--accent); }}
  .export-menu .divider {{ height: 1px; background: var(--border); margin: 4px 0; }}

  /* STATS */
  .stats-section {{ padding: 24px 40px 0; }}
  .stats-bar {{ display: flex; gap: 14px; flex-wrap: wrap; margin-bottom: 12px; }}
  .stat-card {{
    background: var(--surface); border: 1px solid var(--border);
    border-radius: 10px; padding: 18px 24px; flex: 1; min-width: 140px;
    box-shadow: 0 1px 4px rgba(0,0,0,0.04);
  }}
  .stat-card.db-stat {{ padding: 12px 16px; min-width: 110px; }}
  .stat-label {{ font-size: 0.7rem; color: var(--muted); text-transform: uppercase;
                 letter-spacing: 1.5px; margin-bottom: 8px; }}
  .stat-value {{ font-size: 2rem; font-weight: 600; color: var(--accent); }}
  .db-stat .stat-value {{ font-size: 1.5rem; }}

  /* DATABASE SELECTION */
  .mini-btn {{
    background: transparent; border: 1px solid var(--border); color: var(--muted);
    padding: 4px 12px; border-radius: 14px; cursor: pointer;
    font-family: var(--font-main); font-size: 0.68rem; letter-spacing: 0.5px;
    transition: border-color 0.2s, color 0.2s;
  }}
  .mini-btn:hover {{ border-color: var(--accent); color: var(--accent); }}
  .db-stat.db-toggle {{
    display: flex; align-items: flex-start; gap: 10px;
    cursor: pointer; user-select: none;
    transition: opacity 0.2s, border-color 0.2s;
  }}
  .db-stat.db-toggle:hover {{ border-color: var(--accent); }}
  .db-stat.db-toggle.db-off {{ opacity: 0.35; }}
  .db-toggle-checkbox {{
    margin-top: 2px; accent-color: var(--accent); width: 15px; height: 15px; cursor: pointer;
  }}
  .db-toggle-body {{ display: flex; flex-direction: column; }}

  /* LEGEND */
  .legend-bar {{
    background: var(--surface); border: 1px solid var(--border); border-radius: 10px;
    padding: 12px 20px; display: flex; align-items: center; gap: 10px;
    flex-wrap: wrap; margin-bottom: 20px; box-shadow: 0 1px 4px rgba(0,0,0,0.04);
  }}
  .legend-title {{ font-size: 0.68rem; color: var(--muted); text-transform: uppercase;
                   letter-spacing: 1.5px; font-weight: 600; white-space: nowrap; }}

  /* CONTROLS */
  .controls {{
    padding: 0 40px 8px; display: flex; gap: 12px; flex-wrap: wrap; align-items: center;
  }}
  .filter-group {{ display: flex; flex-direction: column; gap: 4px; }}
  .filter-label {{ font-size: 0.65rem; color: var(--muted);
                   text-transform: uppercase; letter-spacing: 1px; }}
  select, input[type=text] {{
    background: var(--surface); border: 1px solid var(--border); color: var(--text);
    padding: 8px 14px; border-radius: 6px; font-family: var(--font-main);
    font-size: 0.875rem; outline: none; transition: border-color 0.2s; min-width: 180px;
  }}
  select:focus, input:focus {{ border-color: var(--accent); }}

  /* MULTI-SELECT */
  .multi-wrap {{ position: relative; min-width: 240px; }}
  .multi-trigger {{
    background: var(--surface); border: 1px solid var(--border); color: var(--text);
    padding: 8px 14px; border-radius: 6px; font-family: var(--font-main);
    font-size: 0.875rem; cursor: pointer; display: flex;
    justify-content: space-between; align-items: center; gap: 8px;
    user-select: none; transition: border-color 0.2s;
  }}
  .multi-trigger:hover, .multi-trigger.open {{ border-color: var(--accent); }}
  .multi-trigger .arrow {{ font-size: 0.6rem; color: var(--muted); transition: transform 0.2s; }}
  .multi-trigger.open .arrow {{ transform: rotate(180deg); }}
  .multi-dropdown {{
    display: none; position: absolute; top: calc(100% + 4px); left: 0;
    background: var(--surface2); border: 1px solid var(--border); border-radius: 8px;
    z-index: 200; min-width: 100%; max-height: 280px; overflow-y: auto;
    box-shadow: 0 8px 24px rgba(0,0,0,0.4);
  }}
  .multi-dropdown.open {{ display: block; }}
  .multi-search {{
    padding: 8px 12px; border-bottom: 1px solid var(--border);
    position: sticky; top: 0; background: var(--surface2); z-index: 1;
  }}
  .multi-search input {{
    width: 100%; background: var(--surface); border: 1px solid var(--border);
    color: var(--text); padding: 6px 10px; border-radius: 4px;
    font-family: var(--font-main); font-size: 0.8rem; outline: none; min-width: unset;
  }}
  .multi-search input:focus {{ border-color: var(--accent); }}
  .multi-actions {{
    display: flex; gap: 8px; padding: 6px 12px;
    border-bottom: 1px solid var(--border);
    position: sticky; top: 45px; background: var(--surface2); z-index: 1;
  }}
  .multi-actions button {{
    background: transparent; border: none; color: var(--accent);
    font-family: var(--font-main); font-size: 0.65rem; cursor: pointer; padding: 2px 4px;
  }}
  .multi-actions button:hover {{ text-decoration: underline; }}
  .multi-option {{
    display: flex; align-items: center; gap: 10px; padding: 9px 14px;
    cursor: pointer; transition: background 0.12s; font-size: 0.85rem;
  }}
  .multi-option:hover {{ background: var(--surface); }}
  .multi-option input[type=checkbox] {{
    accent-color: var(--accent); width: 14px; height: 14px;
    min-width: unset; cursor: pointer;
  }}

  /* TAGS */
  .multi-tag-list {{ display: flex; flex-wrap: wrap; gap: 4px; padding: 4px 40px 12px; }}
  .multi-tag {{
    background: rgba(10,158,116,0.12); color: var(--accent);
    border: 1px solid rgba(10,158,116,0.3); border-radius: 20px;
    padding: 3px 10px; font-size: 0.7rem; font-family: var(--font-main);
    display: flex; align-items: center; gap: 6px;
  }}
  .multi-tag .remove {{ cursor: pointer; opacity: 0.6; font-size: 0.8rem; }}
  .multi-tag .remove:hover {{ opacity: 1; color: var(--accent2); }}
  .reset-btn {{
    background: transparent; border: 1px solid var(--border); color: var(--muted);
    padding: 8px 18px; border-radius: 6px; cursor: pointer;
    font-family: var(--font-main); font-size: 0.7rem; letter-spacing: 1px;
    transition: border-color 0.2s, color 0.2s;
  }}
  .reset-btn:hover {{ border-color: var(--accent2); color: var(--accent2); }}

  /* TABLE */
  .table-wrap {{
    padding: 0 40px 40px; overflow-x: auto; cursor: grab;
    scrollbar-gutter: stable; touch-action: pan-y;
    scrollbar-color: #5f6368 var(--bg); scrollbar-width: thin;
  }}
  .table-wrap.dragging {{ cursor: grabbing; user-select: none; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 0.875rem; }}
  thead tr.header-row {{ background: var(--accent); }}
  thead tr.filter-row {{ background: var(--surface2); border-bottom: 2px solid var(--border); }}
  th {{
    padding: 12px 14px; text-align: left; font-family: var(--font-main);
    font-size: 0.9rem; letter-spacing: 1px; text-transform: uppercase; color: #fff;
    cursor: pointer; user-select: none; white-space: nowrap;
    border-right: 1px solid rgba(255,255,255,0.15);
  }}
  th:last-child {{ border-right: none; }}
  th:hover {{ background: rgba(255,255,255,0.1); }}
  th .sort-icon {{ margin-left: 4px; opacity: 0.4; }}
  th.sorted .sort-icon {{ opacity: 1; }}
  td.filter-cell {{ padding: 6px 10px; }}
  .col-filter {{
    width: 100%; background: var(--surface2); border: 1px solid var(--border);
    color: var(--text); padding: 6px 10px; border-radius: 5px;
    font-family: var(--font-main); font-size: 0.78rem; outline: none;
    transition: border-color 0.2s; min-width: unset;
  }}
  .col-filter:focus {{ border-color: var(--accent); }}
  .col-filter::placeholder {{ color: var(--muted); font-size: 0.72rem; }}
  tbody tr {{ border-bottom: 1px solid var(--border); transition: background 0.15s; }}
  tbody tr:nth-child(even) {{ background: var(--surface2); }}
  tbody tr:hover {{ background: rgba(10,158,116,0.06); }}
  td {{ padding: 12px 16px; vertical-align: middle; }}
  td.sample-cell {{ font-weight: 700; font-size: 1rem; white-space: nowrap; }}
  td.count-cell  {{ text-align: center; font-weight: 600; color: var(--accent); font-size: 1rem; }}
  td.genes-cell  {{ line-height: 2.2; }}

  /* GENE BADGES */
  .gene-badge {{
    display: inline-block; padding: 3px 10px; border-radius: 20px; margin: 2px;
    font-size: 0.9rem; font-family: var(--font-main); font-weight: 600;
    background: rgba(230,126,0,0.1); color: #1a1f2e;
    border: 1px solid rgba(230,126,0,0.3); cursor: default;
    transition: background 0.15s, border-color 0.15s;
  }}
  .gene-badge:hover {{ background: rgba(230,126,0,0.22); border-color: rgba(230,126,0,0.55); }}
  .gene-badge .sups {{ font-size: 0.68em; vertical-align: super; letter-spacing: 0; margin-left: 1px; }}

  /* TOOLTIP */
  #gtooltip {{
    position: fixed; z-index: 9999; display: none;
    background: var(--surface); border: 1px solid var(--border); border-radius: 8px;
    padding: 10px 14px; font-size: 0.78rem; font-family: var(--font-main);
    box-shadow: 0 8px 24px rgba(0,0,0,0.3); max-width: 260px;
    line-height: 1.8; pointer-events: none;
  }}
  #gtooltip .tt-title {{ font-weight: 700; margin-bottom: 4px; font-size: 0.85rem; color: var(--accent); }}
  #gtooltip .tt-meta  {{ font-size: 0.72rem; color: var(--muted); margin-bottom: 6px; }}
  #gtooltip .tt-chip {{
    display: inline-block; padding: 2px 8px; border-radius: 12px; margin: 2px;
    font-size: 0.72rem; font-weight: 600; border: 1px solid;
  }}

  .no-data {{ text-align: center; padding: 60px; color: var(--muted); font-size: 1rem; letter-spacing: 1px; }}

  /* PAGINATION */
  .pagination {{ display: flex; align-items: center; justify-content: flex-end; gap: 8px; padding: 0 40px 40px; flex-wrap: wrap; }}
  .page-info {{ font-size: 0.75rem; color: var(--muted); margin-right: 12px; }}
  .page-btn {{
    background: var(--surface); border: 1px solid var(--border); color: var(--text);
    padding: 6px 14px; border-radius: 6px; cursor: pointer;
    font-family: var(--font-main); font-size: 0.75rem; transition: border-color 0.2s;
  }}
  .page-btn:hover:not(:disabled) {{ border-color: var(--accent); color: var(--accent); }}
  .page-btn:disabled {{ opacity: 0.35; cursor: default; }}
  .page-btn.active {{ background: var(--accent); color: #fff; border-color: var(--accent); font-weight: 700; }}
  ::-webkit-scrollbar {{ width: 6px; height: 6px; }}
  ::-webkit-scrollbar-track {{ background: var(--bg); }}
  ::-webkit-scrollbar-thumb {{ background: #5f6368; border-radius: 3px; }}
  ::-webkit-scrollbar-thumb:hover {{ background: #4b5563; }}
</style>
</head>
<body>

<div id="gtooltip"></div>

<header>
  <div>
    <div class="logo">BASIL <span>//</span> AMR Consensus Report</div>
    <div style="font-size:0.78rem;color:var(--muted);margin-top:4px;font-family:var(--font-main);">
      Minimum Identity ({minid_genes}%) / Minimum Coverage ({mincov_genes}%) / abritAMR 90% Minimum Coverage
    </div>
  </div>
  <div class="export-wrap" tabindex="0">
    <button class="export-btn" id="exportBtn">&#11015; Export &#9662;</button>
    <div class="export-menu">
      <button onclick="exportCSV()"><span>&#128196;</span> CSV</button>
      <button onclick="exportTSV()"><span>&#128203;</span> TSV</button>
      <button onclick="exportExcel()"><span>&#128202;</span> Excel (.xlsx)</button>
      <button onclick="exportPDF()"><span>&#128213;</span> PDF</button>
      <div class="divider"></div>
      <button onclick="copyToClipboard()"><span>&#128203;</span> Copy to Clipboard</button>
    </div>
  </div>
</header>

<div class="stats-section">
  <!-- Summary stats (dynamically updated) -->
  <div class="stats-bar" id="statsBar"></div>
  <!-- Database selection (doubles as the per-database gene-count legend) -->
  <div class="legend-bar" style="justify-content:space-between;">
    <span class="legend-title">
      Databases &nbsp;
      <span style="font-weight:400;text-transform:none;letter-spacing:0;color:var(--muted);">
        (click a card to include/exclude it from the consensus below)
      </span>
    </span>
    <div style="display:flex;gap:6px;">
      <button class="mini-btn" onclick="selectAllDbs()">Select All</button>
      <button class="mini-btn" onclick="selectNoneDbs()">Select None</button>
    </div>
  </div>
  <div class="stats-bar" id="dbStatsBar">
    {db_stat_cards}
  </div>
</div>

<div class="controls">
  <div class="filter-group">
    <div class="filter-label">Samples (multi-select)</div>
    <div class="multi-wrap">
      <div class="multi-trigger" id="sampleTrigger" onclick="toggleMulti()">
        <span id="sampleTriggerText">All Samples</span>
        <span class="arrow">&#9660;</span>
      </div>
      <div class="multi-dropdown" id="sampleDropdown">
        <div class="multi-search">
          <input type="text" id="sampleSearch" placeholder="Search samples..."
                 oninput="filterMultiOptions(this.value)"/>
        </div>
        <div class="multi-actions">
          <button onclick="selectAll()">Select All</button>
          <button onclick="clearAll()">Clear</button>
        </div>
        <div id="sampleOptions">{sample_options}</div>
      </div>
    </div>
  </div>

  <div class="filter-group">
    <div class="filter-label">Rows per page</div>
    <select onchange="setRowsPerPage(Number(this.value))">
      <option value="20" selected>20</option>
      <option value="50">50</option>
      <option value="100">100</option>
    </select>
  </div>

  <button class="reset-btn" onclick="resetFilters()">&#10005; Reset All</button>
</div>

<div class="multi-tag-list" id="sampleTags"></div>

<div class="table-wrap">
  <table id="mainTable">
    <thead>
      <tr class="header-row">
        <th onclick="sortBy('sample_id')"    data-col="sample_id">
          Sample ID <span class="sort-icon">&#8645;</span>
        </th>
        <th onclick="sortBy('total_gene_count')"   data-col="total_gene_count" style="text-align:center;min-width:80px">
          Total AMR genes detected <span class="sort-icon">&#8645;</span>
        </th>
        <th onclick="sortBy('multi_gene_str')" data-col="multi_gene_str">
          Genes detected by &ge;2 databases <span class="sort-icon">&#8645;</span>
        </th>
        <th onclick="sortBy('single_gene_str')" data-col="single_gene_str">
          Genes detected by only one database <span class="sort-icon">&#8645;</span>
        </th>
      </tr>
      <tr class="filter-row">
        <td class="filter-cell">
          <input class="col-filter" id="colSample" placeholder="Filter sample..." oninput="applyFilters()"/>
        </td>
        <td class="filter-cell">
          <input class="col-filter" id="colCount" placeholder="e.g. 5"
                 oninput="applyFilters()" style="text-align:center"/>
        </td>
        <td class="filter-cell">
          <input class="col-filter" id="colMulti" placeholder="Filter consensus genes..." oninput="applyFilters()"/>
        </td>
        <td class="filter-cell">
          <input class="col-filter" id="colSingle" placeholder="Filter single genes..." oninput="applyFilters()"/>
        </td>
      </tr>
    </thead>
    <tbody id="tableBody"></tbody>
  </table>
  <div class="no-data" id="noData" style="display:none">
    No results match the current filters.
  </div>
</div>

<div class="pagination" id="pagination"></div>

<script>
/* Data injected from Python */
const RAW_DATA = {data_json};
const DB_META  = {db_meta_json};

/* State */
let filtered        = [];
let sortCol         = 'sample_id';
let sortAsc         = true;
let currentPage     = 1;
let pageSize        = 20;
let selectedSamples = new Set();
let selectedDbs     = new Set(Object.keys(DB_META).map(Number));

/* Tooltip */
const tooltip = document.getElementById('gtooltip');

function showGeneTooltip(event, gene, dbIds, dbNames) {{
  const chips = dbIds.map((id, i) => {{
    const m = DB_META[id];
    return `<span class="tt-chip"
                  style="background:${{m.color}}22;border-color:${{m.color}}55;color:${{m.color}}">
              ${{m.sup}} ${{dbNames[i]}}
            </span>`;
  }}).join(' ');
  tooltip.innerHTML =
    `<div class="tt-title">${{gene}}</div>` +
    `<div class="tt-meta">Detected by ${{dbIds.length}} of 5 database(s)</div>` +
    `<div>${{chips}}</div>`;
  tooltip.style.display = 'block';
  _positionTip(event);
}}

function _positionTip(e) {{
  const tw = tooltip.offsetWidth, th = tooltip.offsetHeight;
  const x  = e.clientX + 14, y = e.clientY + 14;
  tooltip.style.left = (x + tw > window.innerWidth  ? x - tw - 28 : x) + 'px';
  tooltip.style.top  = (y + th > window.innerHeight ? y - th - 28 : y) + 'px';
}}

document.addEventListener('mousemove', e => {{
  if (tooltip.style.display !== 'none') _positionTip(e);
}});

function hideTooltip() {{ tooltip.style.display = 'none'; }}

/* Multi-select */
function toggleMulti() {{
  const dd   = document.getElementById('sampleDropdown');
  const trig = document.getElementById('sampleTrigger');
  const open = dd.classList.contains('open');
  dd.classList.toggle('open', !open);
  trig.classList.toggle('open', !open);
  if (!open) document.getElementById('sampleSearch').focus();
}}

document.addEventListener('click', e => {{
  if (!e.target.closest('.multi-wrap')) {{
    document.getElementById('sampleDropdown').classList.remove('open');
    document.getElementById('sampleTrigger').classList.remove('open');
  }}
}});

function filterMultiOptions(q) {{
  q = q.toLowerCase();
  document.querySelectorAll('#sampleOptions .multi-option').forEach(opt => {{
    opt.style.display = opt.querySelector('input').value.toLowerCase().includes(q) ? '' : 'none';
  }});
}}

function onMultiChange() {{
  selectedSamples = new Set(
    [...document.querySelectorAll('#sampleOptions input:checked')].map(c => c.value)
  );
  _updateTrigger(); _renderTags(); applyFilters();
}}

function selectAll() {{
  document.querySelectorAll('#sampleOptions .multi-option').forEach(opt => {{
    if (opt.style.display !== 'none') opt.querySelector('input').checked = true;
  }});
  onMultiChange();
}}

function clearAll() {{
  document.querySelectorAll('#sampleOptions input[type=checkbox]').forEach(c => c.checked = false);
  onMultiChange();
}}

function _updateTrigger() {{
  const el = document.getElementById('sampleTriggerText');
  if      (selectedSamples.size === 0) el.textContent = 'All Samples';
  else if (selectedSamples.size === 1) el.textContent = [...selectedSamples][0];
  else                                  el.textContent = `${{selectedSamples.size}} samples selected`;
}}

function _renderTags() {{
  const c = document.getElementById('sampleTags');
  c.innerHTML = '';
  [...selectedSamples].forEach(s => {{
    const div = document.createElement('div');
    div.className = 'multi-tag';
    div.appendChild(document.createTextNode(s + ' '));
    const x = document.createElement('span');
    x.className = 'remove'; x.textContent = '\u2715';
    x.addEventListener('click', (val => () => _removeTag(val))(s));
    div.appendChild(x); c.appendChild(div);
  }});
}}

function _removeTag(s) {{
  selectedSamples.delete(s);
  const cb = document.querySelector(`#sampleOptions input[value='${{s}}']`);
  if (cb) cb.checked = false;
  _updateTrigger(); _renderTags(); applyFilters();
}}

/* Database selection */
function onDbToggle(dbId) {{
  dbId = Number(dbId);
  const cb = document.querySelector(`.db-toggle-checkbox[data-db-id="${{dbId}}"]`);
  if (cb && cb.checked) selectedDbs.add(dbId);
  else selectedDbs.delete(dbId);
  _syncDbCardStyles();
  applyFilters();
}}

function selectAllDbs() {{
  document.querySelectorAll('.db-toggle-checkbox').forEach(cb => cb.checked = true);
  selectedDbs = new Set(Object.keys(DB_META).map(Number));
  _syncDbCardStyles();
  applyFilters();
}}

function selectNoneDbs() {{
  document.querySelectorAll('.db-toggle-checkbox').forEach(cb => cb.checked = false);
  selectedDbs = new Set();
  _syncDbCardStyles();
  applyFilters();
}}

function _syncDbCardStyles() {{
  document.querySelectorAll('.db-toggle').forEach(card => {{
    const id = Number(card.dataset.dbId);
    card.classList.toggle('db-off', !selectedDbs.has(id));
  }});
}}

/* Per-database gene counts, scoped to the current sample filter
   (independent of each database's own on/off state, so a card still
   shows a meaningful number after you switch it off) */
function computeDbCounts() {{
  const colSample = document.getElementById('colSample').value.toLowerCase();
  const counts = {{}};
  Object.keys(DB_META).forEach(k => counts[k] = 0);
  RAW_DATA.forEach(row => {{
    if (selectedSamples.size > 0 && !selectedSamples.has(row.sample_id)) return;
    if (colSample && !row.sample_id.toLowerCase().includes(colSample)) return;
    row.entries.forEach(e => {{
      e.db_ids.forEach(id => {{ counts[id] = (counts[id] || 0) + 1; }});
    }});
  }});
  return counts;
}}

function renderDbCounts() {{
  const counts = computeDbCounts();
  document.querySelectorAll('.db-stat.db-toggle').forEach(card => {{
    const id = card.dataset.dbId;
    const valueEl = card.querySelector('.stat-value');
    if (valueEl) valueEl.textContent = counts[id] ?? 0;
  }});
}}

/* Recompute each sample's gene entries using only the currently selected
   databases -- this is the actual consensus recalculation. A gene keeps
   only the db_ids that are still selected; if none remain it drops out
   entirely; multi/single classification is redone from the surviving count. */
function computeEffectiveRows() {{
  return RAW_DATA.map(row => {{
    const eff = [];
    for (const e of row.entries) {{
      const effIds = e.db_ids.filter(id => selectedDbs.has(id));
      if (effIds.length === 0) continue;
      const sups = effIds.map(id => DB_META[id].sup).join('');
      eff.push({{
        gene: e.gene,
        display: e.gene + sups,
        db_ids: effIds,
        db_names: effIds.map(id => DB_META[id].name),
        db_count: effIds.length,
      }});
    }}
    eff.sort((a, b) => a.gene.toLowerCase().localeCompare(b.gene.toLowerCase()));
    const multi  = eff.filter(e => e.db_count >= 2);
    const single = eff.filter(e => e.db_count === 1);
    return {{
      sample_id:         row.sample_id,
      entries:           eff,
      multi_entries:     multi,
      single_entries:    single,
      total_gene_count:  eff.length,
      multi_gene_count:  multi.length,
      single_gene_count: single.length,
      total_gene_str:    eff.map(e => e.display).join(', '),
      multi_gene_str:    multi.map(e => e.display).join(', '),
      single_gene_str:   single.map(e => e.display).join(', '),
    }};
  }});
}}

/* Filters */
function _geneMatches(entry, query) {{
  if (!query) return true;
  const text = query.toLowerCase();
  return entry.gene.toLowerCase().includes(text) ||
         entry.display.toLowerCase().includes(text);
}}

function applyFilters() {{
  const colSample = document.getElementById('colSample').value.toLowerCase();
  const colCount  = document.getElementById('colCount').value.trim();
  const colMulti  = document.getElementById('colMulti').value.trim();
  const colSingle = document.getElementById('colSingle').value.trim();

  const effRows = computeEffectiveRows();

  filtered = effRows
    .filter(r => {{
      if (selectedSamples.size > 0 && !selectedSamples.has(r.sample_id)) return false;
      if (colSample && !r.sample_id.toLowerCase().includes(colSample))   return false;
      return true;
    }})
    .map(r => {{
      const visMulti  = r.multi_entries.filter(e => _geneMatches(e, colMulti));
      const visSingle = r.single_entries.filter(e => _geneMatches(e, colSingle));
      const vis       = visMulti.concat(visSingle);
      return {{
        ...r,
        visibleEntries: vis,
        visibleMultiEntries: visMulti,
        visibleSingleEntries: visSingle,
        visibleCount: vis.length
      }};
    }})
    .filter(r => r.visibleCount > 0);   // hide samples with no qualifying genes

  // Numeric gene-count filter
  if (colCount !== '') {{
    const n = parseInt(colCount);
    if (!isNaN(n)) filtered = filtered.filter(r => r.visibleCount === n);
  }}

  sortData();
  currentPage = 1;
  renderStats(filtered);
  renderDbCounts();
  renderTable(filtered, currentPage);
}}

function resetFilters() {{
  clearAll();
  selectAllDbs();
  ['colSample','colCount','colMulti','colSingle'].forEach(id => document.getElementById(id).value = '');
  document.getElementById('sampleSearch').value = '';
  filterMultiOptions('');
  applyFilters();
}}

/* Stats bar (dynamic) */
function renderStats(data) {{
  const samples   = new Set(data.map(r => r.sample_id)).size;
  const geneSet   = new Set(data.flatMap(r => r.visibleEntries.map(e => e.gene.toLowerCase())));
  const multiConf = data.flatMap(r => r.visibleEntries).filter(e => e.db_count >= 2).length;
  const avgGenes  = samples > 0
    ? (data.reduce((a, r) => a + r.visibleCount, 0) / samples).toFixed(1)
    : '0';

  document.getElementById('statsBar').innerHTML = `
    <div class="stat-card">
      <div class="stat-label">Samples</div>
      <div class="stat-value">${{samples}}</div>
    </div>
    <div class="stat-card">
      <div class="stat-label">Unique Genes</div>
      <div class="stat-value">${{geneSet.size}}</div>
    </div>
    <div class="stat-card">
      <div class="stat-label">Avg Genes / Sample</div>
      <div class="stat-value">${{avgGenes}}</div>
    </div>
    <div class="stat-card">
      <div class="stat-label">&ge;2 DB Confirmed</div>
      <div class="stat-value">${{multiConf}}</div>
    </div>
  `;
}}

/* Table */
function renderGeneBadges(entries) {{
  if (!entries.length) return '<span style="color:var(--muted)">&#8212;</span>';
  return entries.map(e => {{
    const dbIdsJ   = JSON.stringify(e.db_ids);
    const dbNamesJ = JSON.stringify(e.db_names);
    const sups     = e.db_ids.map(i => DB_META[i].sup).join('');
    return `<span class="gene-badge"
                  onmouseenter="showGeneTooltip(event,${{JSON.stringify(e.gene)}},${{dbIdsJ}},${{dbNamesJ}})"
                  onmouseleave="hideTooltip()">
              ${{e.gene}}<span class="sups">${{sups}}</span>
            </span>`;
  }}).join('');
}}

function renderTable(data, page) {{
  const pages    = Math.max(1, Math.ceil(data.length / pageSize));
  const safe     = Math.min(Math.max(page, 1), pages);
  currentPage    = safe;
  const slice    = data.slice((safe - 1) * pageSize, safe * pageSize);
  const tbody    = document.getElementById('tableBody');
  const noData   = document.getElementById('noData');

  if (!data.length) {{
    tbody.innerHTML = '';
    noData.textContent = selectedDbs.size === 0
      ? 'No databases selected \u2014 enable at least one database above to see results.'
      : 'No results match the current filters.';
    noData.style.display = 'block';
    document.getElementById('pagination').innerHTML = ''; return;
  }}
  noData.style.display = 'none';

  tbody.innerHTML = slice.map(r => {{
    return `<tr>
      <td class="sample-cell">${{r.sample_id}}</td>
      <td class="count-cell">${{r.visibleCount}}</td>
      <td class="genes-cell">${{renderGeneBadges(r.visibleMultiEntries)}}</td>
      <td class="genes-cell">${{renderGeneBadges(r.visibleSingleEntries)}}</td>
    </tr>`;
  }}).join('');

  renderPagination(data.length, safe);
}}

function renderPagination(total, page) {{
  const pages = Math.ceil(total / pageSize);
  const el    = document.getElementById('pagination');
  if (pages <= 1) {{ el.innerHTML = ''; return; }}

  let btns = `<span class="page-info">Showing ${{Math.min((page-1)*pageSize+1,total)}}&ndash;${{Math.min(page*pageSize,total)}} of ${{total}}</span>`;
  btns += `<button class="page-btn" onclick="goPage(${{page-1}})" ${{page===1?'disabled':''}}>&#8249; Prev</button>`;
  const range = [];
  for (let i = Math.max(1,page-2); i<=Math.min(pages,page+2); i++) range.push(i);
  if (range[0]>1) btns+=`<button class="page-btn" onclick="goPage(1)">1</button>${{range[0]>2?'<span style="color:var(--muted);padding:0 4px">&hellip;</span>':''}}`;
  range.forEach(i => btns+=`<button class="page-btn ${{i===page?'active':''}}" onclick="goPage(${{i}})">${{i}}</button>`);
  if (range.at(-1)<pages) btns+=`${{range.at(-1)<pages-1?'<span style="color:var(--muted);padding:0 4px">&hellip;</span>':''}}<button class="page-btn" onclick="goPage(${{pages}})">${{pages}}</button>`;
  btns+=`<button class="page-btn" onclick="goPage(${{page+1}})" ${{page===pages?'disabled':''}}>Next &#8250;</button>`;
  el.innerHTML = btns;
}}

function goPage(p) {{
  currentPage = p;
  renderTable(filtered, currentPage);
  window.scrollTo({{top: 0, behavior: 'smooth'}});
}}

function setRowsPerPage(v) {{ pageSize = v; currentPage = 1; renderTable(filtered, currentPage); }}

/* Sort */
function sortData() {{
  filtered.sort((a, b) => {{
    if (sortCol === 'total_gene_count') {{
      return sortAsc ? a.visibleCount - b.visibleCount : b.visibleCount - a.visibleCount;
    }}
    const va = a[sortCol] || '';
    const vb = b[sortCol] || '';
    return sortAsc ? va.localeCompare(vb) : vb.localeCompare(va);
  }});
}}

function sortBy(col) {{
  if (sortCol === col) sortAsc = !sortAsc;
  else {{ sortCol = col; sortAsc = true; }}
  document.querySelectorAll('th[data-col]').forEach(th => {{
    th.classList.toggle('sorted', th.dataset.col === col);
    th.querySelector('.sort-icon').innerHTML =
      th.dataset.col === col ? (sortAsc ? '&#8593;' : '&#8595;') : '&#8645;';
  }});
  sortData();
  renderTable(filtered, currentPage);
}}

/* Drag-scroll */
function initDragScroll() {{
  document.querySelectorAll('.table-wrap').forEach(wrap => {{
    let drag=false, startX=0, startL=0, moved=false;
    wrap.addEventListener('pointerdown', e => {{
      if (e.button!==0 || wrap.scrollWidth<=wrap.clientWidth) return;
      if (!e.target.closest('tbody td')) return;
      drag=true; moved=false; startX=e.clientX; startL=wrap.scrollLeft;
      wrap.classList.add('dragging');
      if (wrap.setPointerCapture) wrap.setPointerCapture(e.pointerId);
    }});
    wrap.addEventListener('pointermove', e => {{
      if (!drag) return;
      const d = e.clientX - startX;
      if (Math.abs(d)>4) moved=true;
      wrap.scrollLeft = startL - d;
    }});
    const stop = e => {{
      if (!drag) return; drag=false; wrap.classList.remove('dragging');
      try {{ if (wrap.hasPointerCapture?.(e?.pointerId)) wrap.releasePointerCapture(e.pointerId); }} catch(_) {{}}
    }};
    wrap.addEventListener('pointerup',     stop);
    wrap.addEventListener('pointercancel', stop);
    wrap.addEventListener('click', e => {{ if (!moved) return; e.preventDefault(); e.stopPropagation(); moved=false; }}, true);
  }});
}}

/* Export helpers */
function _exportRows() {{
  return filtered.map(r => [
    r.sample_id,
    r.visibleCount,
    r.visibleMultiEntries.map(e => e.display).join(', '),
    r.visibleSingleEntries.map(e => e.display).join(', '),
  ]);
}}

function exportCSV() {{
  const h = ['Sample_ID','Total_AMR_Genes','Genes_2plus_DB','Genes_1_DB'];
  const body = _exportRows().map(r => r.map(v => `"${{String(v).replace(/"/g,'""')}}"`).join(','));
  // Prepend a UTF-8 BOM so Excel correctly detects superscript characters.
  const csv = '\ufeff' + [h.map(v=>`"${{v}}"`).join(','), ...body].join('\\r\\n');
  downloadBlob(csv, 'amr_consensus_report.csv', 'text/csv;charset=utf-8');
}}

function exportTSV() {{
  const h = ['Sample_ID','Total_AMR_Genes','Genes_2plus_DB','Genes_1_DB'];
  const body = _exportRows().map(r => r.join('\\t'));
  const tsv = '\ufeff' + [h.join('\\t'), ...body].join('\\r\\n');
  downloadBlob(tsv, 'amr_consensus_report.tsv', 'text/tab-separated-values;charset=utf-8');
}}

function exportExcel() {{
  const header  = ['Sample ID', 'Total AMR Genes', 'Genes >=2 DBs', 'Genes 1 DB'];
  const wb = XLSX.utils.book_new();
  const ws = XLSX.utils.aoa_to_sheet([header, ..._exportRows()]);
  ws['!cols'] = [{{wch:30}}, {{wch:14}}, {{wch:70}}, {{wch:70}}];
  XLSX.utils.book_append_sheet(wb, ws, 'Consensus AMR');
  XLSX.writeFile(wb, 'amr_consensus_report.xlsx');
}}

function exportPDF() {{
  const {{ jsPDF }} = window.jspdf;
  const doc = new jsPDF({{ orientation:'landscape', unit:'mm', format:'a4' }});
  doc.setFont('helvetica','bold'); doc.setFontSize(14); doc.setTextColor(10,158,116);
  doc.text('BASIL \u2014 AMR Consensus Report', 14, 16);
  doc.setFont('helvetica','normal'); doc.setFontSize(8); doc.setTextColor(108,117,125);
  doc.text(`Generated: ${{new Date().toLocaleString()}} \u00b7 Samples: ${{filtered.length}}`, 14, 22);
  doc.autoTable({{
    startY: 27,
    head: [['Sample ID','Total AMR Genes','Genes >=2 DBs','Genes 1 DB']],
    body: _exportRows(),
    styles: {{ font:'helvetica', fontSize:7, cellPadding:3, overflow:'linebreak' }},
    headStyles: {{ fillColor:[10,158,116], textColor:[255,255,255], fontStyle:'bold', fontSize:8 }},
    alternateRowStyles: {{ fillColor:[245,247,250] }},
    columnStyles: {{ 0:{{cellWidth:45}}, 1:{{cellWidth:22}}, 2:{{cellWidth:'auto'}}, 3:{{cellWidth:'auto'}} }},
    theme: 'grid',
  }});
  doc.save('amr_consensus_report.pdf');
}}

async function copyToClipboard() {{
  const h    = 'Sample_ID\\tTotal_AMR_Genes\\tGenes_2plus_DB\\tGenes_1_DB';
  const body = _exportRows().map(r => r.join('\\t'));
  const text = [h, ...body].join('\\n');
  const btn  = document.getElementById('exportBtn');
  const orig = btn.textContent;
  try {{
    await navigator.clipboard.writeText(text);
  }} catch (_) {{
    const ta = document.createElement('textarea');
    ta.value = text; document.body.appendChild(ta);
    ta.select(); document.execCommand('copy'); document.body.removeChild(ta);
  }}
  btn.textContent = '\u2713 Copied!';
  setTimeout(() => btn.textContent = orig, 1600);
}}

function downloadBlob(content, filename, mime) {{
  const a   = document.createElement('a');
  a.href    = URL.createObjectURL(new Blob([content], {{type: mime}}));
  a.download = filename; a.click();
}}

/* Init */
initDragScroll();
applyFilters();
</script>
</body>
</html>"""

    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write(html)
    print(f"  HTML report : {output_path}")

# Main

def main():
    parser = argparse.ArgumentParser(
        description="BASIL — AMR Consensus Integration Report generator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--abricate_dir", required=True,
        help="Directory containing Abricate TSV files "
             "(ncbi.tsv, card.tsv, resfinder.tsv, argannot.tsv)",
    )
    parser.add_argument(
        "--abritamr", required=True,
        help="AbritAMR summary_matches.tsv file (AMRFinderPlus output)",
    )
    parser.add_argument(
        "--output", required=True,
        help="Output HTML report path",
    )
    parser.add_argument(
        "--output_tsv", default=None,
        help="Optional: long-format TSV summary output path",
    )
    parser.add_argument(
        "--minid_genes", type=float, default=90.0,
        help="Minimum identity threshold for gene calls (default: 90.0)"
    )
    parser.add_argument(
        "--mincov_genes", type=float, default=90.0,
        help="Minimum coverage threshold for gene calls (default: 90.0)"
    )
    args = parser.parse_args()

    print("\n\u2500\u2500 AMR Consensus Integration \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")
    print(f"  Abricate dir : {args.abricate_dir}")
    print(f"  AbritAMR     : {args.abritamr}")
    print(f"  HTML output  : {args.output}")
    if args.output_tsv:
        print(f"  TSV output   : {args.output_tsv}")
    print()

    # Parse all databases
    all_db_genes = parse_abricate_dir(args.abricate_dir)
    all_db_genes[5] = parse_abritamr_tsv(args.abritamr)

    # Build consensus
    print("\n  Building non-redundant consensus ...")
    consensus = build_consensus(all_db_genes)
    total_calls = sum(len(v) for v in consensus.values())
    print(f"  \u2192 {len(consensus)} samples | {total_calls} consensus gene calls")

    # Write outputs
    print()
    if args.output_tsv:
        write_tsv(consensus, args.output_tsv)
    generate_html(consensus, args.output, minid_genes=args.minid_genes, mincov_genes=args.mincov_genes)

    print("\u2500\u2500 Done \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")


if __name__ == "__main__":
    main()