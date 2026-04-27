#!/usr/bin/env python3

import os
import sys
import glob
import json
import re


def parse_assembly_summary(filepath):
    d = {
        "mean_read_length": None, "coverage": None,
        "kmers": None, "assembly_mode": None,
        "asm1_contigs": None, "asm1_length": None, "asm1_n50": None,
        "asm2_contigs": None, "asm2_length": None, "asm2_n50": None,
        "delta_contigs": None, "delta_length": None, "delta_n50": None,
        "best_assembly": None, "best_reasons": [],
        "filt_min_length": None,
    }
    if not filepath or not os.path.exists(filepath):
        return d

    section = None
    with open(filepath) as f:
        for line in f:
            s = line.strip()
            if "ASSEMBLY CONFIGURATION" in s:    section = "config"
            elif "ASSEMBLY 1:" in s:             section = "asm1"
            elif "ASSEMBLY 2:" in s:             section = "asm2"
            elif "COMPARISON" in s:              section = "cmp"
            elif "BEST ASSEMBLY:" in s:
                section = "best"
                m = re.search(r"Assembly (\d+)", s)
                if m: d["best_assembly"] = int(m.group(1))
            elif s.startswith("==="):
                pass
            else:
                if section == "config":
                    if s.startswith("MEAN_READ_LENGTH="):
                        d["mean_read_length"] = float(s.split("=",1)[1])
                    elif s.startswith("COVERAGE="):
                        d["coverage"] = float(s.split("=",1)[1])
                    elif s.startswith("KMERS_SPADES="):
                        d["kmers"] = s.split("=",1)[1]
                    elif s.startswith("ASSEMBLY_MODE="):
                        d["assembly_mode"] = s.split("=",1)[1]
                elif section == "asm1":
                    if s.startswith("Number of contigs:"):
                        d["asm1_contigs"] = int(s.split(":",1)[1].strip().replace(",",""))
                    elif s.startswith("Total assembly length:"):
                        d["asm1_length"] = int(s.split(":",1)[1].strip().replace(",","").replace(" bp",""))
                    elif s.startswith("N50:"):
                        d["asm1_n50"] = int(s.split(":",1)[1].strip().replace(",","").replace(" bp",""))
                elif section == "asm2":
                    if s.startswith("Number of contigs:"):
                        d["asm2_contigs"] = int(s.split(":",1)[1].strip().replace(",",""))
                    elif s.startswith("Total assembly length:"):
                        d["asm2_length"] = int(s.split(":",1)[1].strip().replace(",","").replace(" bp",""))
                    elif s.startswith("N50:"):
                        d["asm2_n50"] = int(s.split(":",1)[1].strip().replace(",","").replace(" bp",""))
                elif section == "cmp":
                    if s.startswith("\u0394 Contigs:") or s.startswith("Contigs:"):
                        d["delta_contigs"] = s.split(":",1)[1].strip()
                    elif s.startswith("\u0394 Length:") or s.startswith("Length:"):
                        d["delta_length"] = s.split(":",1)[1].strip()
                    elif s.startswith("\u0394 N50:") or s.startswith("N50:"):
                        d["delta_n50"] = s.split(":",1)[1].strip()
                elif section == "best":
                    if s.startswith("\u2705") or s.startswith("-"):
                        d["best_reasons"].append(s)
                    elif s == "Similar quality metrics":
                        d["best_reasons"].append(s)
                    elif "min_length=" in s:
                        m = re.search(r"min_length=(\d+)", s)
                        if m: d["filt_min_length"] = int(m.group(1))
    return d


def parse_filtering_report(filepath):
    d = {
        "filt_min_length": None, "filt_min_cov": None,
        "filt_input_contigs": None, "filt_output_contigs": None,
        "filt_input_bp": None, "filt_output_bp": None,
        "filt_retention": None,
    }
    if not filepath or not os.path.exists(filepath):
        return d

    with open(filepath) as f:
        for line in f:
            s = line.strip()
            if s.startswith("Minimum length:"):
                d["filt_min_length"] = int(s.split(":",1)[1].strip().replace(" bp",""))
            elif s.startswith("Minimum coverage:"):
                d["filt_min_cov"] = float(s.split(":",1)[1].strip().replace(" X",""))
            elif s.startswith("Total input contigs:"):
                d["filt_input_contigs"] = int(s.split(":",1)[1].strip())
            elif s.startswith("Total output contigs:"):
                d["filt_output_contigs"] = int(s.split(":",1)[1].strip())
            elif s.startswith("Input:"):
                d["filt_input_bp"] = int(s.split(":",1)[1].strip().replace(" bp",""))
            elif s.startswith("Output:"):
                d["filt_output_bp"] = int(s.split(":",1)[1].strip().replace(" bp",""))
            elif s.startswith("Retention:"):
                d["filt_retention"] = float(s.split(":",1)[1].strip().replace(" %",""))
    return d


def parse_polishing_report(filepath):
    d = {
        # Basic stats
        "pol_output_contigs":           None,
        "pol_output_length":            None,
        "pol_mean_depth":               None,
        "pol_success":                  False,
        "pol_total_reads":              None,
        "pol_total_polished_contigs":   None,
        "pol_total_contigs":            None,
        "pol_total_polished":           None,
        "pol_total_sub":                None,
        "pol_total_ins":                None,
        "pol_total_del":                None,
        # Gap filling
        "pol_contigs_gap_filled":       None,
        "pol_total_gaps_filled":        None,
        "pol_total_gap_bases_original": None,
        "pol_total_bases_filled":       None,
        "pol_total_bases_filled_pair":  None,
        # Local misassembly
        "pol_contigs_rearranged":       None,
        "pol_total_local_fixes":        None,
    }
    if not filepath or not os.path.exists(filepath):
        return d

    section = None
    with open(filepath) as f:
        for line in f:
            s = line.strip()

            if s.startswith("Assembly stats:"):                       section = "asm"
            elif s.startswith("Input Assembly:"):                     section = "input_asm"
            elif s.startswith("Polished Assembly:"):                  section = "polished_asm"
            elif s.startswith("Mapping stats:") or s.startswith("Mapping Statistics:"):
                section = "mapping"
            elif s.startswith("Coverage:") or s.startswith("Coverage Statistics:"):
                section = "cov"
            elif s.startswith("Correction summary:"):                 section = "correction"
            elif s.startswith("Gap filling summary:"):                section = "gap"
            elif s.startswith("Local misassembly rearrangements:"):   section = "local"
            elif "Pilon polishing completed" in s or "Polishing completed successfully" in s:
                d["pol_success"] = True

            if section in {"asm", "polished_asm"}:
                if s.startswith("Contigs:"):
                    d["pol_output_contigs"] = int(s.split(":",1)[1].strip())
                elif s.startswith("Total length:"):
                    d["pol_output_length"] = int(s.split(":",1)[1].strip().replace(" bp",""))

            elif section == "mapping":
                m = re.match(r"^(\d+)\s+\+\s+\d+\s+in total", s)
                if m: d["pol_total_reads"] = int(m.group(1))

            elif section == "cov":
                if s.startswith("Mean depth:"):
                    d["pol_mean_depth"] = float(s.split(":",1)[1].strip().replace(" X",""))

            elif section == "correction":
                if s.startswith("Polished contigs:"):
                    parts = s.split(":",1)[1].strip().split("/")
                    d["pol_total_polished_contigs"] = int(parts[0].strip())
                    d["pol_total_contigs"] = int(parts[1].strip()) if len(parts) > 1 else None
                elif s.startswith("Total bases changed:"):
                    d["pol_total_polished"] = int(s.split(":",1)[1].strip())
                elif s.startswith("Substituted bases:"):
                    d["pol_total_sub"] = int(s.split(":",1)[1].strip())
                elif s.startswith("Inserted bases:"):
                    d["pol_total_ins"] = int(s.split(":",1)[1].strip())
                elif s.startswith("Deleted bases:"):
                    d["pol_total_del"] = int(s.split(":",1)[1].strip())

            elif section == "gap":
                if s.startswith("Contigs with gaps filled:"):
                    d["pol_contigs_gap_filled"] = int(s.split(":",1)[1].strip())
                elif s.startswith("Total gaps filled:"):
                    d["pol_total_gaps_filled"] = int(s.split(":",1)[1].strip())
                elif s.startswith("Total bases filled:"):
                    raw_val = s.split(":",1)[1].strip()
                    if "/" in raw_val:
                        left, right = [x.strip() for x in raw_val.split("/", 1)]
                        if left:
                            try:
                                d["pol_total_gap_bases_original"] = int(left)
                            except ValueError:
                                pass
                        if right:
                            try:
                                d["pol_total_bases_filled"] = int(right)
                            except ValueError:
                                pass
                        if (
                            d["pol_total_gap_bases_original"] is not None
                            and d["pol_total_bases_filled"] is not None
                        ):
                            d["pol_total_bases_filled_pair"] = (
                                f"{d['pol_total_gap_bases_original']}/{d['pol_total_bases_filled']}"
                            )
                    else:
                        d["pol_total_bases_filled"] = int(raw_val)

            elif section == "local":
                if s.startswith("Contigs rearranged:"):
                    d["pol_contigs_rearranged"] = int(s.split(":",1)[1].strip())
                elif s.startswith("Total local fixes:"):
                    d["pol_total_local_fixes"] = int(s.split(":",1)[1].strip())

    return d


def collect_samples(input_dir):
    summaries = sorted(glob.glob(os.path.join(input_dir, "*_assembly_summary.txt")))
    records = []
    for sf in summaries:
        sample_id = os.path.basename(sf).replace("_assembly_summary.txt", "")
        ff = os.path.join(input_dir, f"{sample_id}_filtering_report.txt")
        pf = os.path.join(input_dir, f"{sample_id}_polishing_report.txt")

        rec = {"sample_id": sample_id}
        rec.update(parse_assembly_summary(sf))
        rec.update(parse_filtering_report(ff if os.path.exists(ff) else None))

        pol_path = pf if (os.path.exists(pf) and os.path.basename(pf) != "NO_FILE") else None
        rec.update(parse_polishing_report(pol_path))

        records.append(rec)
        print(f"  Parsed: {sample_id}")
    return records


def main(input_dir, output_html):
    records = collect_samples(input_dir)
    if not records:
        print(f"No *_assembly_summary.txt files found in: {input_dir}")
        sys.exit(1)
    print(f"  Total samples: {len(records)}")
    data_json = json.dumps(records)

    css = """
  :root{--bg:#f8f9fa;--surface:#fff;--surface2:#f0f2f5;--border:#dee2e6;
    --accent:#0a9e74;--accent2:#e53935;--accent3:#e67e00;--accent4:#1976d2;
    --text:#1a1f2e;--muted:#6c757d;
    --font-main:Calibri,'Calibri',Arial,sans-serif}
  *,*::before,*::after{box-sizing:border-box;margin:0;padding:0}
  body{background:var(--bg);color:var(--text);font-family:var(--font-main);min-height:100vh}
  header{background:var(--surface);border-bottom:2px solid var(--accent);padding:20px 40px;
    display:flex;align-items:center;justify-content:space-between;gap:16px;flex-wrap:wrap;
    box-shadow:0 2px 8px rgba(0,0,0,.06)}
  .logo{font-family:var(--font-main);font-size:1.2rem;font-weight:700;color:var(--accent);letter-spacing:2px;text-transform:uppercase}
  .logo span{color:var(--muted)}
  .export-wrap{position:relative}
  .export-btn{background:var(--accent);color:#fff;border:none;padding:10px 22px;border-radius:6px;
    font-family:var(--font-main);font-size:.85rem;font-weight:700;letter-spacing:1px;cursor:pointer;transition:opacity .2s}
  .export-btn:hover{opacity:.85}
  .export-menu{display:none;position:absolute;top:calc(100% + 6px);right:0;background:var(--surface);
    border:1px solid var(--border);border-radius:8px;overflow:hidden;min-width:160px;z-index:100;
    box-shadow:0 8px 24px rgba(0,0,0,.12)}
  .export-wrap:hover .export-menu,.export-wrap:focus-within .export-menu{display:block}
  .export-menu button{display:flex;align-items:center;gap:10px;width:100%;background:transparent;
    border:none;color:var(--text);padding:11px 18px;font-family:var(--font-main);font-size:.85rem;
    cursor:pointer;transition:background .15s;text-align:left}
  .export-menu button:hover{background:var(--surface2);color:var(--accent)}
  .stats-bar{display:flex;gap:16px;padding:24px 40px;flex-wrap:wrap}
  .stat-card{background:var(--surface);border:1px solid var(--border);border-radius:10px;
    padding:16px 24px;flex:1;min-width:140px;box-shadow:0 1px 4px rgba(0,0,0,.05)}
  .stat-label{font-size:.7rem;font-family:var(--font-main);color:var(--muted);text-transform:uppercase;
    letter-spacing:1.5px;margin-bottom:6px;font-weight:600}
  .stat-value{font-size:1.8rem;font-weight:700;color:var(--accent);font-family:var(--font-main)}
  .tab-bar{display:flex;gap:0;padding:0 40px;border-bottom:2px solid var(--border)}
  .tab-btn{background:transparent;border:none;border-bottom:3px solid transparent;
    padding:12px 22px;font-family:var(--font-main);font-size:.85rem;font-weight:600;letter-spacing:.5px;
    text-transform:uppercase;cursor:pointer;color:var(--muted);margin-bottom:-2px;transition:all .15s}
  .tab-btn:hover{color:var(--text)}
  .tab-btn.active{color:var(--accent);border-bottom-color:var(--accent);font-weight:700}
  .tab-panel{display:none}
  .tab-panel.active{display:block}
  .controls{padding:16px 40px 8px;display:flex;gap:12px;flex-wrap:wrap;align-items:center;justify-content:space-between}
  .control-left,.control-right{display:flex;gap:12px;flex-wrap:wrap;align-items:center}
  .reset-btn{background:transparent;border:1px solid var(--border);color:var(--muted);padding:8px 18px;
    border-radius:6px;cursor:pointer;font-family:var(--font-main);font-size:.8rem;letter-spacing:.5px;
    font-weight:600;transition:border-color .2s,color .2s}
  .reset-btn:hover{border-color:var(--accent2);color:var(--accent2)}
  .rows-wrap{display:flex;align-items:center;gap:8px;font-family:var(--font-main);font-size:.8rem;color:var(--muted)}
  .rows-wrap select{background:var(--surface);border:1px solid var(--border);color:var(--text);padding:7px 10px;border-radius:6px;font-family:var(--font-main);font-size:.8rem;cursor:pointer}
  .rows-wrap select:focus{outline:none;border-color:var(--accent)}
  .table-wrap{padding:0 40px 16px;overflow-x:auto;cursor:grab;scrollbar-gutter:stable;touch-action:pan-y;scrollbar-color:#5f6368 var(--bg);scrollbar-width:thin}
  .table-wrap.dragging{cursor:grabbing;user-select:none}
  .pager{padding:0 40px 32px;display:flex;justify-content:flex-end;align-items:center;gap:10px;flex-wrap:wrap}
  .pager-info{font-family:var(--font-main);font-size:.8rem;color:var(--muted)}
  .pager button{background:var(--surface);border:1px solid var(--border);color:var(--text);padding:7px 12px;border-radius:6px;font-family:var(--font-main);font-size:.8rem;cursor:pointer;transition:border-color .2s,color .2s}
  .pager button:hover:not(:disabled){border-color:var(--accent);color:var(--accent)}
  .pager button:disabled{opacity:.45;cursor:not-allowed}
  table{width:100%;border-collapse:collapse;font-size:.85rem}
  thead tr.header-row{background:var(--accent)}
  thead tr.filter-row{background:var(--surface2);border-bottom:2px solid var(--border)}
  th{padding:11px 13px;text-align:left;font-family:var(--font-main);font-size:.75rem;font-weight:700;
    letter-spacing:.8px;text-transform:uppercase;color:#fff;cursor:pointer;user-select:none;
    white-space:nowrap;border-right:1px solid rgba(255,255,255,.15)}
  th:last-child{border-right:none}
  th:hover{background:rgba(255,255,255,.1)}
  th .si{margin-left:4px;opacity:.6}
  th.sorted .si{opacity:1}
  td.fc{padding:4px 5px;vertical-align:top}
  .ms-wrap{position:relative;min-width:90px}
  .ms-trigger{background:var(--surface);border:1px solid var(--border);color:var(--text);
    padding:4px 8px;border-radius:4px;font-family:var(--font-main);font-size:.75rem;
    cursor:pointer;display:flex;justify-content:space-between;align-items:center;
    gap:4px;user-select:none;transition:border-color .2s;white-space:nowrap;min-width:90px}
  .ms-trigger:hover,.ms-trigger.open{border-color:var(--accent)}
  .ms-trigger .arr{font-size:.55rem;color:var(--muted);transition:transform .2s;flex-shrink:0}
  .ms-trigger.open .arr{transform:rotate(180deg)}
  .ms-dropdown{display:none;position:absolute;top:calc(100% + 2px);left:0;
    background:var(--surface);border:1px solid var(--border);border-radius:6px;
    z-index:500;min-width:160px;max-height:220px;overflow-y:auto;
    box-shadow:0 6px 20px rgba(0,0,0,.12)}
  .ms-dropdown.open{display:block}
  .ms-search-row{padding:6px 8px;border-bottom:1px solid var(--border);
    position:sticky;top:0;background:var(--surface);z-index:1}
  .ms-search-row input{width:100%;background:var(--surface2);border:1px solid var(--border);
    color:var(--text);padding:4px 8px;border-radius:4px;font-size:.75rem;
    outline:none;min-width:unset}
  .ms-search-row input:focus{border-color:var(--accent)}
  .ms-actions{display:flex;gap:8px;padding:4px 8px;border-bottom:1px solid var(--border);
    position:sticky;top:37px;background:var(--surface);z-index:1}
  .ms-actions button{background:transparent;border:none;color:var(--accent);
    font-family:var(--font-main);font-size:.7rem;font-weight:600;cursor:pointer;padding:1px 3px}
  .ms-actions button:hover{text-decoration:underline}
  .ms-opt{display:flex;align-items:center;gap:8px;padding:7px 10px;
    cursor:pointer;transition:background .12s;font-size:.8rem;font-family:var(--font-main)}
  .ms-opt:hover{background:var(--surface2)}
  .ms-opt input[type=checkbox]{accent-color:var(--accent);width:13px;height:13px;
    min-width:unset;cursor:pointer;flex-shrink:0}
  .ms-badge-list{display:flex;flex-wrap:wrap;gap:2px;padding:3px 40px 4px;min-height:0}
  .ms-badge{background:rgba(10,158,116,.1);color:var(--accent);border:1px solid rgba(10,158,116,.25);
    border-radius:12px;padding:2px 8px;font-size:.7rem;font-family:var(--font-main);
    display:flex;align-items:center;gap:4px}
  .ms-badge .x{cursor:pointer;opacity:.6;font-size:.65rem}
  .ms-badge .x:hover{opacity:1;color:var(--accent2)}
  tbody tr{border-bottom:1px solid var(--border);transition:background .12s}
  tbody tr:nth-child(even){background:var(--surface2)}
  tbody tr:hover{background:rgba(10,158,116,.06)}
  td{padding:9px 13px;vertical-align:middle;border-right:1px solid var(--border);white-space:nowrap}
  td:last-child{border-right:none}
  td.sc{font-family:var(--font-main);font-size:.85rem;font-weight:700;color:var(--text)}
  .badge{display:inline-block;padding:2px 9px;border-radius:12px;font-size:.75rem;
    font-family:var(--font-main);font-weight:600}
  .b-green{background:rgba(10,158,116,.12);color:#000;border:1px solid rgba(10,158,116,.3)}
  .b-blue{background:rgba(25,118,210,.1);color:#000;border:1px solid rgba(25,118,210,.25)}
  .b-orange{background:rgba(230,126,0,.1);color:#000;border:1px solid rgba(230,126,0,.25)}
  .b-grey{background:rgba(222,226,230,.6);color:#000;border:1px solid var(--border)}
  .b-yes{background:rgba(10,158,116,.12);color:#000;border:1px solid rgba(10,158,116,.3)}
  .b-no{background:rgba(222,226,230,.6);color:#000;border:1px solid var(--border)}
  .b-sub{background:rgba(25,118,210,.1);color:#000;border:1px solid rgba(25,118,210,.25)}
  .b-ins{background:rgba(10,158,116,.1);color:#000;border:1px solid rgba(10,158,116,.25)}
  .b-del{background:rgba(229,57,53,.1);color:#000;border:1px solid rgba(229,57,53,.25)}
  .bar-wrap{display:flex;align-items:center;gap:6px;min-width:110px}
  .bar-track{flex:1;height:5px;background:var(--border);border-radius:3px;min-width:50px}
  .bar-fill{height:100%;border-radius:3px}
  .bar-num{font-family:var(--font-main);font-size:.78rem;color:var(--text);min-width:48px;text-align:right}
  .no-data{text-align:center;padding:60px;color:var(--muted);font-family:var(--font-main);font-size:.9rem}
  ::-webkit-scrollbar{width:6px;height:6px}
  ::-webkit-scrollbar-track{background:var(--bg)}
  ::-webkit-scrollbar-thumb{background:#5f6368;border-radius:3px}
  ::-webkit-scrollbar-thumb:hover{background:#4b5563}
"""

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1.0"/>
<title>Assembly Report</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/xlsx/0.18.5/xlsx.full.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/jspdf/2.5.1/jspdf.umd.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/jspdf-autotable/3.8.2/jspdf.plugin.autotable.min.js"></script>
<style>{css}</style>
</head>
<body>

<header>
  <div>
    <div class="logo">BASIL <span>//</span> Assembly Report</div>
    <div style="font-size:.85rem;color:var(--muted);margin-top:4px;font-family:var(--font-main)">
      Assembly &middot; Polishing &middot; Filtering Summary
    </div>
  </div>
  <div class="export-wrap" tabindex="0">
    <button class="export-btn">&#11015; Export &#9662;</button>
    <div class="export-menu">
      <button onclick="exportCSV()">&#128196; CSV</button>
      <button onclick="exportTSV()">&#128203; TSV</button>
      <button onclick="exportExcel()">&#128202; Excel (.xlsx)</button>
      <button onclick="exportPDF()">&#128213; PDF</button>
    </div>
  </div>
</header>

<div class="stats-bar" id="statsBar"></div>

<div class="tab-bar">
  <button class="tab-btn active" onclick="switchTab('assembly')">Assembly</button>
  <button class="tab-btn"        onclick="switchTab('polishing')">Polishing</button>
  <button class="tab-btn"        onclick="switchTab('filtering')">Filtering</button>
</div>

<!-- ASSEMBLY -->
<div class="tab-panel active" id="tab-assembly">
  <div class="controls">
    <div class="control-left">
      <button class="reset-btn" onclick="resetTab('assembly')">&#10005; Reset All</button>
    </div>
    <div class="control-right">
      <label class="rows-wrap">
        <span>Rows per page</span>
        <select onchange="setRowsPerPage('assembly', Number(this.value))">
          <option value="20" selected>20</option>
          <option value="50">50</option>
          <option value="100">100</option>
        </select>
      </label>
    </div>
  </div>
  <div class="ms-badge-list" id="badges-assembly"></div>
  <div class="table-wrap">
    <table>
      <thead>
        <tr class="header-row">
          <th onclick="sortTab('assembly','sample_id')"    data-col="sample_id">Sample ID <span class="si">&#8645;</span></th>
          <th onclick="sortTab('assembly','coverage')"     data-col="coverage">Coverage <span class="si">&#8645;</span></th>
          <th onclick="sortTab('assembly','assembly_mode')" data-col="assembly_mode">Mode <span class="si">&#8645;</span></th>
          <th onclick="sortTab('assembly','kmers')"        data-col="kmers">K-mers (spades) <span class="si">&#8645;</span></th>
          <th onclick="sortTab('assembly','asm1_contigs')" data-col="asm1_contigs">Asm1 Contigs <span class="si">&#8645;</span></th>
          <th onclick="sortTab('assembly','asm1_length')"  data-col="asm1_length">Asm1 Length <span class="si">&#8645;</span></th>
          <th onclick="sortTab('assembly','asm1_n50')"     data-col="asm1_n50">Asm1 N50 <span class="si">&#8645;</span></th>
          <th onclick="sortTab('assembly','asm2_contigs')" data-col="asm2_contigs">Asm2 Contigs <span class="si">&#8645;</span></th>
          <th onclick="sortTab('assembly','asm2_length')"  data-col="asm2_length">Asm2 Length <span class="si">&#8645;</span></th>
          <th onclick="sortTab('assembly','asm2_n50')"     data-col="asm2_n50">Asm2 N50 <span class="si">&#8645;</span></th>
          <th onclick="sortTab('assembly','best_assembly')" data-col="best_assembly">Best Assembly<span class="si">&#8645;</span></th>
          <th>Reason</th>
        </tr>
        <tr class="filter-row" id="frow-assembly"></tr>
      </thead>
      <tbody id="asm-body"></tbody>
    </table>
    <div class="no-data" id="asm-nodata" style="display:none">No results match the current filters.</div>
  </div>
  <div class="pager" id="pager-assembly"></div>
</div>

<!-- FILTERING -->
<div class="tab-panel" id="tab-filtering">
  <div class="controls">
    <div class="control-left">
      <button class="reset-btn" onclick="resetTab('filtering')">&#10005; Reset All</button>
    </div>
    <div class="control-right">
      <label class="rows-wrap">
        <span>Rows per page</span>
        <select onchange="setRowsPerPage('filtering', Number(this.value))">
          <option value="20" selected>20</option>
          <option value="50">50</option>
          <option value="100">100</option>
        </select>
      </label>
    </div>
  </div>
  <div class="ms-badge-list" id="badges-filtering"></div>
  <div class="table-wrap">
    <table>
      <thead>
        <tr class="header-row">
          <th onclick="sortTab('filtering','sample_id')"           data-col="sample_id">Sample ID <span class="si">&#8645;</span></th>
          <th onclick="sortTab('filtering','filt_min_length')"     data-col="filt_min_length">Min Length <span class="si">&#8645;</span></th>
          <th onclick="sortTab('filtering','filt_min_cov')"        data-col="filt_min_cov">Min Cov <span class="si">&#8645;</span></th>
          <th onclick="sortTab('filtering','filt_input_contigs')"  data-col="filt_input_contigs">Polished Contigs <span class="si">&#8645;</span></th>
          <th onclick="sortTab('filtering','filt_output_contigs')" data-col="filt_output_contigs">Filtered Contigs <span class="si">&#8645;</span></th>
          <th onclick="sortTab('filtering','filt_input_bp')"       data-col="filt_input_bp">Polished Size <span class="si">&#8645;</span></th>
          <th onclick="sortTab('filtering','filt_output_bp')"      data-col="filt_output_bp">Filtered Size <span class="si">&#8645;</span></th>
          <th onclick="sortTab('filtering','filt_retention')"      data-col="filt_retention">Retention <span class="si">&#8645;</span></th>
        </tr>
        <tr class="filter-row" id="frow-filtering"></tr>
      </thead>
      <tbody id="filt-body"></tbody>
    </table>
    <div class="no-data" id="filt-nodata" style="display:none">No results match the current filters.</div>
  </div>
  <div class="pager" id="pager-filtering"></div>
</div>

<!-- POLISHING -->
<div class="tab-panel" id="tab-polishing">
  <div class="controls">
    <div class="control-left">
      <button class="reset-btn" onclick="resetTab('polishing')">&#10005; Reset All</button>
    </div>
    <div class="control-right">
      <label class="rows-wrap">
        <span>Rows per page</span>
        <select onchange="setRowsPerPage('polishing', Number(this.value))">
          <option value="20" selected>20</option>
          <option value="50">50</option>
          <option value="100">100</option>
        </select>
      </label>
    </div>
  </div>
  <div class="ms-badge-list" id="badges-polishing"></div>
  <div class="table-wrap">
    <table>
      <thead>
        <tr class="header-row">
          <th onclick="sortTab('polishing','sample_id')"                  data-col="sample_id">Sample ID <span class="si">&#8645;</span></th>
          <th onclick="sortTab('polishing','pol_output_contigs')"         data-col="pol_output_contigs">Initial Contigs <span class="si">&#8645;</span></th>
          <th onclick="sortTab('polishing','pol_output_length')"          data-col="pol_output_length">Polished Length <span class="si">&#8645;</span></th>
          <th onclick="sortTab('polishing','pol_total_reads')"            data-col="pol_total_reads">Mapped Reads <span class="si">&#8645;</span></th>
          <th onclick="sortTab('polishing','pol_mean_depth')"             data-col="pol_mean_depth">Mean Depth <span class="si">&#8645;</span></th>
          <th onclick="sortTab('polishing','pol_total_polished_contigs')" data-col="pol_total_polished_contigs">Polished Contigs <span class="si">&#8645;</span></th>
          <th onclick="sortTab('polishing','pol_total_polished')"         data-col="pol_total_polished">Bases Changed <span class="si">&#8645;</span></th>
          <th onclick="sortTab('polishing','pol_total_sub')"              data-col="pol_total_sub">Subs <span class="si">&#8645;</span></th>
          <th onclick="sortTab('polishing','pol_total_ins')"              data-col="pol_total_ins">Ins <span class="si">&#8645;</span></th>
          <th onclick="sortTab('polishing','pol_total_del')"              data-col="pol_total_del">Del <span class="si">&#8645;</span></th>
          <th onclick="sortTab('polishing','pol_contigs_gap_filled')"     data-col="pol_contigs_gap_filled">Gap-Filled Ctgs <span class="si">&#8645;</span></th>
          <th onclick="sortTab('polishing','pol_total_gaps_filled')"      data-col="pol_total_gaps_filled">Gaps Filled <span class="si">&#8645;</span></th>
          <th onclick="sortTab('polishing','pol_total_bases_filled')"     data-col="pol_total_bases_filled">Bases Orig/Filled <span class="si">&#8645;</span></th>
          <th onclick="sortTab('polishing','pol_contigs_rearranged')"     data-col="pol_contigs_rearranged">Misassemblies Contigs <span class="si">&#8645;</span></th>
          <th onclick="sortTab('polishing','pol_total_local_fixes')"      data-col="pol_total_local_fixes">Misassemblies Fixes <span class="si">&#8645;</span></th>
          <th onclick="sortTab('polishing','pol_success')"                data-col="pol_success">Status <span class="si">&#8645;</span></th>
        </tr>
        <tr class="filter-row" id="frow-polishing"></tr>
      </thead>
      <tbody id="pol-body"></tbody>
    </table>
    <div class="no-data" id="pol-nodata" style="display:none">No results match the current filters.</div>
  </div>
  <div class="pager" id="pager-polishing"></div>
</div>

<script>
const RAW = {data_json};
const MAX_COV    = Math.max(...RAW.map(r=>r.coverage||0), 1);
const MAX_DEPTH  = Math.max(...RAW.map(r=>r.pol_mean_depth||0), 1);
const MAX_RETAIN = 100;
const MAX_CHANGES= Math.max(...RAW.map(r=>r.pol_total_polished||0), 1);
const TAB_ORDER  = ['assembly','polishing','filtering'];

const TAB_COLS = {{
  assembly: [
    'sample_id','coverage','assembly_mode','kmers',
    'asm1_contigs','asm1_length','asm1_n50',
    'asm2_contigs','asm2_length','asm2_n50',
    'best_assembly', null
  ],
  filtering: [
    'sample_id','filt_min_length','filt_min_cov',
    'filt_input_contigs','filt_output_contigs',
    'filt_input_bp','filt_output_bp','filt_retention'
  ],
  polishing: [
    'sample_id',
    'pol_output_contigs','pol_output_length',
    'pol_total_reads','pol_mean_depth',
    'pol_total_polished_contigs','pol_total_polished',
    'pol_total_sub','pol_total_ins','pol_total_del',
    'pol_contigs_gap_filled','pol_total_gaps_filled','pol_total_bases_filled',
    'pol_contigs_rearranged','pol_total_local_fixes',
    'pol_success'
  ],
}};

const state = {{
  assembly: {{ filtered:[...RAW], sortCol:'sample_id', sortAsc:true, sel:{{}}, page:1, rowsPerPage:20 }},
  filtering:{{ filtered:[...RAW], sortCol:'sample_id', sortAsc:true, sel:{{}}, page:1, rowsPerPage:20 }},
  polishing:{{ filtered:[...RAW], sortCol:'sample_id', sortAsc:true, sel:{{}}, page:1, rowsPerPage:20 }},
}};
let activeTab = 'assembly';

function getInitialAssemblyMetric(r, asm1Key, asm2Key) {{
  const best = Number(r.best_assembly);
  if (best === 1) return r[asm1Key] ?? r[asm2Key] ?? null;
  if (best === 2) return r[asm2Key] ?? r[asm1Key] ?? null;
  return r[asm1Key] ?? r[asm2Key] ?? null;
}}

function getInitialAssemblyLength(r) {{
  return getInitialAssemblyMetric(r, 'asm1_length', 'asm2_length');
}}

function getInitialAssemblyContigs(r) {{
  return getInitialAssemblyMetric(r, 'asm1_contigs', 'asm2_contigs');
}}

function getFilteringInputContigs(r) {{
  return r.pol_output_contigs ?? r.filt_input_contigs ?? getInitialAssemblyContigs(r);
}}

function getFilteringInputBp(r) {{
  return r.pol_output_length ?? r.filt_input_bp ?? getInitialAssemblyLength(r);
}}

function getRecordValue(r, col) {{
  if (col === 'filt_input_contigs') return getFilteringInputContigs(r);
  if (col === 'filt_input_bp') return getFilteringInputBp(r);
  return r[col];
}}

function colValToDisplay(col, v) {{
  if (col === 'best_assembly') {{
    if (v === 1 || v === '1') return 'Asm 1';
    if (v === 2 || v === '2') return 'Asm 2';
    return '';
  }}
  if (col === 'pol_success') return (v === true || v === 'true') ? 'Success' : 'N/A';
  return (v === null || v === undefined) ? '' : String(v);
}}

function uniqueVals(col) {{
  const vals = new Set();
  RAW.forEach(r => {{
    const v = colValToDisplay(col, getRecordValue(r, col));
    if (v !== '') vals.add(v);
  }});
  return [...vals].sort((a,b) => {{
    const na = Number(a), nb = Number(b);
    if (!isNaN(na) && !isNaN(nb)) return na - nb;
    return a.localeCompare(b);
  }});
}}

function buildFilterRows() {{
  TAB_ORDER.forEach(tab => {{
    const row = document.getElementById('frow-' + tab);
    const cols = TAB_COLS[tab];
    row.innerHTML = '';
    state[tab].sel = {{}};
    cols.forEach(col => {{
      const td = document.createElement('td');
      td.className = 'fc';
      if (!col) {{ row.appendChild(td); return; }}
      const vals = uniqueVals(col);
      if (vals.length === 0) {{ row.appendChild(td); return; }}
      state[tab].sel[col] = new Set();
      const id = `ms-${{tab}}-${{col}}`;
      td.innerHTML = `
        <div class="ms-wrap" id="wrap-${{id}}">
          <div class="ms-trigger" id="trigger-${{id}}" onclick="toggleMS('${{id}}','${{tab}}','${{col}}')">
            <span id="label-${{id}}">All</span><span class="arr">&#9660;</span>
          </div>
          <div class="ms-dropdown" id="dd-${{id}}">
            <div class="ms-search-row">
              <input type="text" placeholder="Search..." oninput="searchMS('${{id}}',this.value)"/>
            </div>
            <div class="ms-actions">
              <button onclick="msSelectAll('${{id}}','${{tab}}','${{col}}')">All</button>
              <button onclick="msClearCol('${{id}}','${{tab}}','${{col}}')">Clear</button>
            </div>
            <div id="opts-${{id}}">
              ${{vals.map(v => `<label class="ms-opt">
                <input type="checkbox" value="${{v}}"
                  onchange="msChange('${{id}}','${{tab}}','${{col}}')"/>${{v}}
              </label>`).join('')}}
            </div>
          </div>
        </div>`;
      row.appendChild(td);
    }});
  }});
}}

document.addEventListener('click', e => {{
  if (!e.target.closest('.ms-wrap')) {{
    document.querySelectorAll('.ms-dropdown.open').forEach(d => {{
      d.classList.remove('open');
      d.previousElementSibling && d.previousElementSibling.classList.remove('open');
    }});
  }}
}});

function toggleMS(id, tab, col) {{
  const dd = document.getElementById('dd-' + id);
  const tr = document.getElementById('trigger-' + id);
  const isOpen = dd.classList.contains('open');
  document.querySelectorAll('.ms-dropdown.open').forEach(d => {{
    d.classList.remove('open');
    const t = document.getElementById(d.id.replace('dd-','trigger-'));
    if (t) t.classList.remove('open');
  }});
  if (!isOpen) {{ dd.classList.add('open'); tr.classList.add('open'); }}
}}

function searchMS(id, q) {{
  q = q.toLowerCase();
  document.querySelectorAll(`#opts-${{id}} .ms-opt`).forEach(opt => {{
    opt.style.display = opt.querySelector('input').value.toLowerCase().includes(q) ? '' : 'none';
  }});
}}

function msChange(id, tab, col) {{
  const checked = [...document.querySelectorAll(`#opts-${{id}} input:checked`)].map(c=>c.value);
  state[tab].sel[col] = new Set(checked);
  updateMSLabel(id, tab, col);
  applyFilters(tab);
}}

function msSelectAll(id, tab, col) {{
  document.querySelectorAll(`#opts-${{id}} .ms-opt`).forEach(opt => {{
    if (opt.style.display !== 'none') opt.querySelector('input').checked = true;
  }});
  msChange(id, tab, col);
}}

function msClearCol(id, tab, col) {{
  document.querySelectorAll(`#opts-${{id}} input`).forEach(c => c.checked = false);
  msChange(id, tab, col);
}}

function updateMSLabel(id, tab, col) {{
  const sel = state[tab].sel[col];
  const lbl = document.getElementById('label-' + id);
  if (!lbl) return;
  if (sel.size === 0) lbl.textContent = 'All';
  else if (sel.size === 1) lbl.textContent = [...sel][0];
  else lbl.textContent = `${{sel.size}} selected`;
}}

function renderBadges(tab) {{
  const container = document.getElementById('badges-' + tab);
  const frags = [];
  for (const [col, selSet] of Object.entries(state[tab].sel)) {{
    if (selSet.size === 0) continue;
    [...selSet].forEach(v => {{
      frags.push(`<div class="ms-badge">${{col}}: ${{v}}
        <span class="x" onclick="removeBadge('${{tab}}','${{col}}','${{v.replace(/'/g,"\\\\'")}}')")>&#10005;</span>
      </div>`);
    }});
  }}
  container.innerHTML = frags.join('');
}}

function removeBadge(tab, col, val) {{
  state[tab].sel[col].delete(val);
  const id = `ms-${{tab}}-${{col}}`;
  const cb = document.querySelector(`#opts-${{id}} input[value="${{val}}"]`);
  if (cb) cb.checked = false;
  updateMSLabel(id, tab, col);
  applyFilters(tab);
}}

function applyFilters(tab) {{
  let data = [...RAW];
  for (const [col, selSet] of Object.entries(state[tab].sel)) {{
    if (selSet.size === 0) continue;
    data = data.filter(r => selSet.has(colValToDisplay(col, getRecordValue(r, col))));
  }}
  state[tab].filtered = data;
  state[tab].page = 1;
  renderBadges(tab);
  sortAndRender(tab);
}}

function resetTab(tab) {{
  const cols = TAB_COLS[tab].filter(Boolean);
  cols.forEach(col => {{
    state[tab].sel[col] = new Set();
    const id = `ms-${{tab}}-${{col}}`;
    document.querySelectorAll(`#opts-${{id}} input`).forEach(c => c.checked = false);
    updateMSLabel(id, tab, col);
  }});
  applyFilters(tab);
}}

function switchTab(name) {{
  activeTab = name;
  document.querySelectorAll('.tab-btn').forEach((b,i) => {{
    b.classList.toggle('active', TAB_ORDER[i]===name);
  }});
  document.querySelectorAll('.tab-panel').forEach(p => p.classList.remove('active'));
  document.getElementById('tab-'+name).classList.add('active');
  renderStats();
}}

function dash() {{ return '<span style="color:var(--muted)">&#8212;</span>'; }}
function fmtN(n, dec=0, suf='') {{
  if (n===null||n===undefined) return dash();
  return Number(n).toLocaleString(undefined,{{minimumFractionDigits:dec,maximumFractionDigits:dec}})+suf;
}}
function bar(val, max, color) {{
  const pct = Math.min(((val||0)/max)*100,100);
  return `<div class="bar-wrap">
    <div class="bar-track"><div class="bar-fill" style="width:${{pct.toFixed(1)}}%;background:${{color}}"></div></div>
    <span class="bar-num">${{fmtN(val,2)}}</span></div>`;
}}

function renderStats() {{
  const total    = RAW.length;
  const polished = RAW.filter(r=>r.pol_success).length;
  const avgCov   = RAW.length ? (RAW.reduce((s,r)=>s+(r.coverage||0),0)/RAW.length).toFixed(1) : 0;
  const retSrc   = RAW.filter(r=>r.filt_retention!==null);
  const avgRet   = retSrc.length
    ? (retSrc.reduce((s,r)=>s+(r.filt_retention||0),0)/retSrc.length).toFixed(1)+'%' : '&#8212;';
  const totalChanges = RAW.reduce((s,r)=>s+(r.pol_total_polished||0),0);
  document.getElementById('statsBar').innerHTML = `
    <div class="stat-card"><div class="stat-label">Samples</div><div class="stat-value">${{total}}</div></div>
    <div class="stat-card"><div class="stat-label">Avg Coverage</div><div class="stat-value">${{avgCov}}X</div></div>
    <div class="stat-card"><div class="stat-label">Polished</div><div class="stat-value">${{polished}}</div></div>
    <div class="stat-card"><div class="stat-label">Avg Retention</div><div class="stat-value">${{avgRet}}</div></div>
    <div class="stat-card"><div class="stat-label">Total Corrections</div><div class="stat-value">${{totalChanges.toLocaleString()}}</div></div>
  `;
}}

function renderAssembly(data) {{
  const tbody  = document.getElementById('asm-body');
  const noData = document.getElementById('asm-nodata');
  if (!data.length) {{ tbody.innerHTML=''; noData.style.display='block'; return; }}
  noData.style.display='none';
  tbody.innerHTML = data.map(r => {{
    const covBar    = bar(r.coverage, MAX_COV, 'var(--accent)');
    const bestBadge = r.best_assembly===1
      ? '<span class="badge b-blue">Asm 1</span>'
      : r.best_assembly===2 ? '<span class="badge b-orange">Asm 2</span>' : dash();
    const modeBadge = r.assembly_mode
      ? `<span class="badge b-grey">${{r.assembly_mode}}</span>` : dash();
    const reasons   = (r.best_reasons||[]).map(s=>s.replace(/^[^a-zA-Z0-9(]+/,'')).join('; ') || '&#8212;';
    return `<tr>
      <td class="sc">${{r.sample_id}}</td>
      <td>${{covBar}}</td>
      <td>${{modeBadge}}</td>
      <td style="font-size:.8rem;color:var(--muted)">${{r.kmers||dash()}}</td>
      <td>${{fmtN(r.asm1_contigs)}}</td><td>${{fmtN(r.asm1_length)}} bp</td><td>${{fmtN(r.asm1_n50)}} bp</td>
      <td>${{fmtN(r.asm2_contigs)}}</td><td>${{fmtN(r.asm2_length)}} bp</td><td>${{fmtN(r.asm2_n50)}} bp</td>
      <td>${{bestBadge}}</td>
      <td style="white-space:normal;font-size:.8rem;max-width:200px">${{reasons}}</td>
    </tr>`;
  }}).join('');
}}

function renderFiltering(data) {{
  const tbody  = document.getElementById('filt-body');
  const noData = document.getElementById('filt-nodata');
  if (!data.length) {{ tbody.innerHTML=''; noData.style.display='block'; return; }}
  noData.style.display='none';
  tbody.innerHTML = data.map(r => {{
    const inputContigs = getFilteringInputContigs(r);
    const inputBp = getFilteringInputBp(r);
    const retBar  = r.filt_retention!==null ? bar(r.filt_retention, MAX_RETAIN, 'var(--accent3)') : dash();
    const dropped = (inputContigs!==null && r.filt_output_contigs!==null)
      ? inputContigs - r.filt_output_contigs : null;
    const droppedStr = dropped !== null
        ? ` <span style="font-size:.72rem;color:var(--muted)">(&#8722;${{dropped}})</span>`
        : '';
    const outDiff = (r.filt_output_bp != null && inputBp != null)
      ? r.filt_output_bp - inputBp : null;
    const outDiffStr = outDiff !== null
        ? ` <span style="font-size:.72rem;color:var(--muted)">(${{outDiff >= 0 ? '+' : ''}}${{outDiff.toLocaleString()}} bp)</span>`
        : '';
    return `<tr>
      <td class="sc">${{r.sample_id}}</td>
      <td>${{fmtN(r.filt_min_length)}} bp</td>
      <td>${{fmtN(r.filt_min_cov,1)}}X</td>
      <td>${{fmtN(inputContigs)}}</td>
      <td>${{fmtN(r.filt_output_contigs)}}${{droppedStr}}</td>
      <td>${{fmtN(inputBp)}} bp</td>
      <td>${{fmtN(r.filt_output_bp)}} bp${{outDiffStr}}</td>
      <td style="min-width:130px">${{retBar}}</td>
    </tr>`;
  }}).join('');
}}

function renderPolishing(data) {{
  const tbody  = document.getElementById('pol-body');
  const noData = document.getElementById('pol-nodata');
  if (!data.length) {{ tbody.innerHTML=''; noData.style.display='block'; return; }}
  noData.style.display='none';
  tbody.innerHTML = data.map(r => {{
    const depBar  = r.pol_mean_depth!==null ? bar(r.pol_mean_depth, MAX_DEPTH, 'var(--accent4)') : dash();
    const chgBar  = r.pol_total_polished
      ? `<div class="bar-wrap">
          <div class="bar-track"><div class="bar-fill" style="width:${{Math.min(((r.pol_total_polished||0)/MAX_CHANGES)*100,100).toFixed(1)}}%;background:var(--accent3)"></div></div>
          <span class="bar-num">${{fmtN(r.pol_total_polished,0)}}</span></div>`
      : '<span class="badge b-grey">0</span>';
    const status  = r.pol_success
      ? '<span class="badge b-yes">&#10003; Success</span>'
      : '<span class="badge b-no">&#8212; N/A</span>';
    const corrStr = (r.pol_total_polished_contigs!==null && r.pol_total_contigs!==null)
      ? `${{r.pol_total_polished_contigs}} / ${{r.pol_total_contigs}}`
      : fmtN(r.pol_total_polished_contigs);
    const subBadge = r.pol_total_sub!=null ? `<span class="badge b-sub">S: ${{r.pol_total_sub}}</span>`  : dash();
    const insBadge = r.pol_total_ins!=null ? `<span class="badge b-ins">I: ${{r.pol_total_ins}}</span>`  : dash();
    const delBadge = r.pol_total_del!=null ? `<span class="badge b-del">D: ${{r.pol_total_del}}</span>`  : dash();
    const initialLength = getInitialAssemblyLength(r);
    const lenDiff   = (r.pol_output_length != null && initialLength != null)
      ? r.pol_output_length - initialLength : null;
    const diffStr   = lenDiff !== null
      ? ` <span style="font-size:.72rem;color:var(--muted)">(${{lenDiff >= 0 ? '+' : ''}}${{lenDiff.toLocaleString()}} bp)</span>`
      : '';
    const gapFilled  = r.pol_contigs_gap_filled!=null
      ? `<span class="badge b-ins">${{r.pol_contigs_gap_filled}}</span>` : dash();
    const totalGaps  = fmtN(r.pol_total_gaps_filled);
    const basesFilled= r.pol_total_bases_filled_pair ?? fmtN(r.pol_total_bases_filled);
    const rearranged = r.pol_contigs_rearranged!=null
      ? `<span class="badge b-orange">${{r.pol_contigs_rearranged}}</span>` : dash();
    const localFixes = fmtN(r.pol_total_local_fixes);
    return `<tr>
      <td class="sc">${{r.sample_id}}</td>
      <td>${{fmtN(r.pol_output_contigs)}}</td>
      <td>${{fmtN(r.pol_output_length)}} bp${{diffStr}}</td>
      <td>${{fmtN(r.pol_total_reads)}}</td>
      <td style="min-width:130px">${{depBar}}</td>
      <td>${{corrStr}}</td>
      <td style="min-width:130px">${{chgBar}}</td>
      <td>${{subBadge}}</td>
      <td>${{insBadge}}</td>
      <td>${{delBadge}}</td>
      <td>${{gapFilled}}</td>
      <td>${{totalGaps}}</td>
      <td>${{basesFilled}}</td>
      <td>${{rearranged}}</td>
      <td>${{localFixes}}</td>
      <td>${{status}}</td>
    </tr>`;
  }}).join('');
}}

const RENDERERS = {{ assembly:renderAssembly, filtering:renderFiltering, polishing:renderPolishing }};

function sortTab(tab, col) {{
  const s = state[tab];
  if (s.sortCol===col) s.sortAsc=!s.sortAsc; else {{ s.sortCol=col; s.sortAsc=true; }}
  document.querySelectorAll('#tab-'+tab+' th[data-col]').forEach(th => {{
    th.classList.toggle('sorted', th.dataset.col===col);
    th.querySelector('.si').innerHTML = th.dataset.col===col?(s.sortAsc?'&#8593;':'&#8595;'):'&#8645;';
  }});
  sortAndRender(tab);
}}

function setRowsPerPage(tab, value) {{
  state[tab].rowsPerPage = value;
  state[tab].page = 1;
  sortAndRender(tab);
}}

function changePage(tab, step) {{
  const s = state[tab];
  const totalPages = Math.max(1, Math.ceil(s.filtered.length / s.rowsPerPage));
  s.page = Math.min(totalPages, Math.max(1, s.page + step));
  sortAndRender(tab);
}}

function renderPager(tab) {{
  const s = state[tab];
  const totalRows = s.filtered.length;
  const pager = document.getElementById('pager-' + tab);
  if (!pager) return;
  if (!totalRows) {{ pager.innerHTML = ''; return; }}

  const totalPages = Math.max(1, Math.ceil(totalRows / s.rowsPerPage));
  if (s.page > totalPages) s.page = totalPages;

  const startRow = ((s.page - 1) * s.rowsPerPage) + 1;
  const endRow = Math.min(s.page * s.rowsPerPage, totalRows);

  pager.innerHTML = `
    <span class="pager-info">Rows ${{startRow}}-${{endRow}} of ${{totalRows}}</span>
    <button type="button" onclick="changePage('${{tab}}', -1)" ${{s.page === 1 ? 'disabled' : ''}}>Prev</button>
    <span class="pager-info">Page ${{s.page}} / ${{totalPages}}</span>
    <button type="button" onclick="changePage('${{tab}}', 1)" ${{s.page === totalPages ? 'disabled' : ''}}>Next</button>
  `;
}}

function initDragScroll() {{
  document.querySelectorAll('.table-wrap').forEach(wrap => {{
    let isDragging = false;
    let startX = 0;
    let startScrollLeft = 0;
    let dragged = false;

    wrap.addEventListener('pointerdown', e => {{
      if (e.button !== 0) return;
      if (wrap.scrollWidth <= wrap.clientWidth) return;
      if (!e.target.closest('tbody td')) return;
      if (e.target.closest('.ms-dropdown, input, select, button, textarea, a, label')) return;

      isDragging = true;
      dragged = false;
      startX = e.clientX;
      startScrollLeft = wrap.scrollLeft;
      wrap.classList.add('dragging');
      if (wrap.setPointerCapture) wrap.setPointerCapture(e.pointerId);
    }});

    wrap.addEventListener('pointermove', e => {{
      if (!isDragging) return;
      const delta = e.clientX - startX;
      if (Math.abs(delta) > 4) dragged = true;
      wrap.scrollLeft = startScrollLeft - delta;
    }});

    const stopDrag = e => {{
      if (!isDragging) return;
      isDragging = false;
      wrap.classList.remove('dragging');
      if (e && wrap.releasePointerCapture && wrap.hasPointerCapture && wrap.hasPointerCapture(e.pointerId)) {{
        wrap.releasePointerCapture(e.pointerId);
      }}
    }};

    wrap.addEventListener('pointerup', stopDrag);
    wrap.addEventListener('pointercancel', stopDrag);

    wrap.addEventListener('click', e => {{
      if (!dragged) return;
      e.preventDefault();
      e.stopPropagation();
      dragged = false;
    }}, true);
  }});
}}

function sortAndRender(tab) {{
  const s = state[tab];
  s.filtered.sort((a,b) => {{
    let va=getRecordValue(a, s.sortCol), vb=getRecordValue(b, s.sortCol);
    if (va===null||va===undefined) va=s.sortAsc? Infinity:-Infinity;
    if (vb===null||vb===undefined) vb=s.sortAsc? Infinity:-Infinity;
    if (typeof va==='boolean') return s.sortAsc?(va===vb?0:va?1:-1):(va===vb?0:va?-1:1);
    if (typeof va==='number'&&typeof vb==='number') return s.sortAsc?va-vb:vb-va;
    return s.sortAsc?String(va).localeCompare(String(vb)):String(vb).localeCompare(String(va));
  }});

  const totalPages = Math.max(1, Math.ceil(s.filtered.length / s.rowsPerPage));
  if (s.page > totalPages) s.page = totalPages;

  const start = (s.page - 1) * s.rowsPerPage;
  const pageRows = s.filtered.slice(start, start + s.rowsPerPage);

  RENDERERS[tab](pageRows);
  renderPager(tab);
  renderStats();
}}

const HEADS = {{
  assembly:  ['Sample ID','Coverage (X)','Mode','K-mers (spades)',
               'Asm1 Contigs','Asm1 Length (bp)','Asm1 N50 (bp)',
               'Asm2 Contigs','Asm2 Length (bp)','Asm2 N50 (bp)',
               'Best Assembly','Reasons'],
  filtering: ['Sample ID','Min Length (bp)','Min Cov (X)',
               'Polished Contigs','Filtered Contigs',
               'Polished Size (bp)','Filtered Size (bp)','Retention (%)'],
  polishing: ['Sample ID','Initial Contigs','Polished Length (bp)',
               'Mapped Reads','Mean Depth (X)',
               'Polished Contigs','Bases Changed',
               'Subs','Ins','Del',
               'Gap-Filled Ctgs','Gaps Filled','Bases Orig/Filled',
               'Misassemblies Contigs','Misassemblies Fixes',
               'Status'],
}};

function getRows(tab) {{
  return state[tab].filtered.map(r => {{
    if (tab==='assembly') return [
      r.sample_id,r.coverage??'',r.assembly_mode??'',r.kmers??'',
      r.asm1_contigs??'',r.asm1_length??'',r.asm1_n50??'',
      r.asm2_contigs??'',r.asm2_length??'',r.asm2_n50??'',
      r.best_assembly??'',(r.best_reasons||[]).join('; '),
    ];
    if (tab==='filtering') {{
      const inputContigs = getFilteringInputContigs(r);
      const inputBp = getFilteringInputBp(r);
      const outDiffExp = (r.filt_output_bp != null && inputBp != null)
        ? r.filt_output_bp - inputBp : null;
      const outBpExport = r.filt_output_bp != null
        ? `${{r.filt_output_bp}}${{outDiffExp !== null ? ` (${{outDiffExp >= 0 ? '+' : ''}}${{outDiffExp}})` : ''}}`
        : '';
      return [
        r.sample_id,r.filt_min_length??'',r.filt_min_cov??'',
        inputContigs??'',r.filt_output_contigs??'',
        inputBp??'', outBpExport, r.filt_retention??'',
      ];
    }}
    // polishing
    const corrStr = (r.pol_total_polished_contigs!=null && r.pol_total_contigs!=null)
      ? `${{r.pol_total_polished_contigs}} / ${{r.pol_total_contigs}}`
      : (r.pol_total_polished_contigs??'');
    const initialLength = getInitialAssemblyLength(r);
    const lenDiff = (r.pol_output_length != null && initialLength != null)
      ? r.pol_output_length - initialLength : null;
    const lenExport = r.pol_output_length != null
      ? (lenDiff !== null
          ? `${{r.pol_output_length}} (${{lenDiff >= 0 ? '+' : ''}}${{lenDiff}})`
          : String(r.pol_output_length))
      : '';
    const basesFilledExport = r.pol_total_bases_filled_pair ?? (r.pol_total_bases_filled ?? '');
    return [
      r.sample_id,
      r.pol_output_contigs??'', lenExport,
      r.pol_total_reads??'', r.pol_mean_depth??'',
      corrStr, r.pol_total_polished??'',
      r.pol_total_sub??'', r.pol_total_ins??'', r.pol_total_del??'',
      r.pol_contigs_gap_filled??'', r.pol_total_gaps_filled??'', basesFilledExport,
      r.pol_contigs_rearranged??'', r.pol_total_local_fixes??'',
      r.pol_success?'Success':'N/A',
    ];
  }});
}}

function exportCSV() {{
  const rows=[HEADS[activeTab],...getRows(activeTab)];
  downloadBlob(rows.map(r=>r.map(v=>`"${{String(v).replace(/"/g,'""')}}"`).join(',')).join('\\n'),
    `assembly_${{activeTab}}_report.csv`,'text/csv');
}}
function exportTSV() {{
  const rows=[HEADS[activeTab],...getRows(activeTab)];
  downloadBlob(rows.map(r=>r.join('\\t')).join('\\n'),
    `assembly_${{activeTab}}_report.tsv`,'text/tab-separated-values');
}}
function exportExcel() {{
  const wb=XLSX.utils.book_new();
  TAB_ORDER.forEach(tab=>{{
    const ws=XLSX.utils.aoa_to_sheet([HEADS[tab],...getRows(tab)]);
    ws['!cols']=HEADS[tab].map(()=>({{wch:20}}));
    XLSX.utils.book_append_sheet(wb,ws,tab.charAt(0).toUpperCase()+tab.slice(1));
  }});
  XLSX.writeFile(wb,'assembly_report.xlsx');
}}
function exportPDF() {{
  const {{jsPDF}}=window.jspdf;
  const doc=new jsPDF({{orientation:'landscape',unit:'mm',format:'a4'}});
  doc.setFont('helvetica','bold');doc.setFontSize(14);doc.setTextColor(10,158,116);
  doc.text('Assembly Report',14,16);
  doc.setFont('helvetica','normal');doc.setFontSize(8);doc.setTextColor(108,117,125);
  doc.text(`Generated: ${{new Date().toLocaleString()}} | Tab: ${{activeTab}} | Rows: ${{state[activeTab].filtered.length}}`,14,22);
  doc.autoTable({{startY:27,head:[HEADS[activeTab]],body:getRows(activeTab),
    styles:{{font:'helvetica',fontSize:6,cellPadding:2,overflow:'linebreak'}},
    headStyles:{{fillColor:[10,158,116],textColor:[255,255,255],fontStyle:'bold',fontSize:7}},
    alternateRowStyles:{{fillColor:[245,247,250]}},theme:'grid'}});
  doc.save(`assembly_${{activeTab}}_report.pdf`);
}}
function downloadBlob(c,f,m){{const b=new Blob([c],{{type:m}});const a=document.createElement('a');a.href=URL.createObjectURL(b);a.download=f;a.click();}}

buildFilterRows();
initDragScroll();
TAB_ORDER.forEach(t=>applyFilters(t));
renderStats();
</script>
</body>
</html>"""

    with open(output_html, "w") as fh:
        fh.write(html)
    print(f"Report written to: {output_html}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: assembly_report.py <input_dir> [output.html]", file=sys.stderr)
        sys.exit(1)
    input_dir   = sys.argv[1]
    output_html = sys.argv[2] if len(sys.argv) > 2 else "assembly_report.html"
    main(input_dir, output_html)
