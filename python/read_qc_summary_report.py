#!/usr/bin/env python3

import os
import sys
import glob
import json


def parse_qc_summary(filepath):
    data = {
        "sample_id":          "",
        "total_reads":        None,
        "total_bases":        None,
        "mean_read_length":   None,
        "genome_size_method": "",
        "genome_size":        None,
        "estimated_coverage": None,
        "downsampled":        False,
        "downsample_factor":  None,
        "target_coverage":    None,
        "warnings":           [],
        "qc_status":          None,
        "qc_fail_reason":     None,
    }

    with open(filepath) as f:
        content = f.read()

    section = None
    genome_metrics_complete = False

    for line in content.splitlines():
        s = line.strip()

        if not s:
            continue

        if s.startswith("====") and s.endswith("===="):
            if "FASTP AFTER FILTERING" in s:
                section = "fastp"
            elif "GENOME SIZE" in s:
                section = "genome_size"
                genome_metrics_complete = False
            elif "DOWNSAMPLING" in s:
                section = "downsampling"
            else:
                section = None
            continue

        if s.startswith("Sample ID:"):
            data["sample_id"] = s.split(":", 1)[1].strip()
        elif section == "fastp":
            if s.startswith("Total reads:"):
                data["total_reads"] = int(s.split(":", 1)[1].strip().replace(",", ""))
            elif s.startswith("Total bases:"):
                data["total_bases"] = int(s.split(":", 1)[1].strip().replace(",", ""))
            elif s.startswith("Mean read length:"):
                data["mean_read_length"] = float(s.split(":", 1)[1].strip())
        elif section == "genome_size":
            if s.startswith("Method:"):
                data["genome_size_method"] = s.split(":", 1)[1].strip()
            elif s.startswith("Genome size:"):
                val = s.split(":", 1)[1].strip().replace(",", "").replace(" bp", "").strip()
                data["genome_size"] = int(float(val))
            elif s.startswith("Estimated coverage:") and not genome_metrics_complete:
                val = s.split(":", 1)[1].strip().split()[0].replace("X", "").strip()
                data["estimated_coverage"] = float(val)
                genome_metrics_complete = True
            elif genome_metrics_complete or s.startswith("\u26a0") or s.startswith("\u2192"):
                data["warnings"].append(s)
        elif section == "downsampling":
            if s.startswith("Read coverage exceeds maximum threshold"):
                data["downsampled"] = True
            elif s.startswith("Downsampling factor:"):
                data["downsample_factor"] = float(s.split(":", 1)[1].strip())
            elif s.startswith("Target coverage:"):
                val = s.split(":", 1)[1].strip().replace("X", "").strip()
                data["target_coverage"] = float(val)
        elif s.startswith("\u26a0") or s.startswith("\u2192"):
            data["warnings"].append(s)

    # Inspect qc_status.txt for PASS/FAIL status and reason
    base_dir      = os.path.dirname(filepath)
    basename      = os.path.basename(filepath)
    sample_prefix = basename.replace("_QC_summary.txt", "")
    status_path   = os.path.join(base_dir, f"{sample_prefix}_qc_status.txt")

    if os.path.isfile(status_path):
        with open(status_path) as sf:
            status_line = sf.readline().strip()

        if status_line.startswith("PASS"):
            data["qc_status"]     = "PASS"
            data["qc_fail_reason"] = None
        elif status_line.startswith("FAIL"):
            data["qc_status"]     = "FAIL"
            parts = status_line.split("\t", 1)
            data["qc_fail_reason"] = parts[1].strip() if len(parts) > 1 else "Insufficient coverage"
        else:
            data["qc_status"]     = "UNKNOWN"
            data["qc_fail_reason"] = status_line or None
    else:
        data["qc_status"]     = "UNKNOWN"
        data["qc_fail_reason"] = "qc_status.txt not found"

    return data


def main(input_dir, output_html):
    pattern = os.path.join(input_dir, "*_QC_summary.txt")
    files   = sorted(glob.glob(pattern))

    if not files:
        print(f"No *_QC_summary.txt files found in: {input_dir}")
        sys.exit(1)

    records = []
    for f in files:
        rec = parse_qc_summary(f)
        records.append(rec)
        print(f"  Parsed: {rec['sample_id']} | cov={rec['estimated_coverage']}X | "
              f"downsampled={rec['downsampled']} | qc_status={rec['qc_status']}")

    print(f"  Total samples: {len(records)}")
    data_json = json.dumps(records)

    # header label in html report — QC Status inserted after Sample ID
    col_defs = [
        ("sample_id",          "Sample ID"),
        ("qc_status",          "QC Status"),
        ("total_reads",        "Total Reads"),
        ("total_bases",        "Total Bases"),
        ("mean_read_length",   "Mean Read Length"),
        ("genome_size",        "Genome Size"),
        ("genome_size_method", "Method"),
        ("estimated_coverage", "Coverage"),
        ("downsampled",        "Downsampled"),
        ("downsample_factor",  "Downsample Factor"),
        ("target_coverage",    "Target Cov"),
        (None,                 "Fail Reason"),
        (None,                 "Warnings"),
    ]

    # build header cells
    header_cells = ""
    for key, label in col_defs:
        if key:
            header_cells += (
                f'<th onclick="sortBy(\'{key}\')" data-col="{key}">'
                f'{label} <span class="si">&#8597;</span></th>\n'
            )
        else:
            header_cells += f'<th>{label}</th>\n'

    # build filter row with multi-select dropdowns (only for keyed columns)
    filter_cells = ""
    for key, label in col_defs:
        if key is None:
            filter_cells += '<td class="fc"></td>\n'
        else:
            ms_id = f"ms-{key}"
            filter_cells += (
                f'<td class="fc">'
                f'<div class="ms-wrap" id="wrap-{ms_id}">'
                f'<div class="ms-trigger" id="trigger-{ms_id}" onclick="toggleMS(\'{ms_id}\',\'{key}\')">'
                f'<span id="label-{ms_id}">All</span><span class="arr">&#9660;</span>'
                f'</div>'
                f'<div class="ms-dropdown" id="dd-{ms_id}">'
                f'<div class="ms-search-row"><input type="text" placeholder="Search..." oninput="searchMS(\'{ms_id}\',this.value)"/></div>'
                f'<div class="ms-actions">'
                f'<button onclick="msAll(\'{ms_id}\',\'{key}\')">All</button>'
                f'<button onclick="msClear(\'{ms_id}\',\'{key}\')">Clear</button>'
                f'</div>'
                f'<div id="opts-{ms_id}"></div>'
                f'</div></div></td>\n'
            )

    css = """
  :root {
    --bg:#f8f9fa; --surface:#ffffff; --surface2:#f0f2f5; --border:#dee2e6;
    --accent:#0a9e74; --accent2:#e53935; --warn:#f59e0b;
    --pass:#0a9e74; --fail:#e53935; --unknown:#6c757d;
    --text:#1a1f2e; --muted:#6c757d;
    --font-main: Calibri, 'Calibri', Arial, sans-serif;
  }
  *,*::before,*::after{box-sizing:border-box;margin:0;padding:0}
  body{background:var(--bg);color:var(--text);font-family:var(--font-main);min-height:100vh}
  header{background:var(--surface);border-bottom:2px solid var(--accent);padding:20px 40px;
    display:flex;align-items:center;justify-content:space-between;gap:16px;flex-wrap:wrap;
    box-shadow:0 2px 8px rgba(0,0,0,.06)}
  .logo{font-family:var(--font-main);font-size:1.2rem;font-weight:700;color:var(--accent);
    letter-spacing:2px;text-transform:uppercase}
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
  .stat-label{font-size:.7rem;font-family:var(--font-main);font-weight:600;color:var(--muted);
    text-transform:uppercase;letter-spacing:1.5px;margin-bottom:6px}
  .stat-value{font-size:1.8rem;font-weight:700;color:var(--accent);font-family:var(--font-main)}
  .stat-value.fail{color:var(--fail)}
  .controls{padding:16px 40px 8px;display:flex;gap:12px;flex-wrap:wrap;align-items:center}
  .reset-btn{background:transparent;border:1px solid var(--border);color:var(--muted);padding:8px 18px;
    border-radius:6px;cursor:pointer;font-family:var(--font-main);font-size:.82rem;font-weight:600;
    letter-spacing:.5px;transition:border-color .2s,color .2s}
  .reset-btn:hover{border-color:var(--accent2);color:var(--accent2)}
  .rows-wrap{display:flex;align-items:center;gap:8px;font-family:var(--font-main);font-size:.8rem;color:var(--muted)}
  .rows-wrap select{background:var(--surface);border:1px solid var(--border);color:var(--text);padding:8px 12px;border-radius:6px;font-family:var(--font-main);font-size:.82rem;cursor:pointer}
  .rows-wrap select:focus{outline:none;border-color:var(--accent)}
  /* ── multi-select ── */
  .ms-wrap{position:relative;min-width:110px}
  .ms-trigger{background:var(--surface);border:1px solid var(--border);color:var(--text);
    padding:4px 8px;border-radius:4px;font-family:var(--font-main);font-size:.78rem;
    cursor:pointer;display:flex;justify-content:space-between;align-items:center;
    gap:4px;user-select:none;transition:border-color .2s;white-space:nowrap;min-width:110px}
  .ms-trigger:hover,.ms-trigger.open{border-color:var(--accent)}
  .ms-trigger .arr{font-size:.55rem;color:var(--muted);transition:transform .2s;flex-shrink:0}
  .ms-trigger.open .arr{transform:rotate(180deg)}
  .ms-dropdown{display:none;position:absolute;top:calc(100% + 2px);left:0;
    background:var(--surface);border:1px solid var(--border);border-radius:6px;
    z-index:500;min-width:180px;max-height:220px;overflow-y:auto;
    box-shadow:0 6px 20px rgba(0,0,0,.12)}
  .ms-dropdown.open{display:block}
  .ms-search-row{padding:6px 8px;border-bottom:1px solid var(--border);
    position:sticky;top:0;background:var(--surface);z-index:1}
  .ms-search-row input{width:100%;background:var(--surface2);border:1px solid var(--border);
    color:var(--text);padding:4px 8px;border-radius:4px;font-size:.78rem;font-family:var(--font-main);
    outline:none;min-width:unset}
  .ms-search-row input:focus{border-color:var(--accent)}
  .ms-actions{display:flex;gap:8px;padding:4px 8px;border-bottom:1px solid var(--border);
    position:sticky;top:37px;background:var(--surface);z-index:1}
  .ms-actions button{background:transparent;border:none;color:var(--accent);
    font-family:var(--font-main);font-size:.72rem;font-weight:600;cursor:pointer;padding:1px 3px}
  .ms-actions button:hover{text-decoration:underline}
  .ms-opt{display:flex;align-items:center;gap:8px;padding:7px 10px;
    cursor:pointer;transition:background .12s;font-size:.82rem;font-family:var(--font-main)}
  .ms-opt:hover{background:var(--surface2)}
  .ms-opt input[type=checkbox]{accent-color:var(--accent);width:13px;height:13px;
    min-width:unset;cursor:pointer;flex-shrink:0}
  .ms-badge-list{display:flex;flex-wrap:wrap;gap:3px;padding:4px 40px 8px;min-height:0}
  .ms-badge{background:rgba(10,158,116,.1);color:var(--accent);
    border:1px solid rgba(10,158,116,.25);border-radius:12px;padding:2px 8px;
    font-size:.72rem;font-family:var(--font-main);display:flex;align-items:center;gap:4px}
  .ms-badge .x{cursor:pointer;opacity:.6;font-size:.65rem}
  .ms-badge .x:hover{opacity:1;color:var(--accent2)}
  /*  table  */
  .table-wrap{padding:0 40px 40px;overflow-x:auto;cursor:grab;scrollbar-gutter:stable;touch-action:pan-y;scrollbar-color:#5f6368 var(--bg);scrollbar-width:thin}
  .table-wrap.dragging{cursor:grabbing;user-select:none}
  table{width:100%;border-collapse:collapse;font-size:.85rem}
  thead tr.header-row{background:var(--accent)}
  thead tr.filter-row{background:var(--surface2);border-bottom:2px solid var(--border)}
  th{padding:11px 13px;text-align:left;font-family:var(--font-main);font-size:.78rem;font-weight:700;
    letter-spacing:.8px;text-transform:uppercase;color:#fff;cursor:pointer;user-select:none;
    white-space:nowrap;border-right:1px solid rgba(255,255,255,.15)}
  th:last-child{border-right:none}
  th:hover{background:rgba(255,255,255,.1)}
  th .si{margin-left:5px;opacity:.75;font-style:normal;display:inline-block}
  th.sorted .si{opacity:1}
  td.fc{padding:4px 5px;vertical-align:top}
  tbody tr{border-bottom:1px solid var(--border);transition:background .12s}
  tbody tr:nth-child(even){background:var(--surface2)}
  tbody tr:hover{background:rgba(10,158,116,.06)}
  tbody tr.row-fail{background:rgba(229,57,53,.04) !important}
  tbody tr.row-fail:hover{background:rgba(229,57,53,.09) !important}
  td{padding:10px 13px;vertical-align:middle;border-right:1px solid var(--border);white-space:nowrap}
  td:last-child{border-right:none}
  td.sample-cell{font-family:var(--font-main);font-size:.85rem;font-weight:700;color:var(--text)}
  /*  QC status badges  */
  .badge{display:inline-block;padding:2px 10px;border-radius:12px;font-size:.75rem;
    font-family:var(--font-main);font-weight:600}
  .badge-pass{background:rgba(10,158,116,.12);color:#056b4f;border:1px solid rgba(10,158,116,.35)}
  .badge-fail{background:rgba(229,57,53,.12);color:#b71c1c;border:1px solid rgba(229,57,53,.35)}
  .badge-unknown{background:rgba(108,117,125,.12);color:#495057;border:1px solid rgba(108,117,125,.3)}
  .badge-yes{background:rgba(10,158,116,.12);color:#000;border:1px solid rgba(10,158,116,.3)}
  .badge-no{background:rgba(222,226,230,.6);color:#000;border:1px solid var(--border)}
  /*  fail reason pill  */
  .fail-reason{display:inline-block;padding:2px 8px;border-radius:4px;font-size:.75rem;
    font-family:var(--font-main);background:rgba(229,57,53,.08);color:#b71c1c;
    border:1px solid rgba(229,57,53,.25);white-space:normal;max-width:260px}
  .cov-bar-wrap{display:flex;align-items:center;gap:8px;min-width:140px}
  .cov-track{flex:1;height:5px;background:var(--border);border-radius:3px;min-width:60px}
  .cov-fill{height:100%;border-radius:3px;background:var(--accent)}
  .cov-num{font-family:var(--font-main);font-size:.78rem;color:var(--text);min-width:58px;text-align:right}
  .warn-list{display:flex;flex-direction:column;gap:3px}
  .warn-badge{display:inline-block;padding:2px 8px;border-radius:4px;font-size:.75rem;
    font-family:var(--font-main);background:rgba(245,158,11,.1);color:#000;
    border:1px solid rgba(245,158,11,.3);white-space:normal;max-width:280px}
  .no-data{text-align:center;padding:60px;color:var(--muted);font-family:var(--font-main);font-size:.9rem}
  .pagination{display:flex;align-items:center;justify-content:flex-end;gap:8px;
    padding:0 40px 40px;flex-wrap:wrap}
  .page-info{font-size:.78rem;font-family:var(--font-main);color:var(--muted);margin-right:12px}
  .page-btn{background:var(--surface);border:1px solid var(--border);color:var(--text);
    padding:6px 14px;border-radius:6px;cursor:pointer;font-family:var(--font-main);
    font-size:.78rem;transition:all .15s}
  .page-btn:hover:not(:disabled){border-color:var(--accent);color:var(--accent)}
  .page-btn:disabled{opacity:.35;cursor:default}
  .page-btn.active{background:var(--accent);color:#fff;border-color:var(--accent);font-weight:700}
  ::-webkit-scrollbar{width:6px;height:6px}
  ::-webkit-scrollbar-track{background:var(--bg)}
  ::-webkit-scrollbar-thumb{background:#5f6368;border-radius:3px}
  ::-webkit-scrollbar-thumb:hover{background:#4b5563}
"""

    filter_keys_json = json.dumps([k for k, _ in col_defs if k is not None])

    html_out = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1.0"/>
<title>QC Summary Report</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/xlsx/0.18.5/xlsx.full.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/jspdf/2.5.1/jspdf.umd.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/jspdf-autotable/3.8.2/jspdf.plugin.autotable.min.js"></script>
<style>{css}</style>
</head>
<body>

<header>
  <div>
    <div class="logo">BASIL <span>//</span> Paired-end QC Report</div>
    <div style="font-size:.85rem;color:var(--muted);margin-top:4px;font-family:var(--font-main)">
      Paired-end Read Quality Control
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

<div class="controls">
  <button class="reset-btn" onclick="resetFilters()">&#10005; Reset All</button>
  <label class="rows-wrap">
    <span>Rows per page</span>
    <select onchange="setRowsPerPage(Number(this.value))">
      <option value="20" selected>20</option>
      <option value="50">50</option>
      <option value="100">100</option>
    </select>
  </label>
</div>
<div class="ms-badge-list" id="globalBadges"></div>

<div class="table-wrap">
  <table>
    <thead>
      <tr class="header-row">
        {header_cells}
      </tr>
      <tr class="filter-row">
        {filter_cells}
      </tr>
    </thead>
    <tbody id="tableBody"></tbody>
  </table>
  <div class="no-data" id="noData" style="display:none">No results match the current filters.</div>
</div>
<div class="pagination" id="pagination"></div>

<script>
const RAW_DATA  = {data_json};
const MAX_COV   = Math.max(...RAW_DATA.map(r => r.estimated_coverage || 0), 1);
const COLS      = {filter_keys_json};

let filtered    = [];
let sortCol     = 'sample_id';
let sortAsc     = true;
let currentPage = 1;
let pageSize    = 20;
const selSets   = {{}};

//  Build dropdowns from data 
function colValToDisplay(col, v) {{
  if (col === 'downsampled') return v === true || v === 'true' ? 'Yes' : 'No';
  if (col === 'qc_status')   return v || 'UNKNOWN';
  return (v === null || v === undefined) ? '' : String(v);
}}

function buildDropdowns() {{
  COLS.forEach(col => {{
    selSets[col] = new Set();
    const vals = [...new Set(RAW_DATA.map(r => colValToDisplay(col, r[col])))].sort((a, b) => {{
      const na = Number(a), nb = Number(b);
      return (!isNaN(na) && !isNaN(nb)) ? na - nb : a.localeCompare(b);
    }});
    const container = document.getElementById('opts-ms-' + col);
    if (!container) return;
    container.innerHTML = vals.map(v =>
      `<label class="ms-opt">
        <input type="checkbox" value="${{v}}"
          onchange="msChange('ms-${{col}}','${{col}}')"/>
        ${{v === '' ? '(empty)' : v}}
      </label>`
    ).join('');
  }});
}}

//  Close dropdowns on outside click 
document.addEventListener('click', e => {{
  if (!e.target.closest('.ms-wrap')) {{
    document.querySelectorAll('.ms-dropdown.open').forEach(d => {{
      d.classList.remove('open');
      const t = document.getElementById(d.id.replace('dd-', 'trigger-'));
      if (t) t.classList.remove('open');
    }});
  }}
}});

function toggleMS(id, col) {{
  const dd = document.getElementById('dd-' + id);
  const tr = document.getElementById('trigger-' + id);
  const isOpen = dd.classList.contains('open');
  document.querySelectorAll('.ms-dropdown.open').forEach(d => {{
    d.classList.remove('open');
    const t = document.getElementById(d.id.replace('dd-', 'trigger-'));
    if (t) t.classList.remove('open');
  }});
  if (!isOpen) {{ dd.classList.add('open'); tr.classList.add('open'); }}
}}

function searchMS(id, q) {{
  q = q.toLowerCase();
  document.querySelectorAll('#opts-' + id + ' .ms-opt').forEach(opt => {{
    opt.style.display = opt.querySelector('input').value.toLowerCase().includes(q) ? '' : 'none';
  }});
}}

function msChange(id, col) {{
  selSets[col] = new Set(
    [...document.querySelectorAll('#opts-' + id + ' input:checked')].map(c => c.value)
  );
  updateLabel(id, col);
  renderBadges();
  applyFilters();
}}

function msAll(id, col) {{
  document.querySelectorAll('#opts-' + id + ' .ms-opt').forEach(opt => {{
    if (opt.style.display !== 'none') opt.querySelector('input').checked = true;
  }});
  msChange(id, col);
}}

function msClear(id, col) {{
  document.querySelectorAll('#opts-' + id + ' input').forEach(c => c.checked = false);
  msChange(id, col);
}}

function updateLabel(id, col) {{
  const lbl = document.getElementById('label-' + id);
  if (!lbl) return;
  const sz = selSets[col] ? selSets[col].size : 0;
  lbl.textContent = sz === 0 ? 'All' : sz === 1 ? [...selSets[col]][0] : sz + ' selected';
}}

//  Badge list 
function renderBadges() {{
  const container = document.getElementById('globalBadges');
  container.innerHTML = '';
  COLS.forEach(col => {{
    if (!selSets[col] || selSets[col].size === 0) return;
    [...selSets[col]].forEach(v => {{
      const div = document.createElement('div');
      div.className = 'ms-badge';
      const txt = document.createTextNode(col + ': ' + v + ' ');
      const x   = document.createElement('span');
      x.className = 'x';
      x.innerHTML = '&#10005;';
      x.addEventListener('click', (function(c, val) {{
        return function() {{ removeBadge(c, val); }};
      }})(col, v));
      div.appendChild(txt);
      div.appendChild(x);
      container.appendChild(div);
    }});
  }});
}}

function removeBadge(col, val) {{
  if (selSets[col]) selSets[col].delete(val);
  const cb = document.querySelector('#opts-ms-' + col + ' input[value="' + val + '"]');
  if (cb) cb.checked = false;
  updateLabel('ms-' + col, col);
  renderBadges();
  applyFilters();
}}

//  Stats 
function renderStats(data) {{
  const total       = data.length;
  const passed      = data.filter(r => r.qc_status === 'PASS').length;
  const failed      = data.filter(r => r.qc_status === 'FAIL').length;
  const downsampled = data.filter(r => r.downsampled).length;
  const warned      = data.filter(r => r.warnings && r.warnings.length > 0).length;
  const avgCov      = data.length > 0
    ? (data.reduce((s, r) => s + (r.estimated_coverage || 0), 0) / data.length).toFixed(1) : 0;

  document.getElementById('statsBar').innerHTML =
    `<div class="stat-card"><div class="stat-label">Total Samples</div><div class="stat-value">${{total}}</div></div>` +
    `<div class="stat-card"><div class="stat-label">QC Pass</div><div class="stat-value">${{passed}}</div></div>` +
    `<div class="stat-card"><div class="stat-label">QC Fail</div><div class="stat-value fail">${{failed}}</div></div>` +
    `<div class="stat-card"><div class="stat-label">Downsampled</div><div class="stat-value">${{downsampled}}</div></div>` +
    `<div class="stat-card"><div class="stat-label">Avg Coverage</div><div class="stat-value">${{avgCov}}X</div></div>` +
    `<div class="stat-card"><div class="stat-label">With Warnings</div><div class="stat-value">${{warned}}</div></div>`;
}}

//  QC status badge helper 
function qcBadge(status) {{
  if (status === 'PASS')    return '<span class="badge badge-pass">&#10003; PASS</span>';
  if (status === 'FAIL')    return '<span class="badge badge-fail">&#10007; FAIL</span>';
  return '<span class="badge badge-unknown">&#8212; UNKNOWN</span>';
}}

//  Table render 
function fmtNum(n, decimals, suffix) {{
  decimals = decimals || 0; suffix = suffix || '';
  if (n === null || n === undefined) return '<span style="color:var(--muted)">&#8212;</span>';
  return Number(n).toLocaleString(undefined,
    {{minimumFractionDigits:decimals,maximumFractionDigits:decimals}}) + suffix;
}}

function renderTable(data, page) {{
  const pages  = Math.max(1, Math.ceil(data.length / pageSize));
  const safePage = Math.min(Math.max(page, 1), pages);
  currentPage = safePage;
  const start  = (safePage - 1) * pageSize;
  const slice  = data.slice(start, start + pageSize);
  const tbody  = document.getElementById('tableBody');
  const noData = document.getElementById('noData');
  if (data.length === 0) {{
    tbody.innerHTML = ''; noData.style.display = 'block';
    document.getElementById('pagination').innerHTML = ''; return;
  }}
  noData.style.display = 'none';
  tbody.innerHTML = slice.map(r => {{
    const covPct = Math.min(((r.estimated_coverage || 0) / MAX_COV) * 100, 100);
    const covBar =
      '<div class="cov-bar-wrap">' +
      '<div class="cov-track"><div class="cov-fill" style="width:' + covPct.toFixed(1) + '%"></div></div>' +
      '<span class="cov-num">' + fmtNum(r.estimated_coverage, 2) + 'X</span></div>';
    const dsBadge  = r.downsampled
      ? '<span class="badge badge-yes">Yes</span>'
      : '<span class="badge badge-no">No</span>';
    const dsFactor = r.downsampled
      ? '<span style="font-size:.8rem">' + fmtNum(r.downsample_factor, 4) + '</span>'
      : '<span style="color:var(--muted)">&#8212;</span>';
    const tgtCov = r.downsampled
      ? '<span style="font-size:.8rem">' + fmtNum(r.target_coverage, 0) + 'X</span>'
      : '<span style="color:var(--muted)">&#8212;</span>';

    //  NEW: fail reason cell 
    const failCell = (r.qc_status === 'FAIL' && r.qc_fail_reason)
      ? '<span class="fail-reason">&#9888; ' + r.qc_fail_reason + '</span>'
      : '<span style="color:var(--muted)">&#8212;</span>';

    const warnCell = r.warnings && r.warnings.length > 0
      ? '<div class="warn-list">' +
        r.warnings.map(w =>
          '<span class="warn-badge">&#9888; ' +
          w.replace(/^[\u26a0\ufe0f\u2192\u0020\t]+/, '') +
          '</span>'
        ).join('') + '</div>'
      : '<span style="color:var(--muted)">&#8212;</span>';

    const rowClass = r.qc_status === 'FAIL' ? ' class="row-fail"' : '';

    return '<tr' + rowClass + '>' +
      '<td class="sample-cell">' + r.sample_id + '</td>' +
      '<td>' + qcBadge(r.qc_status) + '</td>' +
      '<td>' + fmtNum(r.total_reads) + '</td>' +
      '<td>' + fmtNum(r.total_bases) + '</td>' +
      '<td>' + fmtNum(r.mean_read_length, 2) + ' bp</td>' +
      '<td>' + fmtNum(r.genome_size) + ' bp</td>' +
      '<td style="font-size:.8rem;color:var(--muted)">' + (r.genome_size_method || '&#8212;') + '</td>' +
      '<td>' + covBar + '</td>' +
      '<td>' + dsBadge + '</td>' +
      '<td>' + dsFactor + '</td>' +
      '<td>' + tgtCov + '</td>' +
      '<td style="white-space:normal;min-width:200px">' + failCell + '</td>' +
      '<td style="white-space:normal;min-width:200px">' + warnCell + '</td>' +
      '</tr>';
  }}).join('');
  renderPagination(data.length, safePage);
}}

function renderPagination(total, page) {{
  const pages = Math.ceil(total / pageSize);
  const el = document.getElementById('pagination');
  if (pages <= 1) {{ el.innerHTML = ''; return; }}
  let btns = '<span class="page-info">Showing ' +
    Math.min((page-1)*pageSize+1, total) + '&#8211;' + Math.min(page*pageSize, total) +
    ' of ' + total + '</span>';
  btns += '<button class="page-btn" onclick="goPage(' + (page-1) + ')"' + (page===1?' disabled':'') + '>&#8249; Prev</button>';
  const range = [];
  for (let i = Math.max(1, page-2); i <= Math.min(pages, page+2); i++) range.push(i);
  if (range[0] > 1) btns += '<button class="page-btn" onclick="goPage(1)">1</button>';
  range.forEach(i => {{
    btns += '<button class="page-btn' + (i===page?' active':'') + '" onclick="goPage(' + i + ')">' + i + '</button>';
  }});
  if (range[range.length-1] < pages)
    btns += '<button class="page-btn" onclick="goPage(' + pages + ')">' + pages + '</button>';
  btns += '<button class="page-btn" onclick="goPage(' + (page+1) + ')"' + (page===pages?' disabled':'') + '>Next &#8250;</button>';
  el.innerHTML = btns;
}}

function goPage(p) {{
  currentPage = p;
  renderTable(filtered, currentPage);
  window.scrollTo({{top: 0, behavior: 'smooth'}});
}}

function setRowsPerPage(value) {{
  pageSize = value;
  currentPage = 1;
  renderTable(filtered, currentPage);
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

//  Filters 
function applyFilters() {{
  filtered = RAW_DATA.filter(r => {{
    for (const col of COLS) {{
      const sel = selSets[col];
      if (!sel || sel.size === 0) continue;
      const vs = colValToDisplay(col, r[col]);
      if (!sel.has(vs)) return false;
    }}
    return true;
  }});
  sortData();
  currentPage = 1;
  renderStats(filtered);
  renderTable(filtered, currentPage);
}}

function resetFilters() {{
  COLS.forEach(col => {{
    selSets[col] = new Set();
    const id = 'ms-' + col;
    document.querySelectorAll('#opts-' + id + ' input').forEach(c => c.checked = false);
    updateLabel(id, col);
  }});
  document.getElementById('globalBadges').innerHTML = '';
  applyFilters();
}}

//  Sort 
function sortData() {{
  filtered.sort((a, b) => {{
    let va = a[sortCol], vb = b[sortCol];
    if (va === null || va === undefined) va = sortAsc ?  Infinity : -Infinity;
    if (vb === null || vb === undefined) vb = sortAsc ?  Infinity : -Infinity;
    if (typeof va === 'boolean') return sortAsc ? (va===vb?0:va?1:-1) : (va===vb?0:va?-1:1);
    if (typeof va === 'number' && typeof vb === 'number') return sortAsc ? va-vb : vb-va;
    return sortAsc ? String(va).localeCompare(String(vb)) : String(vb).localeCompare(String(va));
  }});
}}

function sortBy(col) {{
  if (sortCol === col) sortAsc = !sortAsc; else {{ sortCol = col; sortAsc = true; }}
  document.querySelectorAll('th[data-col]').forEach(th => {{
    const active = th.dataset.col === col;
    th.classList.toggle('sorted', active);
    th.querySelector('.si').innerHTML = active ? (sortAsc ? '&#8593;' : '&#8595;') : '&#8597;';
  }});
  sortData();
  renderTable(filtered, currentPage);
}}

//  Export 
const EXPORT_HEADERS = ['Sample ID','QC Status','Fail Reason',
  'Total Reads','Total Bases','Mean Read Length (bp)',
  'Genome Size (bp)','Method','Coverage (X)','Downsampled',
  'Downsample Factor','Target Coverage (X)','Warnings'];

function getExportRows() {{
  return filtered.map(r => [
    r.sample_id,
    r.qc_status   ?? '',
    r.qc_fail_reason ?? '',
    r.total_reads ?? '', r.total_bases ?? '', r.mean_read_length ?? '',
    r.genome_size ?? '', r.genome_size_method ?? '', r.estimated_coverage ?? '',
    r.downsampled ? 'Yes' : 'No', r.downsample_factor ?? '', r.target_coverage ?? '',
    (r.warnings || []).join(' | '),
  ]);
}}

function exportCSV() {{
  const rows = [EXPORT_HEADERS, ...getExportRows()];
  downloadBlob(rows.map(r => r.map(v => '"' + String(v).replace(/"/g, '""') + '"').join(',')).join('\\n'),
    'read_qc_summary_report.csv', 'text/csv');
}}
function exportTSV() {{
  downloadBlob([EXPORT_HEADERS, ...getExportRows()].map(r => r.join('\\t')).join('\\n'),
    'read_qc_summary_report.tsv', 'text/tab-separated-values');
}}
function exportExcel() {{
  const wb = XLSX.utils.book_new();
  const ws = XLSX.utils.aoa_to_sheet([EXPORT_HEADERS, ...getExportRows()]);
  ws['!cols'] = EXPORT_HEADERS.map(() => ({{wch: 22}}));
  XLSX.utils.book_append_sheet(wb, ws, 'QC Summary');
  XLSX.writeFile(wb, 'read_qc_summary_report.xlsx');
}}
function exportPDF() {{
  const {{jsPDF}} = window.jspdf;
  const doc = new jsPDF({{orientation:'landscape', unit:'mm', format:'a4'}});
  doc.setFont('helvetica','bold'); doc.setFontSize(14); doc.setTextColor(10,158,116);
  doc.text('QC Summary Report', 14, 16);
  doc.setFont('helvetica','normal'); doc.setFontSize(8); doc.setTextColor(108,117,125);
  doc.text('Generated: ' + new Date().toLocaleString() + ' | Samples: ' + filtered.length, 14, 22);
  doc.autoTable({{startY:27, head:[EXPORT_HEADERS], body:getExportRows(),
    styles:{{font:'helvetica', fontSize:6, cellPadding:2, overflow:'linebreak'}},
    headStyles:{{fillColor:[10,158,116], textColor:[255,255,255], fontStyle:'bold', fontSize:7}},
    alternateRowStyles:{{fillColor:[245,247,250]}}, theme:'grid'}});
  doc.save('read_qc_summary_report.pdf');
}}
function downloadBlob(content, filename, mime) {{
  const blob = new Blob([content], {{type: mime}});
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob); a.download = filename; a.click();
}}

// Init
buildDropdowns();
initDragScroll();
applyFilters();
</script>
</body>
</html>"""

    with open(output_html, "w") as fh:
        fh.write(html_out)
    print(f"Report written to: {output_html}")


if __name__ == "__main__":
    input_dir   = sys.argv[1]
    output_html = sys.argv[2] if len(sys.argv) > 2 else "read_qc_summary_report.html"
    main(input_dir, output_html)
