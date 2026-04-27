#!/usr/bin/env python3

import os
import sys
import glob
import pandas as pd
import json

def parse_abricate_tsv(tsv_file, db_name):
    sample_genes = {}
    with open(tsv_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            sample_id = parts[0].replace(".fasta", "").replace(".fna", "")
            gene      = parts[5]
            if sample_id not in sample_genes:
                sample_genes[sample_id] = []
            if gene not in sample_genes[sample_id]:
                sample_genes[sample_id].append(gene)

    rows = []
    for sample_id, genes in sample_genes.items():
        rows.append({
            "sample_id":      sample_id,
            "database":       db_name,
            "genes_detected": ", ".join(sorted(genes)),
        })
    return rows

def main(input_dir, output_html):
    tsv_files = glob.glob(os.path.join(input_dir, "*.tsv"))
    if not tsv_files:
        print(f"No TSV files found in {input_dir}")
        sys.exit(1)

    all_rows = []
    for tsv in sorted(tsv_files):
        db_name = os.path.splitext(os.path.basename(tsv))[0]
        rows = parse_abricate_tsv(tsv, db_name)
        all_rows.extend(rows)
        print(f"  Parsed {len(rows):>4} samples from {db_name}")

    if not all_rows:
        print("No hits found across all databases.")

    df        = pd.DataFrame(all_rows)
    data_json = df.to_json(orient="records")
    databases = sorted(df["database"].unique().tolist()) if not df.empty else []
    samples   = sorted(df["sample_id"].unique().tolist()) if not df.empty else []

    sample_options = "".join(
        f'<label class="multi-option"><input type="checkbox" value="{s}" onchange="onMultiChange(\'sample\')"/>{s}</label>'
        for s in samples
    )
    db_options = "".join(
        f'<option value="{d}">{d}</option>' for d in databases
    )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1.0"/>
<title>Abricate AMR Report</title>
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
    --font-main:Calibri,'Calibri',Arial,sans-serif;
  }}
  *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ background: var(--bg); color: var(--text); font-family: var(--font-main); min-height: 100vh; }}

  /* HEADER */
  header {{
    background: var(--surface); border-bottom: 2px solid var(--accent); box-shadow: 0 2px 8px rgba(0,0,0,0.06);
    padding: 24px 40px; display: flex; align-items: center;
    justify-content: space-between; gap: 16px; flex-wrap: wrap;
  }}
  .logo {{ font-family: var(--font-main); font-size: 1.4rem; font-weight: 600; color: var(--accent); letter-spacing: 2px; text-transform: uppercase; }}
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
    overflow: hidden; min-width: 160px; z-index: 100; box-shadow: 0 8px 24px rgba(0,0,0,0.4);
  }}
  .export-wrap:hover .export-menu, .export-wrap:focus-within .export-menu {{ display: block; }}
  .export-menu button {{
    display: flex; align-items: center; gap: 10px; width: 100%;
    background: transparent; border: none; color: var(--text); padding: 11px 18px;
    font-family: var(--font-main); font-size: 0.72rem; letter-spacing: 0.5px;
    cursor: pointer; transition: background 0.15s; text-align: left;
  }}
  .export-menu button:hover {{ background: var(--surface); color: var(--accent); }}
  .export-menu button .icon {{ font-size: 1rem; width: 18px; text-align: center; }}

  /* STATS */
  .stats-bar {{ display: flex; gap: 16px; padding: 24px 40px; flex-wrap: wrap; }}
  .stat-card {{
    background: var(--surface); border: 1px solid var(--border);
    border-radius: 10px; padding: 18px 28px; flex: 1; min-width: 160px;
  }}
  .stat-label {{ font-size: 0.7rem; font-family: var(--font-main); color: var(--muted); text-transform: uppercase; letter-spacing: 1.5px; margin-bottom: 8px; }}
  .stat-value {{ font-size: 2rem; font-weight: 600; color: var(--accent); font-family: var(--font-main); }}

  /* CONTROLS */
  .controls {{ padding: 0 40px 8px; display: flex; gap: 12px; flex-wrap: wrap; align-items: center; }}
  .filter-group {{ display: flex; flex-direction: column; gap: 4px; }}
  .filter-label {{ font-size: 0.65rem; font-family: var(--font-main); color: var(--muted); text-transform: uppercase; letter-spacing: 1px; }}

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
    padding: 8px 14px; border-radius: 6px; font-family: var(--font-main); font-size: 0.875rem;
    cursor: pointer; display: flex; justify-content: space-between; align-items: center;
    gap: 8px; user-select: none; transition: border-color 0.2s;
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
    font-family: var(--font-main); font-size: 0.65rem; letter-spacing: 0.5px;
    cursor: pointer; padding: 2px 4px;
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
  .multi-tag-list {{ display: flex; flex-wrap: wrap; gap: 4px; padding: 4px 40px 12px; min-height: 0; }}
  .multi-tag {{
    background: rgba(0,217,163,0.12); color: var(--accent);
    border: 1px solid rgba(0,217,163,0.3); border-radius: 20px;
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
  .table-wrap {{ padding: 0 40px 40px; overflow-x: auto; cursor: grab; scrollbar-gutter: stable; touch-action: pan-y; scrollbar-color: #5f6368 var(--bg); scrollbar-width: thin; }}
  .table-wrap.dragging {{ cursor: grabbing; user-select: none; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 0.875rem; }}
  thead tr.header-row {{ background: var(--accent); }}
  thead tr.filter-row {{ background: var(--surface2); border-bottom: 2px solid var(--border); }}
  th {{
    padding: 12px 14px; text-align: left; font-family: var(--font-main);
    font-size: 0.9rem; letter-spacing: 1px; text-transform: uppercase;
    color: #fff; cursor: pointer; user-select: none; white-space: nowrap;
    border-right: 1px solid rgba(255,255,255,0.15);
  }}
  th:last-child {{ border-right: none; }}
  th:hover {{ background: rgba(255,255,255,0.1); color: #fff; }}
  th .sort-icon {{ margin-left: 4px; opacity: 0.4; }}
  th.sorted .sort-icon {{ opacity: 1; color: var(--accent); }}
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
  .badge {{ display: inline-block; padding: 3px 10px; border-radius: 20px; font-size: 1rem; font-family: var(--font-main); font-weight: 700; letter-spacing: 0.5px; margin: 2px; }}
  .badge-db   {{ background: rgba(10,158,116,0.12); color: #000; border: 1px solid rgba(10,158,116,0.25); }}
  .badge-gene {{ background: rgba(230,126,0,0.1); color: #000; border: 1px solid rgba(230,126,0,0.25); }}
  .no-data {{ text-align: center; padding: 60px; color: var(--muted); font-family: var(--font-main); font-size: 1rem; letter-spacing: 1px; }}

  /* PAGINATION */
  .pagination {{ display: flex; align-items: center; justify-content: flex-end; gap: 8px; padding: 0 40px 40px; flex-wrap: wrap; }}
  .page-info {{ font-size: 0.75rem; font-family: var(--font-main); color: var(--muted); margin-right: 12px; }}
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

<header>
  <div>
    <div class="logo">BASIL <span>//</span> ABRICATE Report</div>
    <div style="font-size:0.78rem;color:var(--muted);margin-top:4px;font-family:var(--font-main);">Cummulative report across all databases</div>
  </div>
  <div class="export-wrap" tabindex="0">
    <button class="export-btn">⬇ Export ▾</button>
    <div class="export-menu">
      <button onclick="exportCSV()"><span class="icon">📄</span> CSV</button>
      <button onclick="exportTSV()"><span class="icon">📋</span> TSV</button>
      <button onclick="exportExcel()"><span class="icon">📊</span> Excel (.xlsx)</button>
      <button onclick="exportPDF()"><span class="icon">📕</span> PDF</button>
    </div>
  </div>
</header>

<div class="stats-bar" id="statsBar"></div>

<div class="controls">
  <div class="filter-group">
    <div class="filter-label">Samples (multi-select)</div>
    <div class="multi-wrap" id="sampleMultiWrap">
      <div class="multi-trigger" id="sampleTrigger" onclick="toggleMulti('sample')">
        <span id="sampleTriggerText">All Samples</span>
        <span class="arrow">▼</span>
      </div>
      <div class="multi-dropdown" id="sampleDropdown">
        <div class="multi-search">
          <input type="text" id="sampleSearch" placeholder="Search samples..." oninput="filterMultiOptions('sample', this.value)"/>
        </div>
        <div class="multi-actions">
          <button onclick="selectAllMulti('sample')">Select All</button>
          <button onclick="clearMulti('sample')">Clear</button>
        </div>
        <div id="sampleOptions">{sample_options}</div>
      </div>
    </div>
  </div>

  <div class="filter-group">
    <div class="filter-label">Database</div>
    <select id="filterDB" onchange="applyFilters()">
      <option value="">All Databases</option>
      {db_options}
    </select>
  </div>

  <div class="filter-group">
    <div class="filter-label">Rows per page</div>
    <select onchange="setRowsPerPage(Number(this.value))">
      <option value="20" selected>20</option>
      <option value="50">50</option>
      <option value="100">100</option>
    </select>
  </div>

  <button class="reset-btn" onclick="resetFilters()">✕ Reset All</button>
</div>

<div class="multi-tag-list" id="sampleTags"></div>

<div class="table-wrap">
  <table id="mainTable">
    <thead>
      <tr class="header-row">
        <th onclick="sortBy('sample_id')"      data-col="sample_id">Sample ID <span class="sort-icon">⇅</span></th>
        <th onclick="sortBy('database')"       data-col="database">Database <span class="sort-icon">⇅</span></th>
        <th onclick="sortBy('genes_detected')" data-col="genes_detected">Genes Detected <span class="sort-icon">⇅</span></th>
      </tr>
      <tr class="filter-row">
        <td class="filter-cell"><input class="col-filter" id="colSample" placeholder="Filter sample..." oninput="applyFilters()"/></td>
        <td class="filter-cell"><input class="col-filter" id="colDB"     placeholder="Filter database..." oninput="applyFilters()"/></td>
        <td class="filter-cell"><input class="col-filter" id="colGene"   placeholder="Filter gene..." oninput="applyFilters()"/></td>
      </tr>
    </thead>
    <tbody id="tableBody"></tbody>
  </table>
  <div class="no-data" id="noData" style="display:none">No results match the current filters.</div>
</div>

<div class="pagination" id="pagination"></div>

<script>
const RAW_DATA = {data_json};
let filtered        = [];
let sortCol         = 'sample_id';
let sortAsc         = true;
let currentPage     = 1;
let pageSize        = 20;
let selectedSamples = new Set();

// ── MULTI-SELECT ─────────────────────────────────────────────
function toggleMulti(key) {{
  const dropdown = document.getElementById(key + 'Dropdown');
  const trigger  = document.getElementById(key + 'Trigger');
  const isOpen   = dropdown.classList.contains('open');
  document.querySelectorAll('.multi-dropdown').forEach(d => d.classList.remove('open'));
  document.querySelectorAll('.multi-trigger').forEach(t => t.classList.remove('open'));
  if (!isOpen) {{
    dropdown.classList.add('open');
    trigger.classList.add('open');
    document.getElementById(key + 'Search').focus();
  }}
}}

document.addEventListener('click', e => {{
  if (!e.target.closest('.multi-wrap')) {{
    document.querySelectorAll('.multi-dropdown').forEach(d => d.classList.remove('open'));
    document.querySelectorAll('.multi-trigger').forEach(t => t.classList.remove('open'));
  }}
}});

function filterMultiOptions(key, query) {{
  const q = query.toLowerCase();
  document.querySelectorAll(`#${{key}}Options .multi-option`).forEach(opt => {{
    opt.style.display = opt.querySelector('input').value.toLowerCase().includes(q) ? '' : 'none';
  }});
}}

function onMultiChange(key) {{
  if (key === 'sample') {{
    selectedSamples = new Set(
      [...document.querySelectorAll('#sampleOptions input:checked')].map(c => c.value)
    );
    updateSampleTrigger();
    renderSampleTags();
  }}
  applyFilters();
}}

function selectAllMulti(key) {{
  document.querySelectorAll(`#${{key}}Options .multi-option`).forEach(opt => {{
    if (opt.style.display !== 'none') opt.querySelector('input').checked = true;
  }});
  onMultiChange(key);
}}

function clearMulti(key) {{
  document.querySelectorAll(`#${{key}}Options input[type=checkbox]`).forEach(c => c.checked = false);
  onMultiChange(key);
}}

function updateSampleTrigger() {{
  const el = document.getElementById('sampleTriggerText');
  if      (selectedSamples.size === 0) el.textContent = 'All Samples';
  else if (selectedSamples.size === 1) el.textContent = [...selectedSamples][0];
  else                                  el.textContent = `${{selectedSamples.size}} samples selected`;
}}

function renderSampleTags() {{
  const container = document.getElementById('sampleTags');
  container.innerHTML = '';
  [...selectedSamples].forEach(s => {{
    const div = document.createElement('div');
    div.className = 'multi-tag';
    div.appendChild(document.createTextNode(s + ' '));
    const x = document.createElement('span');
    x.className = 'remove';
    x.textContent = '✕';
    x.addEventListener('click', (function(val) {{
      return function() {{ removeTag(val); }};
    }})(s));
    div.appendChild(x);
    container.appendChild(div);
  }});
}}

function removeTag(sample) {{
  selectedSamples.delete(sample);
  const cb = document.querySelector(`#sampleOptions input[value='${{sample}}']`);
  if (cb) cb.checked = false;
  updateSampleTrigger();
  renderSampleTags();
  applyFilters();
}}

// ── FILTERS ──────────────────────────────────────────────────
function applyFilters() {{
  const globalDB  = document.getElementById('filterDB').value.toLowerCase();
  const colSample = document.getElementById('colSample').value.toLowerCase();
  const colDB     = document.getElementById('colDB').value.toLowerCase();
  const colGene   = document.getElementById('colGene').value.toLowerCase();

  filtered = RAW_DATA.filter(r => {{
    if (selectedSamples.size > 0 && !selectedSamples.has(r.sample_id)) return false;
    if (globalDB  && !r.database.toLowerCase().includes(globalDB))         return false;
    if (colSample && !r.sample_id.toLowerCase().includes(colSample))       return false;
    if (colDB     && !r.database.toLowerCase().includes(colDB))            return false;
    if (colGene   && !r.genes_detected.toLowerCase().includes(colGene))    return false;
    return true;
  }});

  sortData();
  currentPage = 1;
  renderStats(filtered);
  renderTable(filtered, currentPage);
}}

function resetFilters() {{
  clearMulti('sample');
  document.getElementById('filterDB').value  = '';
  document.getElementById('colSample').value = '';
  document.getElementById('colDB').value     = '';
  document.getElementById('colGene').value   = '';
  document.getElementById('sampleSearch').value = '';
  filterMultiOptions('sample', '');
  applyFilters();
}}

// ── STATS ─────────────────────────────────────────────────────
function renderStats(data) {{
  const samples = new Set(data.map(r => r.sample_id)).size;
  const dbs     = new Set(data.map(r => r.database)).size;
  const genes   = new Set(data.flatMap(r => (r.genes_detected || '').split(', ').filter(Boolean))).size;
  document.getElementById('statsBar').innerHTML = `
    <div class="stat-card"><div class="stat-label">Total Rows</div><div class="stat-value">${{data.length.toLocaleString()}}</div></div>
    <div class="stat-card"><div class="stat-label">Samples</div><div class="stat-value">${{samples}}</div></div>
    <div class="stat-card"><div class="stat-label">Databases</div><div class="stat-value">${{dbs}}</div></div>
    <div class="stat-card"><div class="stat-label">Unique Genes</div><div class="stat-value">${{genes}}</div></div>
  `;
}}

// ── TABLE ─────────────────────────────────────────────────────
function renderTable(data, page) {{
  const pages  = Math.max(1, Math.ceil(data.length / pageSize));
  const safePage = Math.min(Math.max(page, 1), pages);
  currentPage = safePage;
  const start  = (safePage - 1) * pageSize;
  const slice  = data.slice(start, start + pageSize);
  const tbody  = document.getElementById('tableBody');
  const noData = document.getElementById('noData');

  if (data.length === 0) {{
    tbody.innerHTML = '';
    noData.style.display = 'block';
    document.getElementById('pagination').innerHTML = '';
    return;
  }}
  noData.style.display = 'none';

  tbody.innerHTML = slice.map(r => {{
    const genes      = (r.genes_detected || '').split(', ').filter(Boolean);
    const geneBadges = genes.map(g => `<span class="badge badge-gene">${{g}}</span>`).join('');
    return `<tr>
      <td style="font-family:var(--font-main);font-size:1rem;font-weight: 700;">${{r.sample_id}}</td>
      <td><span class="badge badge-db">${{r.database}}</span></td>
      <td style="line-height:2.2">${{geneBadges || '<span style="color:var(--muted)">—</span>'}}</td>
    </tr>`;
  }}).join('');

  renderPagination(data.length, safePage);
}}

function renderPagination(total, page) {{
  const pages = Math.ceil(total / pageSize);
  const el    = document.getElementById('pagination');
  if (pages <= 1) {{ el.innerHTML = ''; return; }}

  let btns = `<span class="page-info">Showing ${{Math.min((page-1)*pageSize+1,total)}}–${{Math.min(page*pageSize,total)}} of ${{total}}</span>`;
  btns += `<button class="page-btn" onclick="goPage(${{page-1}})" ${{page===1?'disabled':''}}>‹ Prev</button>`;
  const range = [];
  for (let i = Math.max(1,page-2); i <= Math.min(pages,page+2); i++) range.push(i);
  if (range[0] > 1) btns += `<button class="page-btn" onclick="goPage(1)">1</button>${{range[0]>2?'<span style="color:var(--muted);padding:0 4px">…</span>':''}}`;
  range.forEach(i => btns += `<button class="page-btn ${{i===page?'active':''}}" onclick="goPage(${{i}})">${{i}}</button>`);
  if (range[range.length-1] < pages) btns += `${{range[range.length-1]<pages-1?'<span style="color:var(--muted);padding:0 4px">…</span>':''}}<button class="page-btn" onclick="goPage(${{pages}})">${{pages}}</button>`;
  btns += `<button class="page-btn" onclick="goPage(${{page+1}})" ${{page===pages?'disabled':''}}>Next ›</button>`;
  el.innerHTML = btns;
}}

function goPage(p) {{
  currentPage = p;
  renderTable(filtered, currentPage);
  window.scrollTo({{top:0, behavior:'smooth'}});
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
      if (e.target.closest('.multi-dropdown, input, select, button, textarea, a, label')) return;

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

// ── SORT ──────────────────────────────────────────────────────
function sortData() {{
  filtered.sort((a, b) => {{
    let va = a[sortCol] || '', vb = b[sortCol] || '';
    const na = parseFloat(va), nb = parseFloat(vb);
    if (!isNaN(na) && !isNaN(nb)) return sortAsc ? na-nb : nb-na;
    return sortAsc ? va.localeCompare(vb) : vb.localeCompare(va);
  }});
}}

function sortBy(col) {{
  if (sortCol === col) sortAsc = !sortAsc;
  else {{ sortCol = col; sortAsc = true; }}
  document.querySelectorAll('th[data-col]').forEach(th => {{
    th.classList.toggle('sorted', th.dataset.col === col);
    th.querySelector('.sort-icon').innerHTML = th.dataset.col === col ? (sortAsc ? '&#8593;' : '&#8595;') : '&#8597;';
  }});
  sortData();
  renderTable(filtered, currentPage);
}}

// ── EXPORT ────────────────────────────────────────────────────
function exportCSV() {{
  const h = ['sample_id','database','genes_detected'];
  const rows = [h.join(','), ...filtered.map(r => h.map(k => `"${{(r[k]||'').toString().replace(/"/g,'""')}}"`).join(','))];
  downloadBlob(rows.join('\\n'), 'abricate_report.csv', 'text/csv');
}}

function exportTSV() {{
  const h = ['sample_id','database','genes_detected'];
  const rows = [h.join('\\t'), ...filtered.map(r => h.map(k => (r[k]||'').toString()).join('\\t'))];
  downloadBlob(rows.join('\\n'), 'abricate_report.tsv', 'text/tab-separated-values');
}}

function exportExcel() {{
  const wsData = [['Sample ID','Database','Genes Detected'], ...filtered.map(r => [r.sample_id, r.database, r.genes_detected])];
  const wb = XLSX.utils.book_new();
  const ws = XLSX.utils.aoa_to_sheet(wsData);
  ws['!cols'] = [{{wch:30}},{{wch:20}},{{wch:80}}];
  XLSX.utils.book_append_sheet(wb, ws, 'AMR Report');
  XLSX.writeFile(wb, 'abricate_report.xlsx');
}}

function exportPDF() {{
  const {{ jsPDF }} = window.jspdf;
  const doc = new jsPDF({{ orientation:'landscape', unit:'mm', format:'a4' }});
  doc.setFont('helvetica','bold'); doc.setFontSize(14); doc.setTextColor(0,217,163);
  doc.text('Abricate AMR Report', 14, 16);
  doc.setFont('helvetica','normal'); doc.setFontSize(8); doc.setTextColor(139,148,158);
  doc.text(`Generated: ${{new Date().toLocaleString()}} | Rows: ${{filtered.length}}`, 14, 22);
  doc.autoTable({{
    startY: 27,
    head: [['Sample ID','Database','Genes Detected']],
    body: filtered.map(r => [r.sample_id, r.database, r.genes_detected]),
    styles: {{ font:'helvetica', fontSize:7, cellPadding:3, overflow:'linebreak' }},
    headStyles: {{ fillColor:[0,217,163], textColor:[0,0,0], fontStyle:'bold', fontSize:8 }},
    alternateRowStyles: {{ fillColor:[245,247,250] }},
    columnStyles: {{ 0:{{cellWidth:50}}, 1:{{cellWidth:30}}, 2:{{cellWidth:'auto'}} }},
    theme: 'grid',
  }});
  doc.save('abricate_report.pdf');
}}

function downloadBlob(content, filename, mime) {{
  const blob = new Blob([content], {{type:mime}});
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob); a.download = filename; a.click();
}}

// Init
initDragScroll();
applyFilters();
</script>
</body>
</html>"""

    with open(output_html, "w") as f:
        f.write(html)
    print(f"Report written to: {output_html}")

if __name__ == "__main__":
    input_dir   = sys.argv[1]
    output_html = sys.argv[2] if len(sys.argv) > 2 else "abricate_report.html"
    main(input_dir, output_html)
