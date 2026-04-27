#!/usr/bin/env python3

import os
import sys
import json
import pandas as pd

def parse_summary_matches(tsv_file):
    df = pd.read_csv(tsv_file, sep="\t", dtype=str).fillna("")
    return df

def main(input_file, output_html):
    if not os.path.exists(input_file):
        print(f"File not found: {input_file}")
        sys.exit(1)

    df = parse_summary_matches(input_file)
    print(f"  Loaded {len(df)} isolates, {len(df.columns)-1} resistance classes")

    isolate_col  = df.columns[0]
    class_cols   = list(df.columns[1:])
    isolates     = sorted(df[isolate_col].unique().tolist())

    # build flat row list for JS
    rows = []
    for _, row in df.iterrows():
        entry = {"isolate": row[isolate_col]}
        for col in class_cols:
            entry[col] = row[col]
        rows.append(entry)

    data_json   = json.dumps(rows)
    cols_json   = json.dumps(class_cols)
    isolate_options = "".join(
        f'<label class="multi-option"><input type="checkbox" value="{s}" onchange="onMultiChange()"/>{s}</label>'
        for s in isolates
    )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1.0"/>
<title>ABRitaMR Report</title>
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
    --accent4:   #1976d2;
    --text:      #1a1f2e;
    --muted:     #6c757d;
    --font-main:Calibri,'Calibri',Arial,sans-serif;
  }}
  *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ background: var(--bg); color: var(--text); font-family: var(--font-main); min-height: 100vh; }}

  /* HEADER */
  header {{
    background: var(--surface); border-bottom: 2px solid var(--accent);
    padding: 20px 40px; display: flex; align-items: center;
    justify-content: space-between; gap: 16px; flex-wrap: wrap;
    box-shadow: 0 2px 8px rgba(0,0,0,0.06);
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
    background: var(--surface); border: 1px solid var(--border); border-radius: 8px;
    overflow: hidden; min-width: 160px; z-index: 100;
    box-shadow: 0 8px 24px rgba(0,0,0,0.12);
  }}
  .export-wrap:hover .export-menu, .export-wrap:focus-within .export-menu {{ display: block; }}
  .export-menu button {{
    display: flex; align-items: center; gap: 10px; width: 100%;
    background: transparent; border: none; color: var(--text); padding: 11px 18px;
    font-family: var(--font-main); font-size: 0.72rem; cursor: pointer;
    transition: background 0.15s; text-align: left;
  }}
  .export-menu button:hover {{ background: var(--surface2); color: var(--accent); }}

  /* STATS */
  .stats-bar {{ display: flex; gap: 16px; padding: 24px 40px; flex-wrap: wrap; }}
  .stat-card {{
    background: var(--surface); border: 1px solid var(--border);
    border-radius: 10px; padding: 16px 24px; flex: 1; min-width: 140px;
    box-shadow: 0 1px 4px rgba(0,0,0,0.05);
  }}
  .stat-label {{ font-size: 0.7rem; font-family: var(--font-main); color: var(--muted); text-transform: uppercase; letter-spacing: 1.5px; margin-bottom: 6px; }}
  .stat-value {{ font-size: 1.8rem; font-weight: 600; color: var(--accent); font-family: var(--font-main); }}

  /* CONTROLS */
  .controls {{ padding: 0 40px 8px; display: flex; gap: 12px; flex-wrap: wrap; align-items: center; }}
  .filter-group {{ display: flex; flex-direction: column; gap: 4px; }}
  .filter-label {{ font-size: 0.65rem; font-family: var(--font-main); color: var(--muted); text-transform: uppercase; letter-spacing: 1px; }}
  input[type=text], select {{
    background: var(--surface); border: 1px solid var(--border); color: var(--text);
    padding: 8px 14px; border-radius: 6px; font-family: var(--font-main);
    font-size: 0.875rem; outline: none; transition: border-color 0.2s; min-width: 180px;
  }}
  input[type=text]:focus, select:focus {{ border-color: var(--accent); }}

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
    background: var(--surface); border: 1px solid var(--border); border-radius: 8px;
    z-index: 200; min-width: 100%; max-height: 280px; overflow-y: auto;
    box-shadow: 0 8px 24px rgba(0,0,0,0.12);
  }}
  .multi-dropdown.open {{ display: block; }}
  .multi-search {{
    padding: 8px 12px; border-bottom: 1px solid var(--border);
    position: sticky; top: 0; background: var(--surface); z-index: 1;
  }}
  .multi-search input {{
    width: 100%; background: var(--surface2); border: 1px solid var(--border);
    color: var(--text); padding: 6px 10px; border-radius: 4px;
    font-size: 0.8rem; outline: none; min-width: unset;
  }}
  .multi-search input:focus {{ border-color: var(--accent); }}
  .multi-actions {{
    display: flex; gap: 8px; padding: 6px 12px;
    border-bottom: 1px solid var(--border);
    position: sticky; top: 45px; background: var(--surface); z-index: 1;
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
  .multi-option:hover {{ background: var(--surface2); }}
  .multi-option input[type=checkbox] {{ accent-color: var(--accent); width: 14px; height: 14px; min-width: unset; cursor: pointer; }}

  /* TAGS */
  .multi-tag-list {{ display: flex; flex-wrap: wrap; gap: 4px; padding: 4px 40px 12px; }}
  .multi-tag {{
    background: rgba(10,158,116,0.1); color: var(--accent);
    border: 1px solid rgba(10,158,116,0.3); border-radius: 20px;
    padding: 3px 10px; font-size: 0.7rem; font-family: var(--font-main);
    display: flex; align-items: center; gap: 6px;
  }}
  .multi-tag .remove {{ cursor: pointer; opacity: 0.6; }}
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
  table {{ width: 100%; border-collapse: collapse; font-size: 0.8rem; white-space: nowrap; }}
  thead tr.header-row {{ background: var(--accent); }}
  thead tr.filter-row {{ background: var(--surface2); border-bottom: 2px solid var(--border); }}
  th {{
    padding: 12px 14px; text-align: left; font-family: var(--font-main);
    font-size: 0.9rem; letter-spacing: 1px; text-transform: uppercase;
    color: #fff; cursor: pointer; user-select: none; white-space: nowrap;
    border-right: 1px solid rgba(255,255,255,0.15);
  }}
  th:last-child {{ border-right: none; }}
  th:hover {{ background: rgba(255,255,255,0.1); }}
  th .sort-icon {{ margin-left: 4px; opacity: 0.6; }}
  th.sorted .sort-icon {{ opacity: 1; }}
  td.filter-cell {{ padding: 5px 6px; }}
  .col-filter {{
    width: 100%; background: var(--surface); border: 1px solid var(--border);
    color: var(--text); padding: 5px 8px; border-radius: 4px;
    font-family: var(--font-main); font-size: 0.78rem; outline: none;
    transition: border-color 0.2s; min-width: unset; white-space: normal;
  }}
  .col-filter:focus {{ border-color: var(--accent); }}
  .col-filter::placeholder {{ color: var(--muted); font-size: 0.68rem; }}
  tbody tr {{ border-bottom: 1px solid var(--border); transition: background 0.12s; }}
  tbody tr:nth-child(even) {{ background: var(--surface2); }}
  tbody tr:hover {{ background: rgba(10,158,116,0.06); }}
  td {{
    padding: 10px 14px; vertical-align: middle;
    border-right: 1px solid var(--border); white-space: normal; max-width: 220px;
  }}
  td:last-child {{ border-right: none; }}
  td.isolate-cell {{ font-family: var(--font-main); font-size: 1rem; font-weight: 700; color: var(--text); white-space: nowrap; }}
  .gene-badge {{
    display: inline-block; padding: 2px 8px; border-radius: 12px; margin: 2px;
    font-size: 1rem; font-family: var(--font-main); font-weight: 600;
    background: rgba(230,126,0,0.1); color: #000000;
    border: 1px solid rgba(230,126,0,0.25);
  }}
  .empty-cell {{ color: var(--border); font-size: 1rem; text-align: center; }}
  .no-data {{ text-align: center; padding: 60px; color: var(--muted); font-family: var(--font-main); font-size: 1rem; }}

  /* PAGINATION */
  .pagination {{ display: flex; align-items: center; justify-content: flex-end; gap: 8px; padding: 0 40px 40px; flex-wrap: wrap; }}
  .page-info {{ font-size: 0.75rem; font-family: var(--font-main); color: var(--muted); margin-right: 12px; }}
  .page-btn {{
    background: var(--surface); border: 1px solid var(--border); color: var(--text);
    padding: 6px 14px; border-radius: 6px; cursor: pointer;
    font-family: var(--font-main); font-size: 0.75rem; transition: all 0.15s;
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
    <div class="logo">BASIL <span>//</span> ABRITAMR Report</div>
    <div style="font-size:0.78rem;color:var(--muted);margin-top:4px;font-family:var(--font-main);">Antimicrobial Resistance Gene Prediction</div>
  </div>
  <div class="export-wrap" tabindex="0">
    <button class="export-btn">⬇ Export ▾</button>
    <div class="export-menu">
      <button onclick="exportCSV()">📄 CSV</button>
      <button onclick="exportTSV()">📋 TSV</button>
      <button onclick="exportExcel()">📊 Excel (.xlsx)</button>
      <button onclick="exportPDF()">📕 PDF</button>
    </div>
  </div>
</header>

<div class="stats-bar" id="statsBar"></div>

<div class="controls">
  <div class="filter-group">
    <div class="filter-label">Isolates (multi-select)</div>
    <div class="multi-wrap">
      <div class="multi-trigger" id="isolateTrigger" onclick="toggleMulti()">
        <span id="isolateTriggerText">All Isolates</span>
        <span class="arrow">▼</span>
      </div>
      <div class="multi-dropdown" id="isolateDropdown">
        <div class="multi-search">
          <input type="text" id="isolateSearch" placeholder="Search isolates..." oninput="filterMultiOptions(this.value)"/>
        </div>
        <div class="multi-actions">
          <button onclick="selectAll()">Select All</button>
          <button onclick="clearAll()">Clear</button>
        </div>
        <div id="isolateOptions">{isolate_options}</div>
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

  <button class="reset-btn" onclick="resetFilters()">✕ Reset All</button>
</div>

<div class="multi-tag-list" id="isolateTags"></div>

<div class="table-wrap">
  <table id="mainTable">
    <thead>
      <tr class="header-row" id="headerRow"></tr>
      <tr class="filter-row" id="filterRow"></tr>
    </thead>
    <tbody id="tableBody"></tbody>
  </table>
  <div class="no-data" id="noData" style="display:none">No results match the current filters.</div>
</div>
<div class="pagination" id="pagination"></div>

<script>
const RAW_DATA  = {data_json};
const COLS      = {cols_json};

let filtered        = [];
let sortCol         = null;
let sortAsc         = true;
let currentPage     = 1;
let pageSize        = 20;
let selectedIsolates = new Set();
let colFilters      = {{}};   // col -> search string

// ── BUILD DYNAMIC HEADERS ─────────────────────────────────────
function buildHeaders() {{
  const hr = document.getElementById('headerRow');
  const fr = document.getElementById('filterRow');

  // Isolate column
  hr.innerHTML = `<th onclick="sortBy('isolate')" data-col="isolate">Isolate <span class="sort-icon">⇅</span></th>`;
  fr.innerHTML = `<td class="filter-cell"><input class="col-filter" placeholder="Filter..." oninput="setColFilter('isolate', this.value)" id="cf_isolate"/></td>`;

  COLS.forEach(col => {{
    const safe = col.replace(/[^a-zA-Z0-9]/g, '_');
    hr.innerHTML += `<th onclick="sortBy('${{col}}')" data-col="${{col}}">${{col}} <span class="sort-icon">⇅</span></th>`;
    fr.innerHTML += `<td class="filter-cell"><input class="col-filter" placeholder="Filter..." oninput="setColFilter('${{col}}', this.value)" id="cf_${{safe}}"/></td>`;
  }});
}}

// ── MULTI-SELECT ──────────────────────────────────────────────
function toggleMulti() {{
  const dd = document.getElementById('isolateDropdown');
  const tr = document.getElementById('isolateTrigger');
  const open = dd.classList.contains('open');
  dd.classList.toggle('open', !open);
  tr.classList.toggle('open', !open);
  if (!open) document.getElementById('isolateSearch').focus();
}}

document.addEventListener('click', e => {{
  if (!e.target.closest('.multi-wrap')) {{
    document.getElementById('isolateDropdown').classList.remove('open');
    document.getElementById('isolateTrigger').classList.remove('open');
  }}
}});

function filterMultiOptions(q) {{
  q = q.toLowerCase();
  document.querySelectorAll('#isolateOptions .multi-option').forEach(opt => {{
    opt.style.display = opt.querySelector('input').value.toLowerCase().includes(q) ? '' : 'none';
  }});
}}

function onMultiChange() {{
  selectedIsolates = new Set(
    [...document.querySelectorAll('#isolateOptions input:checked')].map(c => c.value)
  );
  updateTrigger();
  renderTags();
  applyFilters();
}}

function selectAll() {{
  document.querySelectorAll('#isolateOptions .multi-option').forEach(opt => {{
    if (opt.style.display !== 'none') opt.querySelector('input').checked = true;
  }});
  onMultiChange();
}}

function clearAll() {{
  document.querySelectorAll('#isolateOptions input[type=checkbox]').forEach(c => c.checked = false);
  onMultiChange();
}}

function updateTrigger() {{
  const el = document.getElementById('isolateTriggerText');
  if      (selectedIsolates.size === 0) el.textContent = 'All Isolates';
  else if (selectedIsolates.size === 1) el.textContent = [...selectedIsolates][0];
  else                                   el.textContent = `${{selectedIsolates.size}} isolates selected`;
}}

function renderTags() {{
  const container = document.getElementById('isolateTags');
  container.innerHTML = '';
  [...selectedIsolates].forEach(s => {{
    const div = document.createElement('div');
    div.className = 'multi-tag';
    div.appendChild(document.createTextNode(s + ' '));
    const x = document.createElement('span');
    x.className = 'remove';
    x.textContent = '\u2715';
    x.addEventListener('click', (function(val) {{
      return function() {{ removeTag(val); }};
    }})(s));
    div.appendChild(x);
    container.appendChild(div);
  }});
}}

function removeTag(s) {{
  selectedIsolates.delete(s);
  const cb = document.querySelector(`#isolateOptions input[value='${{s}}']`);
  if (cb) cb.checked = false;
  updateTrigger(); renderTags(); applyFilters();
}}

// ── FILTERS ───────────────────────────────────────────────────
function setColFilter(col, val) {{
  colFilters[col] = val.toLowerCase();
  applyFilters();
}}

function applyFilters() {{
  filtered = RAW_DATA.filter(r => {{
    if (selectedIsolates.size > 0 && !selectedIsolates.has(r.isolate)) return false;
    for (const [col, val] of Object.entries(colFilters)) {{
      if (!val) continue;
      const cell = col === 'isolate' ? r.isolate : (r[col] || '');
      if (!cell.toLowerCase().includes(val)) return false;
    }}
    return true;
  }});

  sortData();
  currentPage = 1;
  renderStats(filtered);
  renderTable(filtered, currentPage);
}}

function resetFilters() {{
  clearAll();
  document.getElementById('isolateSearch').value = '';
  filterMultiOptions('');
  colFilters = {{}};
  document.querySelectorAll('.col-filter').forEach(i => i.value = '');
  applyFilters();
}}

// ── STATS ─────────────────────────────────────────────────────
function renderStats(data) {{
  const isolates = new Set(data.map(r => r.isolate)).size;
  // count isolates that have at least one gene detected
  const withGenes = data.filter(r => COLS.some(c => r[c] && r[c].trim())).length;
  // unique genes across all columns
  const allGenes = new Set();
  data.forEach(r => COLS.forEach(c => {{
    if (r[c]) r[c].split(',').forEach(g => {{ const t = g.trim(); if (t) allGenes.add(t); }});
  }}));
  const classes = COLS.filter(c => data.some(r => r[c] && r[c].trim())).length;

  document.getElementById('statsBar').innerHTML = `
    <div class="stat-card"><div class="stat-label">Isolates</div><div class="stat-value">${{isolates}}</div></div>
    <div class="stat-card"><div class="stat-label">With Resistance</div><div class="stat-value">${{withGenes}}</div></div>
    <div class="stat-card"><div class="stat-label">Unique Genes</div><div class="stat-value">${{allGenes.size}}</div></div>
    <div class="stat-card"><div class="stat-label">Active Classes</div><div class="stat-value">${{classes}}</div></div>
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
    tbody.innerHTML = ''; noData.style.display = 'block';
    document.getElementById('pagination').innerHTML = ''; return;
  }}
  noData.style.display = 'none';

  tbody.innerHTML = slice.map(r => {{
    let cells = `<td class="isolate-cell">${{r.isolate}}</td>`;
    COLS.forEach(col => {{
      const val = r[col] || '';
      if (!val.trim()) {{
        cells += `<td><span class="empty-cell">—</span></td>`;
      }} else {{
        const badges = val.split(',').map(g => `<span class="gene-badge">${{g.trim()}}</span>`).join('');
        cells += `<td style="white-space:normal">${{badges}}</td>`;
      }}
    }});
    return `<tr>${{cells}}</tr>`;
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
  if (range[0]>1) btns+=`<button class="page-btn" onclick="goPage(1)">1</button>${{range[0]>2?'<span style="color:var(--muted);padding:0 4px">…</span>':''}}`;
  range.forEach(i => btns+=`<button class="page-btn ${{i===page?'active':''}}" onclick="goPage(${{i}})">${{i}}</button>`);
  if (range[range.length-1]<pages) btns+=`${{range[range.length-1]<pages-1?'<span style="color:var(--muted);padding:0 4px">…</span>':''}}<button class="page-btn" onclick="goPage(${{pages}})">${{pages}}</button>`;
  btns+=`<button class="page-btn" onclick="goPage(${{page+1}})" ${{page===pages?'disabled':''}}>Next ›</button>`;
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
  if (!sortCol) return;
  filtered.sort((a, b) => {{
    const va = (sortCol === 'isolate' ? a.isolate : (a[sortCol]||''));
    const vb = (sortCol === 'isolate' ? b.isolate : (b[sortCol]||''));
    return sortAsc ? va.localeCompare(vb) : vb.localeCompare(va);
  }});
}}

function sortBy(col) {{
  if (sortCol === col) sortAsc = !sortAsc;
  else {{ sortCol = col; sortAsc = true; }}
  document.querySelectorAll('th[data-col]').forEach(th => {{
    th.classList.toggle('sorted', th.dataset.col === col);
    th.querySelector('.sort-icon').innerHTML = th.dataset.col===col ? (sortAsc?'&#8593;':'&#8595;') : '&#8597;';
  }});
  sortData();
  renderTable(filtered, currentPage);
}}

// ── EXPORT ────────────────────────────────────────────────────
function getHeaders() {{ return ['Isolate', ...COLS]; }}
function getRows()    {{ return filtered.map(r => [r.isolate, ...COLS.map(c => r[c]||'')]); }}

function exportCSV() {{
  const h = getHeaders(), rows = getRows();
  const lines = [h.map(v=>`"${{v.replace(/"/g,'""')}}"`).join(','),
                 ...rows.map(r => r.map(v=>`"${{v.replace(/"/g,'""')}}"`).join(','))];
  downloadBlob(lines.join('\\n'), 'abritamr_report.csv', 'text/csv');
}}

function exportTSV() {{
  const h = getHeaders(), rows = getRows();
  const lines = [h.join('\\t'), ...rows.map(r => r.join('\\t'))];
  downloadBlob(lines.join('\\n'), 'abritamr_report.tsv', 'text/tab-separated-values');
}}

function exportExcel() {{
  const wb = XLSX.utils.book_new();
  const ws = XLSX.utils.aoa_to_sheet([getHeaders(), ...getRows()]);
  ws['!cols'] = [{{wch:35}}, ...COLS.map(() => ({{wch:25}}))];
  XLSX.utils.book_append_sheet(wb, ws, 'AMR Summary');
  XLSX.writeFile(wb, 'abritamr_report.xlsx');
}}

function exportPDF() {{
  const {{ jsPDF }} = window.jspdf;
  const doc = new jsPDF({{ orientation:'landscape', unit:'mm', format:'a4' }});
  doc.setFont('helvetica','bold'); doc.setFontSize(14); doc.setTextColor(10,158,116);
  doc.text('ABRitaMR Report', 14, 16);
  doc.setFont('helvetica','normal'); doc.setFontSize(8); doc.setTextColor(108,117,125);
  doc.text(`Generated: ${{new Date().toLocaleString()}} | Isolates: ${{filtered.length}}`, 14, 22);
  doc.autoTable({{
    startY: 27,
    head: [getHeaders()],
    body: getRows(),
    styles: {{ font:'helvetica', fontSize:6, cellPadding:2, overflow:'linebreak' }},
    headStyles: {{ fillColor:[10,158,116], textColor:[255,255,255], fontStyle:'bold', fontSize:7 }},
    alternateRowStyles: {{ fillColor:[245,247,250] }},
    theme: 'grid',
  }});
  doc.save('abritamr_report.pdf');
}}

function downloadBlob(content, filename, mime) {{
  const blob = new Blob([content], {{type:mime}});
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob); a.download = filename; a.click();
}}

// Init
buildHeaders();
initDragScroll();
applyFilters();
</script>
</body>
</html>"""

    with open(output_html, "w") as f:
        f.write(html)
    print(f"Report written to: {output_html}")

if __name__ == "__main__":
    input_file  = sys.argv[1]
    output_html = sys.argv[2] if len(sys.argv) > 2 else "abritamr_report.html"
    main(input_file, output_html)
