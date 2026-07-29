#!/usr/bin/env python3
"""
generate_igv_session.py

BNDregion BAM 파일들을 patient 단위로 묶어서
  1. IGV session XML 생성 (IGV desktop app용)
  2. IGV.js 웹 뷰어 HTML 생성 (샘플 선택 후 로드)

Usage:
    python generate_igv_session.py \
        --bam_dir /mnt/NAS3/.../BNDregion_bam \
        --genome hg19 \
        --out_dir /mnt/NAS3/.../igv_sessions \
        --base_url "http://서버IP:8889/files"

    # 특정 patient + region
    python generate_igv_session.py \
        --bam_dir ... --genome hg19 --out_dir ... \
        --base_url "http://서버IP:8889/files" \
        --patient ECTRES-ECGI1-0001 \
        --locus chr8:127700000-128000000
"""

import os
import argparse
import json
import xml.etree.ElementTree as ET
from xml.dom import minidom
from collections import defaultdict


def to_url(path, base_url):
    if base_url:
        return base_url.rstrip('/') + path
    return path


# ──────────────────────────────────────────────────────────────────────────────
# BAM 탐색
# ──────────────────────────────────────────────────────────────────────────────

def find_bam_files(bam_dir, patient_filter=None):
    patient_bams = defaultdict(list)
    if not os.path.isdir(bam_dir):
        print(f"ERROR: bam_dir not found: {bam_dir}")
        return patient_bams
    for patient_barcode in sorted(os.listdir(bam_dir)):
        if patient_filter and patient_barcode != patient_filter:
            continue
        patient_dir = os.path.join(bam_dir, patient_barcode)
        if not os.path.isdir(patient_dir):
            continue
        for aliquot in sorted(os.listdir(patient_dir)):
            aliquot_dir = os.path.join(patient_dir, aliquot)
            if not os.path.isdir(aliquot_dir):
                continue
            bam_file = os.path.join(aliquot_dir, f"{aliquot}_ms_BNDregion.bam")
            bai_file = bam_file + ".bai"
            if os.path.exists(bam_file):
                patient_bams[patient_barcode].append((aliquot, bam_file, bai_file))
    return patient_bams


# ──────────────────────────────────────────────────────────────────────────────
# 1. IGV session XML (desktop)
# ──────────────────────────────────────────────────────────────────────────────

def generate_igv_session_xml(patient_barcode, bam_list, genome, locus, out_path):
    session = ET.Element("Session")
    session.set("genome", genome)
    session.set("hasGeneTrack", "true")
    session.set("hasSequenceTrack", "true")
    session.set("version", "8")
    if locus:
        session.set("locus", locus)

    resources = ET.SubElement(session, "Resources")
    for aliquot, bam_path, bai_path in bam_list:
        resource = ET.SubElement(resources, "Resource")
        resource.set("path", bam_path)

    panel = ET.SubElement(session, "Panel")
    panel.set("height", str(max(200, len(bam_list) * 90)))
    panel.set("name", "DataPanel")
    panel.set("width", "1800")

    colors = ["150,200,255","255,180,130","130,220,160","255,220,100",
              "210,150,255","255,130,160","100,220,220","255,200,80"]
    for i, (aliquot, bam_path, bai_path) in enumerate(bam_list):
        track = ET.SubElement(panel, "Track")
        track.set("id", bam_path)
        track.set("name", aliquot)
        track.set("color", colors[i % len(colors)])
        track.set("colorByStrand", "true")
        track.set("displayMode", "SQUISHED")
        track.set("fontSize", "10")
        track.set("visible", "true")

    xml_str = minidom.parseString(
        ET.tostring(session, encoding='unicode')
    ).toprettyxml(indent="    ")
    xml_str = '\n'.join(xml_str.split('\n')[1:])
    with open(out_path, 'w') as f:
        f.write(xml_str)
    print(f"  [XML]  {out_path}")


# ──────────────────────────────────────────────────────────────────────────────
# 2. IGV.js 웹 뷰어 HTML (샘플 선택 후 로드)
# ──────────────────────────────────────────────────────────────────────────────

TRACK_COLORS = [
    "#4dabf7","#ff8c42","#51cf66","#ffd43b",
    "#cc5de8","#ff6b9d","#20c997","#ff9f43",
    "#74c0fc","#ffa94d","#a9e34b","#f783ac",
]


def generate_igv_html(patient_barcode, bam_list, genome, locus, out_path, base_url=None):
    # Build track metadata as JSON for JS
    tracks_data = []
    for i, (aliquot, bam_path, bai_path) in enumerate(bam_list):
        bai = bai_path if os.path.exists(bai_path) else bam_path + ".bai"
        tracks_data.append({
            "name": aliquot,
            "url": to_url(bam_path, base_url),
            "indexURL": to_url(bai, base_url),
            "color": TRACK_COLORS[i % len(TRACK_COLORS)],
        })

    tracks_json = json.dumps(tracks_data, indent=2)
    locus_str = f'"{locus}"' if locus else '"all"'
    n_tracks = len(bam_list)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>IGV — {patient_barcode}</title>
<script src="https://cdn.jsdelivr.net/npm/igv@2.15.8/dist/igv.min.js"></script>
<style>
  @import url('https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;600&family=Syne:wght@700;800&display=swap');
  :root {{
    --bg:#0d1117; --s1:#161b22; --s2:#1c2128;
    --bd:#30363d; --ac:#4dabf7; --tx:#e6edf3; --tm:#8b949e;
  }}
  *{{box-sizing:border-box;margin:0;padding:0;}}
  body{{background:var(--bg);color:var(--tx);font-family:'JetBrains Mono',monospace;min-height:100vh;display:flex;flex-direction:column;}}

  /* ── header ── */
  header{{background:var(--s1);border-bottom:1px solid var(--bd);padding:14px 20px;display:flex;align-items:center;gap:12px;position:sticky;top:0;z-index:200;}}
  .badge{{background:var(--ac);color:#0d1117;font-family:'Syne',sans-serif;font-weight:800;font-size:10px;padding:2px 7px;border-radius:4px;letter-spacing:.06em;text-transform:uppercase;}}
  .htitle{{font-family:'Syne',sans-serif;font-size:16px;font-weight:700;}}
  .hsub{{font-size:11px;color:var(--tm);margin-left:auto;}}

  /* ── selector panel ── */
  #selector{{background:var(--s2);border-bottom:1px solid var(--bd);padding:16px 20px;}}
  .sel-header{{display:flex;align-items:center;justify-content:space-between;margin-bottom:12px;}}
  .sel-title{{font-size:12px;color:var(--tm);letter-spacing:.04em;text-transform:uppercase;}}
  .sel-actions{{display:flex;gap:6px;}}
  .sample-grid{{display:grid;grid-template-columns:repeat(auto-fill,minmax(320px,1fr));gap:6px;margin-bottom:14px;}}
  .sample-item{{display:flex;align-items:center;gap:8px;padding:7px 10px;border-radius:6px;border:1px solid var(--bd);cursor:pointer;transition:border-color .15s;}}
  .sample-item:hover{{border-color:var(--ac);}}
  .sample-item.checked{{border-color:#4dabf750;background:#4dabf708;}}
  .dot{{width:9px;height:9px;border-radius:50%;flex-shrink:0;}}
  .sname{{font-size:11px;flex:1;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;}}
  .cb{{accent-color:var(--ac);width:14px;height:14px;cursor:pointer;flex-shrink:0;}}
  .load-bar{{display:flex;align-items:center;gap:12px;}}
  .btn{{border:none;font-family:'JetBrains Mono',monospace;font-size:12px;font-weight:600;padding:8px 16px;border-radius:6px;cursor:pointer;transition:opacity .2s;}}
  .btn-primary{{background:var(--ac);color:#0d1117;}}
  .btn-primary:hover{{opacity:.85;}}
  .btn-primary:disabled{{opacity:.4;cursor:not-allowed;}}
  .btn-sm{{background:var(--s1);color:var(--tx);border:1px solid var(--bd);font-size:11px;padding:5px 10px;border-radius:5px;cursor:pointer;}}
  .btn-sm:hover{{border-color:var(--ac);color:var(--ac);}}
  #sel-count{{font-size:12px;color:var(--tm);}}

  /* ── controls ── */
  #controls{{background:var(--s1);border-bottom:1px solid var(--bd);padding:10px 20px;display:none;align-items:center;gap:10px;flex-wrap:wrap;}}
  .locus-input{{background:var(--s2);border:1px solid var(--bd);color:var(--tx);font-family:'JetBrains Mono',monospace;font-size:12px;padding:6px 10px;border-radius:6px;width:280px;}}
  .locus-input:focus{{outline:none;border-color:var(--ac);}}
  .btn-ctrl{{background:var(--s2);color:var(--tx);border:1px solid var(--bd);font-size:11px;padding:5px 10px;border-radius:5px;cursor:pointer;font-family:'JetBrains Mono',monospace;}}
  .btn-ctrl:hover{{border-color:var(--ac);color:var(--ac);}}

  /* ── igv ── */
  #igv-wrap{{flex:1;padding:14px 20px;min-height:500px;display:none;}}

  /* ── status ── */
  #status-bar{{background:var(--s1);border-top:1px solid var(--bd);padding:5px 20px;font-size:11px;color:var(--tm);display:flex;gap:16px;}}
  #status-bar span{{color:var(--ac);}}

  /* ── loading overlay ── */
  #loading{{display:none;position:fixed;inset:0;background:#0d111799;z-index:500;align-items:center;justify-content:center;flex-direction:column;gap:12px;}}
  #loading.show{{display:flex;}}
  .spinner{{width:36px;height:36px;border:3px solid var(--bd);border-top-color:var(--ac);border-radius:50%;animation:spin .7s linear infinite;}}
  @keyframes spin{{to{{transform:rotate(360deg);}}}}
  #loading-msg{{font-size:13px;color:var(--tm);}}
</style>
</head>
<body>

<header>
  <div class="badge">IGV</div>
  <div class="htitle">{patient_barcode}</div>
  <div class="hsub">{n_tracks} aliquots · {genome}</div>
</header>

<div id="selector">
  <div class="sel-header">
    <span class="sel-title">샘플 선택</span>
    <div class="sel-actions">
      <button class="btn-sm" onclick="selAll(true)">모두 선택</button>
      <button class="btn-sm" onclick="selAll(false)">모두 해제</button>
    </div>
  </div>
  <div class="sample-grid" id="sample-grid"></div>
  <div class="load-bar">
    <button class="btn btn-primary" id="load-btn" onclick="loadSelected()">▶ 선택 샘플 로드</button>
    <span id="sel-count">0개 선택</span>
    <span style="font-size:11px;color:var(--tm);">* 한 번에 3~5개 권장</span>
  </div>
</div>

<div id="controls">
  <input class="locus-input" id="locus-input" type="text" placeholder="chr8:127700000-128000000" value="{locus or ''}">
  <button class="btn-ctrl" onclick="goLocus()">Go</button>
  <button class="btn-ctrl" onclick="igvBrowser.zoomIn()">＋ Zoom</button>
  <button class="btn-ctrl" onclick="igvBrowser.zoomOut()">－ Zoom</button>
  <button class="btn-ctrl" onclick="igvBrowser.search('all')">View All</button>
  <button class="btn-ctrl" onclick="resetView()" style="margin-left:auto;">↩ 다시 선택</button>
</div>

<div id="igv-wrap"></div>

<div id="status-bar">
  <div>Patient <span>{patient_barcode}</span></div>
  <div>Genome <span>{genome}</span></div>
  <div id="status-loaded">로드됨 <span>—</span></div>
</div>

<div id="loading">
  <div class="spinner"></div>
  <div id="loading-msg">BAM 로딩 중...</div>
</div>

<script>
const ALL_TRACKS = {tracks_json};
const GENOME = "{genome}";
const INIT_LOCUS = {locus_str};
let igvBrowser = null;

const grid = document.getElementById('sample-grid');
ALL_TRACKS.forEach((t, i) => {{
  const item = document.createElement('label');
  item.className = 'sample-item checked';
  item.innerHTML = `
    <input type="checkbox" class="cb" checked onchange="onCheck(this)">
    <span class="dot" style="background:${{t.color}}"></span>
    <span class="sname" title="${{t.name}}">${{t.name}}</span>
  `;
  grid.appendChild(item);
}});

function onCheck(cb) {{
  cb.closest('.sample-item').classList.toggle('checked', cb.checked);
  updateCount();
}}

function selAll(v) {{
  grid.querySelectorAll('input').forEach(cb => {{
    cb.checked = v;
    cb.closest('.sample-item').classList.toggle('checked', v);
  }});
  updateCount();
}}

function updateCount() {{
  const n = grid.querySelectorAll('input:checked').length;
  document.getElementById('sel-count').textContent = n + '개 선택';
  document.getElementById('load-btn').disabled = n === 0;
}}

updateCount();

async function loadSelected() {{
  const selected = [...grid.querySelectorAll('input:checked')].map((cb, _) => {{
    const idx = [...grid.querySelectorAll('input')].indexOf(cb);
    return ALL_TRACKS[idx];
  }});
  if (!selected.length) return;

  showLoading(`${{selected.length}}개 샘플 로드 중...`);

  const tracks = selected.map(t => ({{
    name: t.name,
    url: t.url,
    indexURL: t.indexURL,
    format: 'bam',
    displayMode: 'SQUISHED',
    colorBy: 'strand',
    color: t.color,
    height: 80,
  }}));

  try {{
    if (igvBrowser) {{
      await igvBrowser.loadTrackList(tracks);
    }} else {{
      igvBrowser = await igv.createBrowser(document.getElementById('igv-wrap'), {{
        genome: GENOME,
        locus: INIT_LOCUS,
        tracks: tracks,
      }});
    }}
    document.getElementById('selector').style.display = 'none';
    document.getElementById('controls').style.display = 'flex';
    document.getElementById('igv-wrap').style.display = 'block';
    document.getElementById('status-loaded').innerHTML = `로드됨 <span>${{selected.length}}개</span>`;
  }} catch(e) {{
    alert('로드 실패: ' + e.message);
  }} finally {{
    hideLoading();
  }}
}}

function resetView() {{
  document.getElementById('selector').style.display = 'block';
  document.getElementById('controls').style.display = 'none';
  document.getElementById('igv-wrap').style.display = 'none';
  if (igvBrowser) {{
    igv.removeBrowser(igvBrowser);
    igvBrowser = null;
    document.getElementById('igv-wrap').innerHTML = '';
  }}
}}

function goLocus() {{
  const l = document.getElementById('locus-input').value.trim();
  if (l && igvBrowser) igvBrowser.search(l);
}}
document.getElementById('locus-input')?.addEventListener('keydown', e => {{
  if (e.key === 'Enter') goLocus();
}});

function showLoading(msg) {{
  document.getElementById('loading-msg').textContent = msg;
  document.getElementById('loading').classList.add('show');
}}
function hideLoading() {{
  document.getElementById('loading').classList.remove('show');
}}
</script>
</body>
</html>"""

    with open(out_path, 'w') as f:
        f.write(html)
    print(f"  [HTML] {out_path}")


# ──────────────────────────────────────────────────────────────────────────────
# main
# ──────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Generate IGV session XML + IGV.js HTML per patient")
    parser.add_argument("--bam_dir",  required=True)
    parser.add_argument("--genome",   default="hg19")
    parser.add_argument("--out_dir",  required=True)
    parser.add_argument("--patient",  default=None)
    parser.add_argument("--locus",    default=None)
    parser.add_argument("--base_url", default=None,
                        help="HTTP base URL, e.g. http://10.201.134.195:3275/files")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    print(f"Scanning: {args.bam_dir}")
    if args.base_url:
        print(f"Base URL: {args.base_url}")

    patient_bams = find_bam_files(args.bam_dir, patient_filter=args.patient)
    if not patient_bams:
        print("No BAM files found.")
        return

    print(f"Found {len(patient_bams)} patient(s)\n")
    for patient_barcode, bam_list in sorted(patient_bams.items()):
        print(f"[{patient_barcode}]  {len(bam_list)} aliquot(s)")
        for aliquot, bam_path, _ in bam_list:
            print(f"    {'✔' if os.path.exists(bam_path) else '✘'} {aliquot}")

        xml_path  = os.path.join(args.out_dir, f"{patient_barcode}_igv_session.xml")
        html_path = os.path.join(args.out_dir, f"{patient_barcode}_igv_viewer.html")
        generate_igv_session_xml(patient_barcode, bam_list, args.genome, args.locus, xml_path)
        generate_igv_html(patient_barcode, bam_list, args.genome, args.locus, html_path, base_url=args.base_url)
        print()

    print(f"Done! → {args.out_dir}")


if __name__ == "__main__":
    main()