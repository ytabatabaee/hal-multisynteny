"""Deterministic, dependency-free visualization of stored run checkpoints."""

from __future__ import annotations

import csv
import html
import json
import statistics
from collections import Counter, defaultdict
from pathlib import Path

CLASSIFICATIONS = (
    "shared_consistent", "left_only", "right_only", "orientation_conflict",
    "order_conflict", "duplication_conflict", "ambiguous", "unaligned", "complex",
)
STATUSES = ("unique", "duplicated", "ambiguous", "unaligned")
DEFAULT_MAX_BLOCKS = 5000


def _rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _ints(row, *keys):
    return tuple(int(row[k]) for k in keys)


def _node_dir(run_dir: Path, node: str, downward: bool = False) -> Path:
    return run_dir / ("downward" if downward else "") / "nodes" / node


def load_run(run_dir: str | Path, node: str) -> dict:
    """Load authoritative checkpoint tables into a browser-neutral model."""
    root = Path(run_dir)
    upward = _node_dir(root, node)
    if not (upward / "blocks.tsv").exists():
        raise ValueError(f"incomplete checkpoint: missing {upward / 'blocks.tsv'}")
    manifest = {}
    for candidate in (upward / "manifest.json", root / "manifest.json", root / "downward" / "manifest.json"):
        if candidate.exists():
            manifest.update(json.loads(candidate.read_text(encoding="utf-8")))
    down = _node_dir(root, node, True)
    blocks = _rows(upward / "blocks.tsv")
    atoms = _rows(upward / "atomic_intervals.tsv")
    occurrences = _rows(upward / "node_occurrences.tsv")
    leaves = _rows(upward / "leaf_occurrences.tsv")
    edge_runs = _rows(upward / "left_edge_runs.tsv") + _rows(upward / "right_edge_runs.tsv")
    conflicts = _rows(upward / "conflicts.tsv")
    merges = _rows(upward / "merge_decisions.tsv")
    down_blocks = _rows(down / "blocks.tsv")
    down_occurrences = _rows(down / "node_occurrences.tsv")
    down_leaves = _rows(down / "leaf_occurrences.tsv")
    breakpoints = []
    if (down / "breakpoints.json").exists():
        breakpoints = json.loads((down / "breakpoints.json").read_text(encoding="utf-8"))
    stats = {}
    for candidate in (upward / "summary.json", upward / "node_statistics.json", down / "statistics.json"):
        if candidate.exists():
            stats.update(json.loads(candidate.read_text(encoding="utf-8")))
    children = manifest.get("children") or manifest.get("child_hal_genomes") or []
    if not children:
        children = sorted({r.get("child_node", "") for r in edge_runs if r.get("child_node")})
    parent_chroms = sorted({r["chrom"] for r in blocks if r.get("chrom")})
    if not parent_chroms:
        parent_chroms = sorted({r["parent_chrom"] for r in edge_runs if r.get("parent_chrom")})
    return {
        "run_dir": str(root), "node": node, "manifest": manifest, "children": children,
        "blocks": blocks, "atoms": atoms, "occurrences": occurrences, "leaves": leaves,
        "edge_runs": edge_runs, "conflicts": conflicts, "merges": merges,
        "downward_blocks": down_blocks, "downward_occurrences": down_occurrences,
        "downward_leaves": down_leaves, "breakpoints": breakpoints, "stats": stats,
        "chromosomes": parent_chroms,
    }


def select_region(model: dict, chrom: str | None = None, start: int | None = None,
                  end: int | None = None, max_blocks: int = DEFAULT_MAX_BLOCKS) -> dict:
    blocks = model["blocks"]
    if chrom is None:
        candidates = [(b["chrom"], int(b["start"]), int(b["end"])) for b in blocks]
        if not candidates:
            raise ValueError("checkpoint contains no parent blocks")
        grouped = defaultdict(list)
        for c, s, e in candidates: grouped[c].append((s, e))
        scored = sorted(((len(v), c, min(x[0] for x in v), max(x[1] for x in v)) for c, v in grouped.items()), reverse=True)
        _, chrom, lo, hi = scored[0]
        # Prefer a compact deterministic preview while retaining multiple blocks.
        start, end = lo, min(hi, lo + 1_000_000)
        automatic = True
    else:
        automatic = False
        if chrom not in model["chromosomes"]: raise ValueError(f"unknown chromosome {chrom!r}")
        values = [b for b in blocks if b["chrom"] == chrom]
        if not values: raise ValueError(f"no blocks on chromosome {chrom!r}")
        start = min(int(b["start"]) for b in values) if start is None else start
        end = max(int(b["end"]) for b in values) if end is None else end
    if start is None or end is None or start < 0 or end <= start: raise ValueError("region must be non-empty half-open coordinates")
    selected = [b for b in blocks if b.get("chrom") == chrom and int(b["end"]) > start and int(b["start"]) < end]
    if len(selected) > max_blocks:
        raise ValueError(f"region contains {len(selected)} blocks, exceeding --max-blocks={max_blocks}; select a narrower region")
    return {"chrom": chrom, "start": start, "end": end, "automatic": automatic, "blocks": selected}


def _summary(model: dict) -> dict:
    classes = Counter(b.get("classification", "complex") for b in model["blocks"])
    statuses = Counter(r.get("status", "unique") for r in model["edge_runs"])
    sizes = sorted(int(b["end"]) - int(b["start"]) for b in model["blocks"])
    down = model["downward_blocks"]
    return {"upward_blocks": len(model["blocks"]), "atomic_blocks": len(model["atoms"]),
            "downward_blocks": len(down) if down else None, "classifications": {k: classes[k] for k in CLASSIFICATIONS},
            "mapping_status": {k: statuses[k] for k in STATUSES}, "size": {
                "min": min(sizes) if sizes else 0, "median": statistics.median(sizes) if sizes else 0,
                "max": max(sizes) if sizes else 0}, "split_blocks": sum(1 for b in model["blocks"] if ".down." in b["block_id"])}


def render_html(model: dict, region: dict, output: str | Path) -> None:
    payload = {"model": {k: model[k] for k in ("node", "children", "blocks", "atoms", "leaves", "edge_runs", "conflicts", "merges", "downward_blocks", "downward_leaves", "breakpoints", "manifest")}, "region": region, "summary": _summary(model), "classifications": CLASSIFICATIONS, "statuses": STATUSES}
    # JSON is embedded inside a script element; escape HTML-significant bytes
    # so hostile checkpoint labels cannot terminate the script or create markup.
    data = (json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
            .replace("<", "\\u003c").replace(">", "\\u003e")
            .replace("&", "\\u0026").replace("\\u2028", "\\u2028").replace("\\u2029", "\\u2029"))
    title = f"hal-multisynteny visualization: {model['node']}"
    page = f'''<!doctype html><html><head><meta charset="utf-8"><title>{html.escape(title)}</title>
<style>body{{font:14px system-ui;margin:0;color:#17202a}}header{{padding:12px 18px;background:#17324d;color:white}}main{{padding:14px}}button,select,input{{margin:3px;padding:5px}}#track{{border:1px solid #9aa;overflow:auto}}.row{{height:28px;position:relative;border-bottom:1px solid #ddd}}.item{{position:absolute;height:20px;top:3px;box-sizing:border-box;border:1px solid #222;overflow:hidden;font-size:10px;padding:2px}}#details{{white-space:pre-wrap;background:#f4f6f7;padding:10px;min-height:50px}}svg{{max-width:100%;border:1px solid #ddd}}</style></head><body><header><h2>{html.escape(title)}</h2><div id="meta"></div></header><main>
<label>View <select id="view"><option>upward</option><option>downward</option></select></label><label>Blocks <select id="kind"><option>merged</option><option>atomic</option></select></label><label>Classification <select id="class"><option value="">all</option></select></label><label>Status <select id="status"><option value="">all</option></select></label><button onclick="zoom(1.5)">Zoom in</button><button onclick="zoom(.67)">Zoom out</button><button onclick="reset()">Reset</button><input id="search" placeholder="search block ID" oninput="draw()"><div id="warning"></div><div id="track"></div><h3>Details</h3><div id="details">Select a block or atom.</div><h3>Overview</h3><svg id="plot" width="760" height="180"></svg></main>
<script>const D={data};let scale=1;const colors={{shared_consistent:'#0072B2',left_only:'#E69F00',right_only:'#56B4E9',orientation_conflict:'#D55E00',order_conflict:'#CC79A7',duplication_conflict:'#009E73',ambiguous:'#F0E442',unaligned:'#999999',complex:'#666666'}};const esc=s=>String(s??'').replace(/[&<>"']/g,c=>({{'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}}[c]));function reset(){{scale=1;draw()}}function zoom(x){{scale*=x;draw()}}function val(x){{return Number(x)}}function selected(){{let v=document.getElementById('view').value;let a=v==='downward'?D.model.downward_blocks:D.model.blocks;return (a||[]).filter(b=>(!D.region.chrom||b.chrom===D.region.chrom)&&val(b.end)>D.region.start&&val(b.start)<D.region.end&&(!document.getElementById('class').value||b.classification===document.getElementById('class').value)&&(!document.getElementById('search').value||b.block_id.includes(document.getElementById('search').value)))}}function draw(){{let a=selected(),w=Math.max(700,(D.region.end-D.region.start)*scale/2),h=Math.max(40,a.length*30);document.getElementById('track').style.width=w+'px';document.getElementById('track').innerHTML=a.map((b,i)=>{{let x=(val(b.start)-D.region.start)/(D.region.end-D.region.start)*w,y=i*30;let arrow=b.strand==='-'?'◀':'▶';return `<div class="row"><div class="item" style="left:${{x}}px;width:${{Math.max(3,(val(b.end)-val(b.start))/(D.region.end-D.region.start)*w)}}px;background:${{colors[b.classification]||'#ddd'}}" onclick='detail(${{JSON.stringify(b)}})'>${{esc(b.block_id)}} ${{arrow}}</div></div>`}}).join('');document.getElementById('meta').textContent=`node ${{D.model.node}}; children ${{D.model.children.join(', ')}}; region ${{D.region.chrom}}:${{D.region.start}}-${{D.region.end}}${{D.region.automatic?' (automatic preview)':''}}`;document.getElementById('warning').textContent=D.summary.downward_blocks===null?'Downward refinement unavailable; showing upward evidence only.':'Missing or one-sided evidence is not a confirmed deletion.';document.getElementById('details').textContent=a.length?`Showing ${{a.length}} retained objects; click one for authoritative fields.`:'No objects match the filters.'}}function detail(b){{document.getElementById('details').textContent=JSON.stringify(b,null,2)}}for(const c of D.classifications)document.getElementById('class').insertAdjacentHTML('beforeend',`<option>${{c}}</option>`);for(const s of D.statuses)document.getElementById('status').insertAdjacentHTML('beforeend',`<option>${{s}}</option>`);document.getElementById('view').onchange=draw;document.getElementById('kind').onchange=draw;draw();</script></body></html>'''
    Path(output).write_text(page, encoding="utf-8")


def visualize_run(run_dir: str | Path, node: str, output: str | Path, chrom: str | None = None,
                  start: int | None = None, end: int | None = None, max_blocks: int = DEFAULT_MAX_BLOCKS) -> dict:
    model = load_run(run_dir, node)
    region = select_region(model, chrom, start, end, max_blocks)
    render_html(model, region, output)
    return {"node": node, "output": str(output), "region": region, "summary": _summary(model)}
