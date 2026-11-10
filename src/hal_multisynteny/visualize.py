"""Checkpoint discovery, normalized evidence model, and local HTML viewer.

The viewer is deliberately read-only: it consumes recorded TSV/JSON evidence
and never invokes HAL or changes reconciliation behavior.
"""

from __future__ import annotations

import csv
import json
import statistics
from collections import Counter
from pathlib import Path

CLASSIFICATIONS = ("shared_consistent", "left_only", "right_only", "orientation_conflict",
                   "order_conflict", "duplication_conflict", "ambiguous", "unaligned", "complex")
STATUSES = ("unique", "duplicated", "ambiguous", "unaligned")
DEFAULT_MAX_BLOCKS = 5000


def _read(path: Path, required=()):
    if not path.exists():
        if required:
            raise ValueError(f"incomplete checkpoint: missing {path}")
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if required:
        missing = sorted(set(required) - set(rows[0] if rows else []))
        if missing:
            raise ValueError(f"incompatible checkpoint {path}: missing columns {', '.join(missing)}")
    return rows


def _node_dir(run_dir: Path, node: str, downward=False) -> Path:
    return run_dir / ("downward" if downward else "") / "nodes" / node


def _ints(row, *keys):
    return tuple(int(row[key]) for key in keys)


def _json(path: Path, default):
    if not path.exists():
        return default
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid JSON checkpoint {path}: {exc}") from exc


def _normal_block(row, kind, node):
    start_key, end_key = ("start", "end") if "start" in row else ("parent_start", "parent_end")
    chrom_key = "chrom" if "chrom" in row else "parent_chrom"
    return {"object_type": kind, "block_id": row.get("block_id", ""), "node": node,
            "chrom": row.get(chrom_key, ""), "start": int(row[start_key]), "end": int(row[end_key]),
            "classification": row.get("classification", ""), "mapping_status": row.get("mapping_status", "unique"),
            "record": row}


def load_run(run_dir: str | Path, node: str) -> dict:
    """Discover and validate a run, then normalize authoritative evidence tables."""
    root = Path(run_dir)
    up = _node_dir(root, node)
    blocks_rows = _read(up / "blocks.tsv", ("block_id", "node", "chrom", "start", "end", "classification", "mapping_status"))
    atoms_rows = _read(up / "atomic_intervals.tsv")
    if not blocks_rows:
        raise ValueError(f"incomplete checkpoint: missing {up / 'blocks.tsv'}")
    manifest = {}
    for path in (up / "manifest.json", root / "manifest.json", root / "downward" / "manifest.json"):
        manifest.update(_json(path, {}))
    edge = []
    for side in ("left", "right"):
        for row in _read(up / f"{side}_edge_runs.tsv"):
            row = dict(row); row["side"] = side; row["object_type"] = "mapping"
            row["parent_start"] = int(row["parent_start"]); row["parent_end"] = int(row["parent_end"])
            row["child_start"] = int(row["child_start"]); row["child_end"] = int(row["child_end"])
            edge.append(row)
    occurrences = _read(up / "node_occurrences.tsv")
    leaves = _read(up / "leaf_occurrences.tsv")
    conflicts = _read(up / "conflicts.tsv")
    merges = _read(up / "merge_decisions.tsv")
    down = _node_dir(root, node, True)
    down_blocks = _read(down / "blocks.tsv")
    down_occurrences = _read(down / "node_occurrences.tsv")
    down_leaves = _read(down / "leaf_occurrences.tsv")
    lineage = _json(down / "lineage.json", [])
    breakpoints = _json(down / "breakpoints.json", [])
    children = manifest.get("children") or manifest.get("child_hal_genomes") or sorted({r["child_node"] for r in edge})
    blocks = [_normal_block(row, "upward_block", node) for row in blocks_rows]
    atoms = []
    for row in atoms_rows:
        if not row.get("parent_chrom"):
            continue
        atoms.append({"object_type": "atomic_interval", "atom_id": row.get("atomic_interval_id", ""),
                      "block_id": row.get("final_parent_block_id", ""), "node": node,
                      "chrom": row["parent_chrom"], "start": int(row["parent_start"]),
                      "end": int(row["parent_end"]), "classification": row.get("classification", ""),
                      "mapping_status": "ambiguous" if row.get("classification") == "ambiguous" else "unique",
                      "record": row})
    def occ(row, kind):
        return {"object_type": kind, "block_id": row.get("block_id", ""), "node": row.get("node", node),
                "chrom": row.get("node_chrom", ""), "start": int(row.get("node_start", 0)),
                "end": int(row.get("node_end", 0)), "classification": "", "mapping_status": row.get("status", "unique"),
                "side": row.get("leaf", ""), "occurrence_id": row.get("occurrence_id", ""),
                "copy_id": row.get("copy_id", ""), "record": row}
    occs = [occ(r, "node_occurrence") for r in occurrences]
    leaf_objs = [occ(r, "leaf_occurrence") for r in leaves]
    down_objs = [_normal_block(r, "downward_block", node) for r in down_blocks]
    boundary_objs = []
    for event in breakpoints:
        if event.get("parent_chrom") and event.get("parent_coordinate") is not None:
            boundary_objs.append({"object_type": "downward_boundary", "block_id": event.get("origin_block", ""),
                                  "node": node, "chrom": event["parent_chrom"], "start": int(event["parent_coordinate"]),
                                  "end": int(event["parent_coordinate"]), "classification": "",
                                  "mapping_status": "ambiguous" if event.get("ambiguous") else "unique", "record": event})
    parent_chroms = sorted({b["chrom"] for b in blocks})
    return {"run_dir": str(root), "node": node, "children": children, "manifest": manifest,
            "blocks": blocks, "atoms": atoms, "mappings": edge, "occurrences": occs + leaf_objs,
            "conflicts": conflicts, "merges": merges, "downward_blocks": down_objs,
            "downward_occurrences": down_occurrences, "downward_leaves": down_leaves,
            "lineage": lineage, "breakpoints": breakpoints, "boundaries": boundary_objs,
            "chromosomes": parent_chroms, "stats": _json(up / "summary.json", {})}


def _intersects(obj, chrom, start, end):
    return obj.get("chrom") == chrom and obj.get("end", 0) > start and obj.get("start", 0) < end


def _region_objects(model, chrom, start, end):
    groups = {"atoms": model["atoms"], "blocks": model["blocks"], "mappings": model["mappings"],
              "occurrences": model["occurrences"], "downward_blocks": model["downward_blocks"],
              "boundaries": model["boundaries"]}
    selected = {}
    for key, rows in groups.items():
        if key == "mappings":
            selected[key] = [r for r in rows if r.get("parent_chrom") == chrom and r["parent_end"] > start and r["parent_start"] < end]
        elif key == "boundaries":
            selected[key] = [r for r in rows if r.get("chrom") == chrom and start <= r.get("start", 0) <= end]
        else:
            selected[key] = [r for r in rows if _intersects(r, chrom, start, end)]
    return selected


def select_region(model: dict, chrom: str | None = None, start: int | None = None,
                  end: int | None = None, max_blocks: int = DEFAULT_MAX_BLOCKS) -> dict:
    blocks = model["blocks"]
    if not blocks:
        raise ValueError("checkpoint contains no parent blocks")
    if chrom is None:
        candidates = []
        for c in model["chromosomes"]:
            rows = [b for b in blocks if b["chrom"] == c]
            if not rows: continue
            lo, hi = min(b["start"] for b in rows), max(b["end"] for b in rows)
            region = _region_objects(model, c, lo, min(hi, lo + 1_000_000))
            score = (sum(1 for b in region["blocks"] if b["classification"] not in {"", "left_only", "right_only"}),
                     len(region["boundaries"]), len({b["classification"] for b in region["blocks"]}), len(region["blocks"]))
            candidates.append((score, c, lo, min(hi, lo + 1_000_000)))
        _, chrom, start, end = max(candidates)
        reason = "informative preview (boundary/conflict, then classification diversity, then density)"
        automatic = True
    else:
        if chrom not in model["chromosomes"]: raise ValueError(f"unknown chromosome {chrom!r}")
        rows = [b for b in blocks if b["chrom"] == chrom]
        start = min(b["start"] for b in rows) if start is None else start
        end = max(b["end"] for b in rows) if end is None else end
        reason, automatic = "user-selected region", False
    if start is None or end is None or start < 0 or end <= start:
        raise ValueError("region must be non-empty half-open coordinates")
    selected = _region_objects(model, chrom, start, end)
    counts = {k: len(v) for k, v in selected.items()}
    total = sum(counts.values())
    if total > max_blocks:
        raise ValueError(f"region {chrom}:{start}-{end} contains {total} rendered objects ({counts}), exceeding --max-blocks={max_blocks}; select a narrower region")
    return {"chrom": chrom, "start": start, "end": end, "automatic": automatic, "reason": reason,
            "counts": counts, "objects": selected, "blocks": selected["blocks"]}


def _summary(model):
    classes = Counter(b.get("classification", "complex") for b in model["blocks"])
    statuses = Counter(r.get("status", "unique") for r in model["mappings"])
    sizes = sorted(b["end"] - b["start"] for b in model["blocks"])
    lineage_counts = Counter(row.get("original_block_id") for row in model["lineage"] if row.get("original_block_id"))
    split = sum(count > 1 for count in lineage_counts.values())
    return {"upward_blocks": len(model["blocks"]), "atomic_blocks": len(model["atoms"]),
            "downward_blocks": len(model["downward_blocks"]) or None, "split_blocks": split,
            "classifications": {key: classes[key] for key in CLASSIFICATIONS},
            "mapping_status": {key: statuses[key] for key in STATUSES},
            "size": {"min": min(sizes) if sizes else 0, "median": statistics.median(sizes) if sizes else 0,
                     "max": max(sizes) if sizes else 0}}


def _safe_json(payload):
    return (json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
            .replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
            .replace("\\u2028", "\\u2028").replace("\\u2029", "\\u2029"))


def render_html(model: dict, region: dict, output: str | Path) -> None:
    summary = _summary(model)
    payload = {"node": model["node"], "children": model["children"], "chromosomes": model["chromosomes"], "region": region,
               "summary": summary, "objects": region["objects"], "merges": model["merges"],
               "conflicts": model["conflicts"], "lineage": model["lineage"],
               "manifest": model["manifest"], "classifications": CLASSIFICATIONS, "statuses": STATUSES}
    data = _safe_json(payload)
    page = """<!doctype html><html><head><meta charset="utf-8"><title>hal-multisynteny run viewer</title>
<style>body{font:14px system-ui;margin:0;color:#17202a}header{padding:12px 18px;background:#17324d;color:white}main{padding:14px}label{margin-right:8px}button,select,input{margin:3px;padding:5px}#track{border:1px solid #9aa;overflow:auto}.lane{height:32px;position:relative;border-bottom:1px solid #ddd}.lane b{display:inline-block;width:175px}.item{position:absolute;height:20px;top:5px;box-sizing:border-box;border:1px solid #222;overflow:hidden;font-size:10px;padding:2px;cursor:pointer}.mapping{opacity:.78}.reverse{border-style:dashed}.boundary{position:absolute;height:100%;border-left:2px dotted #111}.gap{background:repeating-linear-gradient(45deg,#fff,#fff 3px,#ddd 3px,#ddd 5px)}#details{white-space:pre-wrap;background:#f4f6f7;padding:10px;min-height:80px;max-height:360px;overflow:auto}svg{max-width:100%;border:1px solid #ddd}.warning{padding:8px;background:#fff8d5}</style></head><body><header><h2>hal-multisynteny checkpoint viewer</h2><div id="meta"></div></header><main>
<div><label>Chromosome <select id="chrom"></select></label><label>Start <input id="start" type="number"></label><label>End <input id="end" type="number"></label><button id="apply">Apply</button><button id="zin">Zoom in</button><button id="zout">Zoom out</button><button id="reset">Reset</button></div>
<div><label>View <select id="view"><option value="upward">Upward</option><option value="downward">Downward</option></select></label><label>Track <select id="kind"><option value="merged">Merged blocks</option><option value="atomic">Atomic intervals</option></select></label><label>Classification <select id="class"><option value="">All</option></select></label><label>Status <select id="status"><option value="">All</option></select></label><label><input id="lanes" type="checkbox" checked> Child evidence</label><label><input id="bounds" type="checkbox" checked> Boundaries</label><input id="search" placeholder="block/atom ID"></div>
<div class="warning" id="warning"></div><div id="track"></div><h3>Details</h3><div id="details">Select any rendered object.</div><h3>Overview (recorded or viewer-computed)</h3><svg id="plot" width="820" height="260" role="img" aria-label="overview plots"></svg></main>
<script>""" + "const D=" + data + r""";
const C={shared_consistent:'#0072B2',left_only:'#E69F00',right_only:'#56B4E9',orientation_conflict:'#D55E00',order_conflict:'#CC79A7',duplication_conflict:'#009E73',ambiguous:'#F0E442',unaligned:'#999999',complex:'#666'};let scale=1, state={...D.region};
const $=id=>document.getElementById(id); const val=(x)=>Number(x); const text=(x)=>String(x??'');
function all(){let k=$("kind").value,v=$("view").value,o=D.objects;let a=k==='atomic'?o.atoms:(v==='downward'?o.downward_blocks:o.blocks);let filt=$("class").value,s=$("status").value,q=$("search").value.toLowerCase();return (a||[]).filter(x=>val(x.end)>state.start&&val(x.start)<state.end&&(!filt||x.classification===filt)&&(!s||x.mapping_status===s)&&(!q||text(x.block_id||x.atom_id).toLowerCase().includes(q)))}
function evidence(){return $("lanes").checked?D.objects.mappings.filter(x=>x.parent_chrom===state.chrom&&x.parent_end>state.start&&x.parent_start<state.end):[]}
function detail(x){$("details").textContent=JSON.stringify(x,null,2)}
function draw(){let w=Math.max(850,(state.end-state.start)*scale/2), rows=[];let add=(name,items,coord)=>{rows.push('<div class="lane"><b>'+name+'</b>'+items.map(x=>{let s=coord?x.parent_start:x.start,e=coord?x.parent_end:x.end;if(e<=state.start||s>=state.end)return '';let left=(s-state.start)/(state.end-state.start)*w,width=Math.max(3,(Math.min(e,state.end)-Math.max(s,state.start))/(state.end-state.start)*w);let cls=x.strand==='-'?' reverse':'';let label=text(x.block_id||x.atom_id||x.child_block_id);return '<div class="item'+cls+'" data-json="'+encodeURIComponent(JSON.stringify(x))+'" style="left:'+left+'px;width:'+width+'px;background:'+(C[x.classification]||'#ddd')+'">'+label+(x.strand==='-'?' ◀':' ▶')+'</div>'}).join('')+'</div>')};let a=all();add('Parent '+($("kind").value==='atomic'?'atoms':'blocks'),a,false);if($("lanes").checked){let by={};for(const x of evidence())(by[x.side]??=[]).push(x);for(const [side,items] of Object.entries(by))add((side||'child')+' evidence',items,true)}if($("view").value==='downward')add('Downward',D.objects.downward_blocks,false);if($("bounds").checked)add('Propagated boundaries',D.objects.boundaries,false);$("track").style.width=w+'px';$("track").innerHTML=rows.join('');$("track").querySelectorAll('.item').forEach(e=>e.addEventListener('click',()=>detail(JSON.parse(decodeURIComponent(e.dataset.json)))));$("meta").textContent='node '+D.node+'; children '+D.children.join(', ')+'; '+state.chrom+':'+state.start+'-'+state.end+(D.region.automatic?' (automatic: '+D.region.reason+')':'');$("warning").textContent='Missing or one-sided evidence is not a confirmed biological deletion.'+(D.summary.downward_blocks===null?' Downward refinement unavailable.':'');drawPlot()}
function drawPlot(){let svg=$("plot"),vals=[['upward',D.summary.upward_blocks],['atomic',D.summary.atomic_blocks],['downward',D.summary.downward_blocks||0],['split',D.summary.split_blocks]],m=Math.max(...vals.map(x=>x[1]),1);svg.innerHTML=vals.map((x,i)=>'<text x="5" y="'+(25+i*42)+'">'+x[0]+' '+x[1]+'</text><rect x="120" y="'+(10+i*42)+'" width="'+(650*x[1]/m)+'" height="22" fill="#0072B2"/>').join('')}
function fill(){D.region.chrom&&$('chrom').insertAdjacentHTML('beforeend','<option>'+D.region.chrom+'</option>');for(const x of D.summary?D.region.chromosomes||[]:[])$('chrom').insertAdjacentHTML('beforeend','<option>'+x+'</option>');$('chrom').value=state.chrom;$('start').value=state.start;$('end').value=state.end;for(const c of D.summary?D.summary.classifications?Object.keys(D.summary.classifications):[]:[])$('class').insertAdjacentHTML('beforeend','<option>'+c+'</option>');for(const s of D.summary?D.summary.mapping_status?Object.keys(D.summary.mapping_status):[]:[])$('status').insertAdjacentHTML('beforeend','<option>'+s+'</option>')}
function apply(){state={chrom:$('chrom').value,start:val($('start').value),end:val($('end').value)};draw()}$("apply").addEventListener('click',apply);$("zin").addEventListener('click',()=>{let c=(state.start+state.end)/2;let d=(state.end-state.start)/3;state.start=Math.max(0,Math.floor(c-d));state.end=Math.ceil(c+d);draw()});$("zout").addEventListener('click',()=>{let c=(state.start+state.end)/2;let d=(state.end-state.start)*1.5;state.start=Math.max(0,Math.floor(c-d));state.end=Math.ceil(c+d);draw()});$("reset").addEventListener('click',()=>{state={chrom:D.region.chrom,start:D.region.start,end:D.region.end};draw()});for(const id of ['view','kind','class','status','lanes','bounds','search'])$(id).addEventListener('input',draw);fill();draw();</script></body></html>"""
    Path(output).write_text(page, encoding="utf-8")


def visualize_run(run_dir, node, output, chrom=None, start=None, end=None, max_blocks=DEFAULT_MAX_BLOCKS):
    model = load_run(run_dir, node)
    region = select_region(model, chrom, start, end, max_blocks)
    render_html(model, region, output)
    return {"node": node, "output": str(output), "region": region, "summary": _summary(model)}
