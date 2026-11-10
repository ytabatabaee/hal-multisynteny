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
    # Attach parent classification to child evidence by recorded parent overlap;
    # this is a display join, never a new biological classification.
    for mapping in edge:
        covering = [b for b in blocks if b["chrom"] == mapping["parent_chrom"]
                    and b["end"] > mapping["parent_start"] and b["start"] < mapping["parent_end"]]
        mapping["classification"] = covering[0]["classification"] if len(covering) == 1 else ""
    for occurrence in occs + leaf_objs:
        covering = [b for b in blocks if b["chrom"] == occurrence["chrom"]
                    and b["end"] > occurrence["start"] and b["start"] < occurrence["end"]]
        occurrence["classification"] = covering[0]["classification"] if len(covering) == 1 else ""
    for block in blocks:
        block["merge_decisions"] = [m for m in merges if m.get("parent_chrom") == block["chrom"] and
                                     int(m.get("left_start", 0)) < block["end"] and int(m.get("right_end", 0)) > block["start"]]
        block["conflicts"] = [c for c in conflicts if c.get("parent_chrom") == block["chrom"] and
                               int(c.get("parent_start", 0)) < block["end"] and int(c.get("parent_end", 0)) > block["start"]]
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
    loaded = [b for b in blocks if b["chrom"] == chrom]
    loaded_start, loaded_end = min(b["start"] for b in loaded), max(b["end"] for b in loaded)
    if start < loaded_start or end > loaded_end:
        raise ValueError(f"requested region {chrom}:{start}-{end} exceeds loaded extent {chrom}:{loaded_start}-{loaded_end}")
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
            "split_fraction": split / len(model["blocks"]) if model["blocks"] else 0.0,
            "left_mapping_runs": sum(r.get("side") == "left" for r in model["mappings"]),
            "right_mapping_runs": sum(r.get("side") == "right" for r in model["mappings"]),
            "occurrences": len(model["occurrences"]), "conflicts": len(model["conflicts"]),
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
<div class="warning" id="warning"></div><div id="track"></div><h3>Details</h3><div id="details">Select any rendered object.</div><h3>Overview (recorded or viewer-computed)</h3><svg id="plot" width="820" height="300" role="img" aria-label="overview plots"></svg></main>
<script>""" + "const D=" + data + r""";
const C={shared_consistent:'#0072B2',left_only:'#E69F00',right_only:'#56B4E9',orientation_conflict:'#D55E00',order_conflict:'#CC79A7',duplication_conflict:'#009E73',ambiguous:'#F0E442',unaligned:'#999999',complex:'#666'};let scale=1,state={...D.region};
const $=id=>document.getElementById(id), text=x=>String(x??''), num=x=>Number(x);
function detail(x){let box=$('details');box.replaceChildren();let lead=document.createElement('div');let id=x.block_id||x.atom_id||x.child_block_id||'';let coords=x.parent_chrom?x.parent_chrom+':'+x.parent_start+'-'+x.parent_end:(x.chrom||'')+':'+x.start+'-'+x.end;lead.textContent=(x.object_type||'object')+' '+id+'; '+coords+'; classification='+(x.classification||'recorded field absent')+'; mapping_status='+(x.mapping_status||x.status||'recorded field absent')+'; strand='+(x.strand||'recorded field absent')+'; occurrence='+(x.occurrence_id||'recorded field absent')+'; copy='+(x.copy_id||'recorded field absent');box.append(lead);let raw=document.createElement('details'),sum=document.createElement('summary');sum.textContent='Raw authoritative record';let pre=document.createElement('pre');pre.textContent=JSON.stringify(x.record||x,null,2);raw.append(sum,pre);box.append(raw)}
function matches(x){let q=$('search').value.toLowerCase(),c=$('class').value,s=$('status').value;return (!c||x.classification===c)&&(!s||x.mapping_status===s)&&(!q||text(x.block_id||x.atom_id||x.child_block_id).toLowerCase().includes(q))}
function items(){let o=D.objects,k=$('kind').value,v=$('view').value;let a=k==='atomic'?o.atoms:(v==='downward'?o.downward_blocks:o.blocks);return (a||[]).filter(x=>num(x.end)>state.start&&num(x.start)<state.end&&matches(x))}
function evidence(){return $('lanes').checked?D.objects.mappings.filter(x=>x.parent_chrom===state.chrom&&x.parent_end>state.start&&x.parent_start<state.end&&matches(x)):[]}
function lane(name,items,coord,w){let row=document.createElement('div');row.className='lane';let label=document.createElement('b');label.textContent=name;row.append(label);for(const x of items){let s=coord?x.parent_start:x.start,e=coord?x.parent_end:x.end;if(e<=state.start||s>=state.end)continue;let el=document.createElement('button');el.type='button';el.className='item'+(x.strand==='-'?' reverse':'');el.textContent=text(x.block_id||x.atom_id||x.child_block_id)+(x.strand==='-'?' ◀':' ▶');el.style.left=((s-state.start)/(state.end-state.start)*w)+'px';el.style.width=Math.max(3,(Math.min(e,state.end)-Math.max(s,state.start))/(state.end-state.start)*w)+'px';el.style.background=C[x.classification]||'#ddd';el.addEventListener('click',()=>detail(x));row.append(el)}return row}
function draw(){let w=Math.max(850,(state.end-state.start)*scale/2),track=$('track');track.replaceChildren();track.style.width=w+'px';let a=items();track.append(lane('Parent '+($('kind').value==='atomic'?'atoms':($('view').value==='downward'?'downward blocks':'blocks')),a,false,w));if($('lanes').checked){let by={};for(const x of evidence())(by[x.side||'child']??=[]).push(x);for(const [side,rows] of Object.entries(by))track.append(lane(side+' evidence',rows,true,w))}if($('view').value==='upward'&&D.objects.downward_blocks.length)track.append(lane('Downward blocks',D.objects.downward_blocks.filter(matches),false,w));if($('bounds').checked)track.append(lane('Propagated boundaries',D.objects.boundaries,true,w));$('meta').textContent='node '+D.node+'; children '+D.children.join(', ')+'; '+state.chrom+':'+state.start+'-'+state.end+(D.region.automatic?' (automatic: '+D.region.reason+')':'');$('warning').textContent='Loaded extent is '+D.region.chrom+':'+D.region.start+'-'+D.region.end+'; edits stay inside it. Missing or one-sided evidence is not a confirmed biological deletion.'+(D.summary.downward_blocks===null?' Downward refinement unavailable.':'');drawPlot()}
function drawPlot(){let svg=$('plot');svg.replaceChildren();let vals=[['atomic intervals',D.summary.atomic_blocks],['reconciled blocks',D.summary.upward_blocks],['downward blocks',D.summary.downward_blocks||0],['split blocks',D.summary.split_blocks],['left mappings',D.summary.left_mapping_runs],['right mappings',D.summary.right_mapping_runs],['occurrences',D.summary.occurrences],['conflicts',D.summary.conflicts]],m=Math.max(...vals.map(x=>x[1]),1);for(let i=0;i<vals.length;i++){let t=document.createElementNS('http://www.w3.org/2000/svg','text');t.setAttribute('x','5');t.setAttribute('y',String(20+i*30));t.textContent=vals[i][0]+' '+vals[i][1];svg.append(t);let r=document.createElementNS('http://www.w3.org/2000/svg','rect');r.setAttribute('x','180');r.setAttribute('y',String(7+i*30));r.setAttribute('width',String(600*vals[i][1]/m));r.setAttribute('height','18');r.setAttribute('fill','#0072B2');svg.append(r)}let c=document.createElementNS('http://www.w3.org/2000/svg','text');c.setAttribute('x','5');c.setAttribute('y','265');c.textContent='classifications: '+Object.entries(D.summary.classifications).map(x=>x[0]+'='+x[1]).join(', ');svg.append(c)}
function addOption(sel,value,label){let o=document.createElement('option');o.value=value;o.textContent=label||value;sel.append(o)}
function fill(){let chrom=$('chrom');chrom.replaceChildren();for(const c of D.chromosomes)addOption(chrom,c);chrom.value=state.chrom;$('start').value=state.start;$('end').value=state.end;for(const c of D.classifications)addOption($('class'),c);for(const s of D.statuses)addOption($('status'),s);$('view').options[1].disabled=!D.objects.downward_blocks.length;if(D.chromosomes.length===1){chrom.disabled=true;chrom.title='Only one chromosome is embedded in this bounded preview.'}}
function apply(){let c=$('chrom').value,s=num($('start').value),e=num($('end').value);if(c!==D.region.chrom||s<D.region.start||e>D.region.end||e<=s){$('warning').textContent='Requested region is outside the embedded bounded preview '+D.region.chrom+':'+D.region.start+'-'+D.region.end+'.';return}state={chrom:c,start:s,end:e};draw()}
$('apply').addEventListener('click',apply);$('zin').addEventListener('click',()=>{let c=(state.start+state.end)/2,d=(state.end-state.start)/3;state.start=Math.max(D.region.start,Math.floor(c-d));state.end=Math.min(D.region.end,Math.ceil(c+d));draw()});$('zout').addEventListener('click',()=>{let c=(state.start+state.end)/2,d=(state.end-state.start)*1.5;state.start=Math.max(D.region.start,Math.floor(c-d));state.end=Math.min(D.region.end,Math.ceil(c+d));draw()});$('reset').addEventListener('click',()=>{state={chrom:D.region.chrom,start:D.region.start,end:D.region.end};draw()});for(const id of ['view','kind','class','status','lanes','bounds','search'])$(id).addEventListener('input',draw);fill();draw();</script></body></html>"""
    Path(output).write_text(page, encoding="utf-8")


def visualize_run(run_dir, node, output, chrom=None, start=None, end=None, max_blocks=DEFAULT_MAX_BLOCKS):
    model = load_run(run_dir, node)
    region = select_region(model, chrom, start, end, max_blocks)
    render_html(model, region, output)
    return {"node": node, "output": str(output), "region": region, "summary": _summary(model)}
