#!/usr/bin/env python3
"""Build the findings graph from each findings/OSTRA-NNN-*/finding.md.

Reads the fields `libraries`, `pattern`, `oracle_class`, `impact_level` and `report`,
plus the vocabulary and library dependencies in tools/graph_config.json.
Writes findings/graph.json and docs/graph.html (self-contained, d3 from cdnjs).
Run from the repository root:  python3 tools/make_graph.py
"""
import json, re, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_index import parse, ROOT, FINDINGS  # noqa: E402

CFG = json.loads((ROOT / "tools" / "graph_config.json").read_text())
ROMAN = {"I": 1, "II": 2, "III": 3, "IV": 4}


def split(v):
    v = (v or "").strip()
    return [] if v in ("", "none", "-") else [x.strip() for x in v.split(",") if x.strip()]


def main():
    errors, findings = [], []
    for md in sorted(FINDINGS.glob("OSTRA-*/finding.md")):
        r = parse(md)
        libs, pats, oras = split(r.get("libraries")), split(r.get("pattern")), split(r.get("oracle_class"))
        for x in libs:
            if x not in CFG["libraries"]: errors.append(f"{r['id']}: unknown library '{x}'")
        for x in pats:
            if x not in CFG["patterns"]: errors.append(f"{r['id']}: unknown pattern '{x}'")
        for x in oras:
            if x not in CFG["oracles"]: errors.append(f"{r['id']}: unknown oracle_class '{x}'")
        if not libs: errors.append(f"{r['id']}: missing 'libraries'")
        m = re.match(r"(IV|I{1,3})", r.get("impact_level") or "")
        rep = (r.get("report") or "").lower()
        findings.append({
            "id": r["id"], "num": r["id"].split("-")[1], "title": r["title"], "path": r["path"],
            "libraries": libs, "patterns": pats, "oracles": oras,
            "impact": ROMAN[m.group(1)] if m else 4, "impact_level": r.get("impact_level", ""),
            "status": r.get("status", ""),
            "reported": rep.startswith(("reported", "fixed", "commented")),
        })
    if errors:
        sys.exit("\n".join(errors))
    data = {"findings": findings, **CFG}
    (FINDINGS / "graph.json").write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")
    tpl = (ROOT / "tools" / "graph_template.html").read_text()
    out = ROOT / "docs" / "graph.html"
    out.parent.mkdir(exist_ok=True)
    out.write_text(tpl.replace("/*__DATA__*/null", json.dumps(data, ensure_ascii=False)))
    print(f"graph: {len(findings)} findings -> findings/graph.json, docs/graph.html")


if __name__ == "__main__":
    main()
