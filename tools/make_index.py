#!/usr/bin/env python3
"""Regenerate findings/INDEX.md and findings/index.json from each findings/OSTRA-NNN-*/finding.md.

finding.md is the single source of truth. Run from the repository root:  python3 tools/make_index.py
"""
import json, re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FINDINGS = ROOT / "findings"
HEADER = """# Findings index

Impact levels: I = changes publishable observables; II = breaks numerical consistency (energy/force/potential); III = niche but severe or silent when hit; IV = no impact. See `../PROTOCOL.md`.

| ID | Detected | Library | Finding | Class | Impact | Status | Report |
|---|---|---|---|---|---|---|---|
"""


def parse(md: Path):
    text = md.read_text()
    m = re.match(r"#\s*(OSTRA-\d+):\s*(.+)", text)
    rec = {"id": m.group(1), "title": m.group(2).strip()}
    for line in text.splitlines():
        f = re.match(r"\|\s*([a-z_]+)\s*\|\s*(.*?)\s*\|\s*$", line)
        if f and f.group(1) not in ("field",):
            rec.setdefault(f.group(1), f.group(2))
    rec["path"] = str(md.parent.relative_to(ROOT))
    rec["slug"] = md.parent.name[len(rec["id"]) + 1:]
    return rec


def cell(s):
    return (s or "").replace("|", "/")


def main():
    recs = [parse(p) for p in sorted(FINDINGS.glob("OSTRA-*/finding.md"))]
    rows = []
    for r in recs:
        lvl = r.get("impact_level") or ""
        m = re.match(r"[IV]+(?:\s*-\s*[IV]+)?", lvl)
        impact = m.group(0).replace(" ", "") if m else lvl.split(" (")[0]
        link = f"[{r['id']}]({Path(r['path']).name}/finding.md)"
        rows.append("| " + " | ".join(cell(x) for x in (
            link, r.get("detected"), r.get("library"), r["title"], r.get("class"),
            impact, r.get("status"), r.get("report"))) + " |")
    (FINDINGS / "INDEX.md").write_text(HEADER + "\n".join(rows) + "\n")
    (FINDINGS / "index.json").write_text(json.dumps(recs, indent=2, ensure_ascii=False) + "\n")
    print(f"{len(recs)} findings indexed")


if __name__ == "__main__":
    main()
