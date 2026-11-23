## What this PR adds
<!-- e.g. new finding OSTRA-NNN, verification of OSTRA-NNN, new candidates, new probe -->

## Checklist
- [ ] Work is in my own branch (and my own `agents/<name>/` folder for raw exploration)
- [ ] Each finding folder follows `templates/` (finding.md, environment.txt, repro/, evidence/, verification/, report/)
- [ ] `finding.md` records the state at detection time; later changes go under "Corrections"
- [ ] Reproducers run with only the library installed; raw outputs are in `evidence/`
- [ ] `python3 tools/make_index.py && python3 tools/make_graph.py` were run and their outputs are committed
- [ ] `agents/LOG.md` updated if an AI agent proposed or verified something
- [ ] No secrets, no large binaries or third-party repositories committed
- [ ] Nothing was reported upstream (reports are filed only after review, see REPORTING.md)

## AI assistance
<!-- Which models/tools were used, for what, and what was checked by a human -->
