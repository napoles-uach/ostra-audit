import json, collections
from mols import MOLS, ATOMS, REACTIONS
H2K = 627.509474
E = {}
for d in map(json.loads, open('res/sc.jsonl')):
    if d.get('E') is not None and d.get('conv'):
        fam = 'r4scan' if 'r4scan' in d['xc'] else 'revtpss'
        E[(fam, d['tag'], d['sp'])] = d['E']
def ae(fam, tag, m):
    from pyscf import gto
    atoms = [a.split()[0] for a in MOLS[m][0].split(';')]
    try: return (sum(E[(fam, tag, a)] for a in atoms) - E[(fam, tag, m)]) * H2K
    except KeyError: return None
def rx(fam, tag, r):
    try: return sum(c * E[(fam, tag, s)] for s, c in REACTIONS[r].items()) * H2K
    except KeyError: return None
out = {}
for fam, tags in (('r4scan', ('bundled', 'patched', 'conv', 'libxc')), ('revtpss', ('patched', 'conv', 'libxc'))):
    rows = []
    for m in MOLS:
        rows.append(('AE', m, {t: ae(fam, t, m) for t in tags}))
    for r in REACTIONS:
        rows.append(('RX', r, {t: rx(fam, t, r) for t in tags}))
    out[fam] = rows
json.dump(out, open('res/derived.json', 'w'), indent=1)
for fam, rows in out.items():
    print('==', fam)
    for kind, name, v in rows:
        print(f"{kind} {name:24s} " + "  ".join(f"{t}={v[t]:9.3f}" if v[t] is not None else f"{t}=   n/a   " for t in v))
