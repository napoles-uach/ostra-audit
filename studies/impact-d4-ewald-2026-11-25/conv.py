# Is the standard (primitive) cell itself converged? n x n x n supercells contain more Ewald images.
import numpy as np, json
from ase.build import bulk, make_supercell
from cells import d4, H2KJ
out = {}
for name, kw, nfa in (("MgO", dict(name="MgO", crystalstructure="rocksalt", a=4.21), 2),
                      ("NaCl", dict(name="NaCl", crystalstructure="rocksalt", a=5.64), 2),
                      ("CaF2", dict(name="CaF2", crystalstructure="fluorite", a=5.46), 3)):
    for cell in ("primitive", "cubic"):
        b = bulk(**kw, cubic=(cell == "cubic"))
        for n in (1, 2, 3, 4):
            at = make_supercell(b, n * np.eye(3))
            if len(at) > 600: continue
            e, q = d4(at)
            out[f"{name} {cell} {n}"] = (e / (len(at) / nfa) * H2KJ, float(q[0]))
            print(f"{name:5s} {cell:9s} {n}x{n}x{n} ({len(at):3d} atoms)  E = {out[f'{name} {cell} {n}'][0]:.5f} kJ/mol/f.u.  q = {q[0]:+.5f}", flush=True)
json.dump(out, open("supercell_conv.json", "w"), indent=1)
