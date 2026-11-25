# OSTRA-005 impact: D4 dispersion energy and EEQ charges for equivalent descriptions of the same crystal.
# Every cell below describes the same infinite crystal; physically all must give identical energy per formula unit.
import json, numpy as np
from ase.build import bulk, make_supercell, niggli_reduce
from dftd4.interface import DispersionModel, DampingParam
AU = 1 / 0.52917721067  # Angstrom -> bohr
H2KJ = 2625.4996
MATS = {  # name: (ase.build.bulk args, atoms per formula unit)
 "NaCl (rocksalt)":   (dict(name="NaCl", crystalstructure="rocksalt", a=5.64), 2),
 "MgO (rocksalt)":    (dict(name="MgO", crystalstructure="rocksalt", a=4.21), 2),
 "LiF (rocksalt)":    (dict(name="LiF", crystalstructure="rocksalt", a=4.03), 2),
 "CsCl":              (dict(name="CsCl", crystalstructure="cesiumchloride", a=4.12), 2),
 "CaF2 (fluorite)":   (dict(name="CaF2", crystalstructure="fluorite", a=5.46), 3),
 "ZnO (wurtzite)":    (dict(name="ZnO", crystalstructure="wurtzite", a=3.25, c=5.21), 2),
 "GaN (wurtzite)":    (dict(name="GaN", crystalstructure="wurtzite", a=3.19, c=5.19), 2),
 "Si (diamond)":      (dict(name="Si", crystalstructure="diamond", a=5.43), 1),
}
def variants(kw):
    prim = bulk(**kw)
    out = {"primitive": prim}
    try: out["conventional cubic"] = bulk(**kw, cubic=True)
    except Exception: pass
    if kw["crystalstructure"] == "wurtzite":   # gamma = 120 deg (default) vs the equivalent 60 deg setting
        out["hexagonal, gamma = 60"] = make_supercell(prim, [[1, 0, 0], [1, 1, 0], [0, 0, 1]])
    n = prim.copy(); niggli_reduce(n); out["Niggli-reduced"] = n
    out["2x2x2 supercell"] = make_supercell(prim, 2 * np.eye(3))
    for k in (1, 2):
        out[f"a' = a + {k}b"] = make_supercell(prim, [[1, k, 0], [0, 1, 0], [0, 0, 1]])
    return out
def d4(at):
    m = DispersionModel(at.numbers, at.positions * AU, lattice=np.array(at.cell) * AU, periodic=np.array([True] * 3))
    q = m.get_properties()["partial charges"]
    e = m.get_dispersion(DampingParam(method="pbe0"), grad=False)["energy"]
    return e, q
def main():
  res = {}
  for name, (kw, nfu_atoms) in MATS.items():
      rows = {}
      for vname, at in variants(kw).items():
          e, q = d4(at)
          nfu = len(at) / nfu_atoms
          ang = at.cell.cellpar()[3:]
          rows[vname] = {"E_fu_kJ": e / nfu * H2KJ, "q_first": float(q[0]), "natoms": len(at),
                         "angles": [round(float(x), 1) for x in ang]}
      ref = rows["Niggli-reduced"]
      for r in rows.values():
          r["dE_kJ"] = r["E_fu_kJ"] - ref["E_fu_kJ"]; r["dq"] = r["q_first"] - ref["q_first"]
      res[name] = rows
      print(f"\n{name}:  E_D4(Niggli) = {ref['E_fu_kJ']:.3f} kJ/mol per f.u., q = {ref['q_first']:+.4f}")
      for v, r in rows.items():
          print(f"   {v:24s} angles {r['angles']}  dE = {r['dE_kJ']:+8.4f} kJ/mol/f.u. ({100*r['dE_kJ']/abs(ref['E_fu_kJ']):+6.2f} %)  dq = {r['dq']:+.4f}")
  json.dump(res, open("ewald_cells.json", "w"), indent=1)

if __name__ == "__main__": main()
