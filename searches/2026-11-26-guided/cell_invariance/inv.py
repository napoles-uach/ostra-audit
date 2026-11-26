# OSTRA guided search, 2026-11-26. Pattern: convergence parameters fixed for the typical case (OSTRA-005).
# Oracle: invariance of the energy per formula unit under equivalent descriptions of the same crystal.
import sys, json, time, numpy as np
sys.path.insert(0, "/home/claude/impact2")
from cells import MATS, variants
AU = 1 / 0.52917721067; H2KJ = 2625.4996
which = sys.argv[1]
def e_tblite(at, method):
    from tblite.interface import Calculator
    c = Calculator(method, at.numbers, at.positions * AU, lattice=np.array(at.cell) * AU, periodic=np.array([True] * 3))
    c.set("verbosity", 0); c.set("accuracy", 1.0); c.set("max-iter", 300)
    r = c.singlepoint()
    return r.get("energy"), float(r.get("charges")[0])
def e_d3(at):
    from dftd3.interface import RationalDampingParam, DispersionModel
    m = DispersionModel(at.numbers, at.positions * AU, lattice=np.array(at.cell) * AU, periodic=np.array([True] * 3))
    return m.get_dispersion(RationalDampingParam(method="pbe0"), grad=False)["energy"], None
def e_pyscf_ewald(at):
    from pyscf.pbc import gto as pgto
    # nuclear (point-charge) Ewald energy with formal ionic charges via a minimal cell: use energy_nuc of a cell
    # whose atoms carry their true Z (neutralised by the uniform background is NOT what energy_nuc does for ions),
    # so instead use pyscf.pbc.tools ewald on the atomic numbers as given (nuclear repulsion of the cell).
    c = pgto.Cell(); c.atom = [(s, p) for s, p in zip(at.get_chemical_symbols(), at.positions)]
    c.a = np.array(at.cell); c.basis = "gth-szv"; c.pseudo = "gth-pade"; c.unit = "A"; c.verbose = 0; c.build()
    return c.energy_nuc(), None
F = {"gfn2": lambda a: e_tblite(a, "GFN2-xTB"), "gfn1": lambda a: e_tblite(a, "GFN1-xTB"),
     "d3": e_d3, "pyscf_ewald": e_pyscf_ewald}[which]
skip = {"pyscf_ewald": ()}
res = {}
for name, (kw, nfa) in MATS.items():
    rows = {}
    for vname, at in variants(kw).items():
        t = time.time()
        try:
            e, q = F(at); rows[vname] = {"E_fu_kJ": e / (len(at) / nfa) * H2KJ, "q": q, "nat": len(at), "t": round(time.time() - t, 1)}
        except Exception as ex:
            rows[vname] = {"err": f"{type(ex).__name__}: {str(ex)[:120]}"}
    ref = rows.get("Niggli-reduced", {})
    for r in rows.values():
        if "E_fu_kJ" in r and "E_fu_kJ" in ref:
            r["dE"] = r["E_fu_kJ"] - ref["E_fu_kJ"]
            if r["q"] is not None and ref.get("q") is not None: r["dq"] = r["q"] - ref["q"]
    res[name] = rows
    print(f"\n[{which}] {name}: E(Niggli) = {ref.get('E_fu_kJ', float('nan')):.4f} kJ/mol/f.u.")
    for v, r in rows.items():
        if "err" in r: print(f"   {v:24s} {r['err']}")
        else: print(f"   {v:24s} dE = {r.get('dE', float('nan')):+10.4f} kJ/mol/f.u.  dq = {r.get('dq', 0) if r.get('dq') is not None else 0:+.4f}  ({r['nat']} at, {r['t']} s)", flush=True)
json.dump(res, open(f"inv_{which}.json", "w"), indent=1)
