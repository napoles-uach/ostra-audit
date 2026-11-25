# Impact study 2: fixed Ewald image range in EEQ charges (OSTRA-005), D4 dispersion of ionic crystals (2026-11-25)

dftd4 4.3.1 (Python), PBE0 damping, structures from `ase.build.bulk` at experimental lattice constants.
Every cell of a material describes the same infinite crystal, so the D4 energy per formula unit must not change.

## 1. Standard cells (what users normally run)
Primitive, conventional cubic, hexagonal gamma = 120 deg vs 60 deg, Niggli-reduced, 2x2x2 supercells:
identical to < 0.001 kJ/mol per formula unit for NaCl, LiF, CsCl, CaF2, ZnO, GaN, Si. Exception: MgO, the most
compact and most charged lattice: the 2-atom primitive cell is off by 0.046 kJ/mol/f.u. (0.15 %), q(Mg) by 0.004,
relative to the converged value (supercells / cubic cells, `supercell_conv.json`, converged to 1e-5 kJ/mol).

## 2. Non-reduced (sheared) but equivalent cells
a' = a + k b. The error grows quickly with the shear and with charge and compactness:

| Material | a' = a + b (angle 30 deg) | a' = a + 2b (angle 19 deg) |
|---|---|---|
| MgO | -0.54 kJ/mol/f.u. (-1.8 %), dq = +0.043 | -10.2 kJ/mol/f.u. (-34 %), dq = +0.61 |
| CaF2 | -0.005 (-0.02 %) | -0.58 (-2.1 %), dq = +0.11 |
| LiF | -0.004 | -0.21 (-1.9 %) |
| NaCl | -0.002 | -0.10 (-0.6 %) |
| CsCl, ZnO, GaN | <= 0.002 | <= 0.05 (<= 0.1 %) |
| Si (neutral) | 0 | 0 |

## Reading
- In routine use (standard, reduced cells) the effect is negligible (<= 0.05 kJ/mol/f.u.): impact level III-IV.
- It becomes large and silent for non-reduced cells of compact ionic solids: cells produced by variable-cell
  MD or relaxations with shear, or taken unreduced from databases/structure generators. Here the D4 energy of MgO
  changes by a third and the Mg charge nearly doubles, with no warning. Level I for that use case.
- A Niggli/LLL reduction before the Ewald sum (or an image range chosen from the cell geometry) removes it; see
  the finding's repro.
