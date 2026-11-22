**Title:** Ewald sums use a fixed ±2 image range: EEQ charges and D4 energies depend on the choice of an equivalent unit cell

### Summary

For periodic systems, EEQ partial charges change when the same crystal is described with a different but equivalent unit cell. For NaCl, cells related by `a' = a + k b` give q(Na) up to 18 % larger than the primitive or conventional cell, and the D4 dispersion energy changes by up to 1.25 %. Coordination numbers are unaffected (1e-15), so the discrepancy comes from the charge model.

The cause appears to be `get_ewald_trans` in `src/multicharge/ewald/cutoff.f90`. It builds the direct and reciprocal translations from a fixed number of images per lattice vector (`rep = [2, 2, 2]`) instead of a cutoff radius. That is enough for compact cells, but for skewed (non-reduced) cells it misses lattice translations that lie well inside the needed radius.

### Reproducer (dftd4 4.3.1 from PyPI; `rep(3) = [2, 2, 2]` is still on multicharge `main`, 4bb6c96)

```python
import numpy as np
from dftd4.interface import DispersionModel, DampingParam

a = 10.66                                     # NaCl rocksalt, bohr
prim = 0.5 * a * np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]])
num = np.array([11, 17]); pos = np.array([[0, 0, 0], [0.5, 0.5, 0.5]]) @ prim
pbc = np.array([True] * 3)

for k in range(0, 6):
    lat = np.array([[1, k, 0], [0, 1, 0], [0, 0, 1]]) @ prim   # same lattice
    m = DispersionModel(num, pos, lattice=lat, periodic=pbc)
    print(k, m.get_properties()["partial charges"][0],
          m.get_dispersion(DampingParam(method="pbe0"), grad=False)["energy"])
```

| Cell | q(Na) | dq | E_D4 (Eh) | dE/E |
| --- | --- | --- | --- | --- |
| primitive (k = 0) | 0.34807009 | — | −6.88678e−3 | — |
| a′ = a + b | 0.34875674 | +6.9e−4 | −6.88771e−3 | −1.4e−4 |
| a′ = a + 2b | 0.37713721 | +2.9e−2 | −6.92652e−3 | −5.8e−3 |
| a′ = a + 3b | 0.41016439 | +6.2e−2 | −6.97227e−3 | −1.2e−2 |
| a′ = a + 5b | 0.41045396 | +6.2e−2 | −6.97268e−3 | −1.25e−2 |

The conventional cubic cell (8 atoms) gives q(Na) = 0.34806642, consistent with the primitive cell. A triclinic test cell with C, N, O, H and S shows the same effect (max |dq| = 1.1e−3 for a single unimodular transform).

### Verification of the cause

Reducing each skewed cell with LLL (a basis of the same lattice) before calling dftd4 restores q(Na) = 0.34807009 to within 2e−8 for every k. The plateau for k ≥ 3 is consistent with translations that require more than two images along the skewed vector being missing.

### Practical relevance

Non-reduced cells are common: monoclinic and triclinic structures as distributed, cells sheared during NPT molecular dynamics, supercells built with transformation matrices (phonopy, ASE), and tilted slab cells. In those cases D4 results depend on an arbitrary choice of cell, which also breaks reproducibility between users. Forces and stresses inherit the charge error through the D4 charge scaling.

### Possible fixes

1. Generate the Ewald translations from a cutoff radius (direct and reciprocal) derived from `alpha` and the tolerance, as `get_lattice_points(lattice, cutoff, ...)` does elsewhere, instead of a fixed `rep`.
2. Alternatively, reduce the lattice internally (LLL or Niggli) before the Ewald sums. Results do not change for cells that are already reduced.

---
*Found during a systematic, AI-assisted audit (Claude) of numerical libraries used in DFT codes. Verified with the reproducer above and by reading the source.*
