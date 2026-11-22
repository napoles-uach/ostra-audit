# Registro de búsqueda (`searched.md`) — Exploración 3

Registro sistemático de qué patrones se buscaron, en qué librerías, con qué método (lectura estática + compilación/ejecución de oráculos) y qué búsquedas **no produjeron hallazgos** o fueron descartadas por el protocolo de 6 pasos.

---

## 1. Librerías auditadas y versiones/commits exactos

| Librería | Lenguaje | Versión / Commit | Fecha commit | Fuente analizada |
|---|---|---|---|---|
| `libecpint` | C++17 | `v1.1` (`d7966f196a5276ac0203c8ff0cca547d840ef51d`) | 2026-11-18 | `include/libecpint/`, `src/lib/`, `tests/units/` |
| `simple-dftd3` | Fortran 2008 | `v1.7.0` (`84cb7cddd8d6669c4ecd59b06d53e0d373656239`) | 2026-10-28 | `src/dftd3/` (incl. `gcp.f90`, `gcp/param.f90`, `api.f90`, `disp.f90`, `accuracy.f90`) |
| `dftd4` | Fortran 2008 | `HEAD` (`4a5bbf0ffbef983f6b225a4a275e61ea265088f6`) | 2026-11-16 | `src/dftd4/` (incl. `model/d4.f90`, `model/d4s.f90`, `disp.f90`) |
| `multicharge` | Fortran 2008 | `HEAD` (`4bb6c9637508bd1d21a28f72fa514fcf4e968fa0`) | 2026-11-17 | `src/multicharge/` (incl. `charge.f90`, `model/eeq.f90`, `model/eeqbc.f90`, `model/type.F90`, `param/`) |
| `spglib` | C | `HEAD` (`5245d5d315fdcaf79dcb2a5468fecb650fdc177f`) | 2026-11-17 | `include/spglib.h`, `src/*.c` (26 archivos C) |
| `PySCF` (`pyscf/lib`) | C | `2.14.0` instalado + `master` (`nr_ecp.c`, `nr_ecp.h`, 2026-11-22) | 2026-11-22 | `gto/`, `dft/`, `pbc/` (123 archivos C/H excl. `deps/`) |

---

## 2. Matriz de búsqueda (Patrón × Librería × Método)

### Patrón 1: Tablas de tamaño fijo desacopladas del índice máximo

- **Método:** minería estática de arreglos de tamaño constante y `parameter`/`#define` de límites (`max_elem`, `NPRIMAX`, `CART_CUM`, `K_TAB_ENTRIES`, `shell(...)`, `emiss_hf_*`), seguida de rastreo de variables índice, verificación contra `master`, y ejecución con `gfortran -fcheck=bounds` / `clang -fsanitize=address`.
- **Resultados positivos confirmados:**
  - `simple-dftd3`: `Cand-01` (`gcp/param.f90` tablas de 35/36 entradas vs `effective_atomic_number` devolviendo `36..39` y `95..118`; confirmado con `-O2` y `-fcheck=bounds`).
  - `PySCF libcgto` (`master`): `Cand-02` (`nr_ecp.c` `z > 16` $\to$ `entry = 400` OOB en tablas de 400 filas, confirmado con ASan; `omega_nuc[CART_CUM]` con `CART_CUM = 456` $\to$ `SIGABRT` en `l = 6` `ECPscalar_ipnuc`).
  - `dftd4 × multicharge`: `Cand-03` (`dftd4` `max_elem = 118` permite `Z = 112..118` y llama a `multicharge` `max_elem = 103`, que devuelve centinela `-1.0` sin error).
- **Resultados negativos / descartados / ya corregidos:**
  - `PySCF` `nr_ecp.c` `_factorial[24]`: **confirmado en el release `2.14.0`** (ASan `global-buffer-overflow` y error `8.54e-5` en `ECPscalar_cart` para `l = 6`), **pero ya corregido en `master`** (ampliado a 65 entradas). Descartado como defecto nuevo de `master` y documentado solo como contexto en `Cand-02`.
  - `PySCF` `grid_ao_drv.h:33` `NPRIMAX = 40` (`deriv1.c:235`, `pbc/grid_ao.c:81`): probado en Python con una capa `s` de 45 primitivas (`mol.eval_gto("GTOval_sph_deriv1")`); no provocó fallo ni discrepancia frente a FD (`1.74e-11`), aunque en `pbc/grid_ao.c:81` (`double logcoeff[NPRIMAX]`) la inspección estática muestra ausencia de guarda si se pasan >40 primitivas en una sola capa PBC (no explotado con ASan en PBC; marcado como pendiente menor, no promovido a hallazgo).
  - `PySCF` `nr_ecp.c:4593` `_angular_moment_matrix[5]` vs `ECP_LMAX = 5`: probado `mol.intor("ECPso")` con `lc = 5`; devolvió `0.0` limpiamente sin OOB (descartado).
  - `libecpint`: barrido de tablas en `bessel.cpp`, `angular.cpp`, `gshell.cpp`, `mathutil.cpp`, `radial.cpp`. `FAC`/`DFAC` y tablas de `BesselFunction` tienen guardas o tamaños consistentes con `LIBECPINT_MAX_L`; sin OOB de Patrón 1 encontrado.
  - `spglib` `src/spacegroup.c:1240` (`change_of_centering_monocli[36]` indexado con `change_of_basis_index` hasta `47` para grupos de capa monoclínicos `-11 <= hall_number <= -3`): verificado que el acceso solo se ejecuta si `centering == C_FACE`, y los grupos de capa `-11..-3` (monoclínicos/oblicuos) son todos primitivos (`PRIMITIVE`), por lo que el índice `36..47` no alcanza `change_of_centering_monocli[36]`. Descartado.

### Patrón 2: Comentarios o documentación que afirman un hecho numérico no verificado

- **Método:** extracción automática de comentarios que mencionan cifras de precisión, términos suficientes, cotas superiores o equivalencias aproximadas (`1e-`, `accurate`, `enough`, `upper bound`, `machine`, `sufficient`, `~`, `approx`), seguida de contraste analítico con `mpmath` (50 dígitos) y FD4 con barrido de paso.
- **Resultados positivos confirmados:**
  - `libecpint`: `Cand-04` (`BesselFunction::upper_bound` documentada como cota superior de $K_L(z)$ devuelve una cota inferior en $(0, z_{\max}(L))$ para todo $L \ge 1$ y subestima hasta `1,958×` en `L = 24`).
  - `simple-dftd3`: `Cand-05` (`gcp.f90:1109` `"12 terms seem accurate enough"` + corte `|x| < 1e-6` + conmutación por `|za - zb| < 0.1` en vez de `|bx|`, provocando cancelación en `baux_all`/`g3s3s`).
  - `PySCF master` (nota en `Cand-02`): nuevo comentario `nr_ecp.c:4464` afirmando soporte hasta `order = 32` mediante la serie alternante $R_l(-z)/(2z)$ (cancela con error `4.58e9×` en `l = 32, z = 16.01`).
- **Resultados negativos / descartados (verificados limpios):**
  - `simple-dftd3/src/dftd3/accuracy.f90:47, 108-116`: factores de superficie esférica (`2*pi/3`, `pi/4`, `1/5`); verificados algebraicamente correctos.
  - `libecpint/src/lib/screening.cpp:76-77`: cotas `logS` y `Tc`; verificadas algebraicamente válidas.
  - `libecpint/src/lib/bessel.cpp:203` (`"5 terms is usually sufficient for machine accuracy"`): probado por `BesselAccuracy` y nuestra malla con `mpmath`; con `N = 1600` (`dz <= 0.005`), Taylor de grado 5 (`TAYLOR_CUT = 5`) es efectivamente preciso a `~1e-15` en `[0, 16]`. Limpio.
  - `PySCF/lib/pbc/cell.c:35-39` (`pgf_rcut`): `rmin = 2*r_peak` está más allá del máximo de $r^l e^{-\alpha r^2}$; verificado seguro.
  - `PySCF` `ECPscalar_ipnuc` vs FD en `l = 4, 5` (`alpha = 4.0`): inicialmente mostró `1.8%` de discrepancia en `h = 1e-5`, pero el **Paso 2 (barrido de paso `h = 1e-3, 1e-4, 1e-5, 1e-6`)** demostró que el error crecía como $O(1/h)$ por saltos de nivel en la cuadratura adaptativa de `ECPscalar`, no por error en el núcleo. **Descartado por el Paso 2.**

### Patrón 3: Contratos de API implícitos (salidas no escritas / postcondiciones rotas)

- **Método:** búsqueda sistemática de parámetros `optional` combinados con `.and.` en Fortran, `intent(out)` en ramas condicionales, funciones C/C++ con salidas por referencia o puntero que tienen retornos tempranos (`if (...) return`) o ramas `if/else` que solo escriben un subconjunto de salidas.
- **Resultados positivos confirmados:**
  - `libecpint`: `Cand-06` (`BesselFunction::calculate(z <= 0)` deja `values[1..maxL]` sin escribir en la sobrecarga vectorial y devuelve `1.0` para $L \ge 1$ en la escalar).
  - `simple-dftd3`: `Cand-07` (`get_geometric_counterpoise` y API C `dftd3_get_counterpoise` con `gradient` presente y `sigma` ausente/`NULL` dejan `gradient` sin escribir / en ceros).
  - `multicharge`: `Cand-08` (`get_charges` / `get_eeq_charges` con `dqdr` presente y `dqdL` ausente dejan `intent(out) :: dqdr` sin inicializar; ídem `solve` con `gradient` sin `sigma`).
- **Resultados negativos / menores (no promovidos a hallazgo principal):**
  - `simple-dftd3/src/dftd3/disp.f90:102`: verificado que usa `grad = present(gradient) .or. present(sigma)` correctamente (limpio).
  - `dftd4/src/dftd4/disp.f90`: verificado que usa `.or.` para `gradient`/`sigma` (limpio).
  - `spglib/src/spglib.c:983-1007` (`spg_niggli_reduce`): usa copia local `vals[9]` y solo sobrescribe `lattice` si tiene éxito (limpio).
  - `spglib/src/spglib.c:485-509` (`spgms_get_symmetry_with_site_tensors`): **observación menor verificada** — cuando `magnetic_symmetry->size > max_size` (retorna `0` con `SPGERR_ARRAY_SIZE_SHORTAGE`), `equivalent_atoms` y `primitive_lattice` ya fueron sobrescritos mientras `rotation`/`translation`/`spin_flips` quedan sin escribir, y además se fugan `cell` (línea 490 y 501) y `permutations` (línea 501); también el comentario de la línea 507 está invertido respecto al código `1 - 2 * timerev[i]`. No se promovió a hallazgo principal porque el código de error es distinto de cero (no es silencioso).

---

## 3. Observaciones secundarias en `spglib` (fuera de los 3 patrones prioritarios)

1. **Caja de búsqueda fija `KPT_NUM_BZ_SEARCH_SPACE = 125` (`{-2..2}^3`) en `spglib/src/kpoint.c:14-41, 554-620` (`spg_relocate_BZ_grid_address`):**
   Análogo conceptual a `multicharge #72` (suma de Ewald fija en `{-2..2}^3`). En una red cúbica simple sometida a cizalla unimodular entera $U = \begin{pmatrix}1 & s & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1\end{pmatrix}$ (red físicamente idéntica, primera zona de Brillouin con $\max |\mathbf{q}|^2 = 0.75$), ejecutado en `spglib` `HEAD` (`5245d5d3`):
   - `s = 1..5`: `0/64` puntos fuera de la 1.ª BZ (`max |q|^2 = 0.750000`).
   - `s = 6`: **`12/64` puntos fuera de la 1.ª BZ** (`max |q|^2 = 1.500000`).
   - `s = 10`: **`22/64` puntos fuera de la 1.ª BZ** (`max |q|^2 = 9.500000`).
   Documentado aquí porque corresponde al patrón ya conocido de celdas no reducidas (como `multicharge #72`), no a los tres patrones prioritarios de un solo caso.
2. **UBSan en `spglib` (`src/overlap.c` / `src/mathfunc.c:29`):** al compilar `spglib` con `-fsanitize=undefined`, se reportan accesos `double` desalineados a 4 bytes (`0x...7c`) en el búfer empaquetado de `overlap.c` (sin impacto en resultados en x86_64/arm64, sí UB estricto en C).
