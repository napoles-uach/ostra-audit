# Cand-01: Lectura fuera de límites en tablas de tamaño fijo de `simple-dftd3` gCP (`effective_atomic_number` vs arreglos de 35/36 entradas)

- **Librería:** `simple-dftd3` (`v1.7.0`, commit `84cb7cddd8d6669c4ecd59b06d53e0d373656239`, 2026-10-28)
- **Archivo fuente:** `src/dftd3/gcp/param.f90` (líneas 121–127, 487–489, 548–929, 1051–1067) y `src/dftd3/gcp.f90` (líneas 770, 1164, 1281)
- **Patrón del catálogo:** **Patrón 1 — Tablas de tamaño fijo desacopladas del índice máximo** (segundo caso independiente respecto a `OSTRA-008` en `libcint`)
- **Oráculo:** Contrato de límites de arreglo (`gfortran -fcheck=bounds`) + positividad/finitud de energía y gradiente gCP (`E`, `g`)
- **Clase:** Error numérico nuevo / acceso fuera de límites (OOB)
- **Nivel de evidencia:** **Defecto confirmado** (reproducido en binario sin parchear con `-O2` y `-fcheck=bounds`, aislado a línea exacta y eliminado con parche mínimo en copia aislada)
- **Estado upstream:** Presente en `v1.7.0` / `HEAD` (`84cb7cdd`); no encontrado en el tracker de `dftd3/simple-dftd3` ni en `ostra-audit` (`findings/`, `candidates/`, `discarded/`). Sin reportar (Paso 6).

---

## 1. Detección (Paso 1)

En `simple-dftd3`, la corrección geométrica *counterpoise* (`get_gcp_param`) mapea el número atómico `Z` a un índice efectivo `param%zeff` mediante `effective_atomic_number(Z)` (`src/dftd3/gcp/param.f90:1051-1067`) y utiliza ese índice directamente para indexar todas las tablas de parámetros (`emiss_hf_*`, `nbas_*`, `slater_exp`) sin comprobación de rango ni parámetro de error.

Todas las tablas de parámetros de gCP tienen longitud fija **36** (`H`–`Kr`), con una excepción aún más corta: `emiss_hf_minis(*)` (`param.f90:121-127`) tiene solo **35** elementos (falta `Kr = 36`). Sin embargo, `effective_atomic_number` devuelve:
1. **`zeff = 37, 38, 39`** para **`Z = 87` (Fr), `88` (Ra) y `89` (Ac)** — explícitamente incluidos dentro de `case(72:89)` (`zeff = number - 50`)!
2. **`zeff = 36`** para **`Z = 36` (Kr), `54` (Xe) y `86` (Rn)**, que excede la longitud **35** de `emiss_hf_minis`.
3. **`zeff = 95..118`** para **`Z >= 95` (Am–Og)** vía `case default: zeff = number`, excediendo tanto las tablas de longitud 36 en `param.f90` como `integer :: shell(72)` en `ssovl`, `gsovl` y `dsovl` (`src/dftd3/gcp.f90:770, 1164, 1281`).

### Evidencia observada (sin parchear)

Compilando los módulos originales sin modificar de `simple-dftd3` (`cutoff.f90`, `partition.f90`, `gcp/param.f90`, `gcp.f90`) con `gfortran 16.2.0`:

- **En compilación de producción (`-O2`, sin `-fcheck=bounds`)** (`evidence/unpatched_O2.txt`):
  - `Fr-O` (`Z = 87, 8`, `r2SCAN-3c`, `R = 4.0 Bohr`): `zeff = 37` lee fuera de `emiss_hf_def2mtzvpp(36)` (toma `0.00793`, valor de H en `emiss_hf_def2mtzvp`) y fuera de `slater_exp(36)` (toma `0.0` de `slater_d`), devolviendo **`E = +Infinity` y `grad_z = NaN`**.
  - `Ra-O` (`Z = 88, 8`, `r2SCAN-3c`): `zeff = 38` lee `emiss = 0.00331` y `slater = 0.0`, devolviendo **`E = +Infinity` y `grad_z = NaN`**.
  - `Ac-O` (`Z = 89, 8`, `r2SCAN-3c`): `zeff = 39` lee `emiss = 0.08676` y `slater = 0.0`, devolviendo **`E = +Infinity` y `grad_z = NaN`**.
  - `Kr-O` (`Z = 36, 8`, `HF/MINIS`): **error completamente silencioso** — `zeff = 36` lee pasado el final de `emiss_hf_minis(35)` y toma `emiss_hf_631gd(1) = 0.010083` (el parámetro de **Hidrógeno en 6-31G\***!) en lugar de `0.0`, devolviendo una energía gCP espuria **`E = +2.991646e-5 Eh`** (lo mismo ocurre para `Xe` y `Rn`, que también mapean a `zeff = 36`).
- **Con verificación de límites (`-fcheck=bounds`)** (`evidence/unpatched_bounds_check.txt`):
  - `Fr`, `Ra`, `Ac`, `Am` abortan en `src/dftd3/gcp/param.f90:899`: `Index '37'/'38'/'39'/'95' of dimension 1 of array 'emiss_hf_def2mtzvpp' outside of expected range (1:36)`.
  - `Kr` aborta en `src/dftd3/gcp/param.f90:652`: `Index '36' of dimension 1 of array 'emiss_hf_minis' outside of expected range (1:35)`.

## 2. Filtrado de artefactos (Paso 2)

- No depende de pasos de diferencias finitas ni de tolerancias numéricas: es una violación determinista de límites de arreglos Fortran verificada directamente por `gfortran -fcheck=bounds`.
- El shim `mctc_shim.f90` solo provee `wp`, `structure_type` y `get_vdw_rad`; todo el código de `cutoff.f90`, `partition.f90`, `gcp/param.f90` y `gcp.f90` se compiló **íntegro y sin modificaciones** desde el repositorio clonado.

## 3. Cruce con documentación y código fuente (Paso 3)

- `effective_atomic_number` (`param.f90:1051-1067`) mapea deliberadamente los periodos 5 y 6 y los actínidos iniciales (`Th–Pu`, `90:94 -> 21`) a los primeros 36 elementos (`H–Kr`), siguiendo el esquema del método gCP original (Kruse & Grimme, *J. Chem. Phys.* **136**, 154101 (2012)).
- El rango `case(72:89)` fue escrito como si `89 - 50` cupiera en `1..36`, pasando por alto que el periodo 6 termina en `Rn (Z = 86, zeff = 36)` y que `Fr (87), Ra (88), Ac (89)` ya pertenecen al periodo 7. Análogamente, en `emiss_hf_minis` la última fila (`param.f90:127`) contiene solo 5 ceros en vez de 6 (faltando `Kr = 36`).

## 4. Aislamiento de la causa y parche mínimo (Paso 4)

- **Líneas exactas:**
  - `src/dftd3/gcp/param.f90:121-127` (`emiss_hf_minis` tiene 35 entradas en lugar de 36).
  - `src/dftd3/gcp/param.f90:1062-1066` (`case(72:89)` produce `37..39` para `87..89`; `case default` deja pasar `95..118`).
- **Parche mínimo verificado en copia aislada (`repro/minimal.patch`):**
  1. Completar la entrada 36 (`0.0_wp` para Kr) en `emiss_hf_minis`.
  2. Ajustar `effective_atomic_number`: `case(72:86) -> number - 50` (`22..36`), `case(87:88) -> number - 68` (`19..20`, análogo a `Cs, Ba -> 19, 20`), `case(89:103) -> 21` (análogo a `La–Lu, Th–Pu -> 21`), y acotar `case(104:) -> 36` (o idealmente devolver un error explícito si `get_gcp_param` incorporara `error_type`).
- **Resultado con parche (`evidence/patched_bounds_check.txt`, compilado con `-fcheck=bounds`):**
  - Cero errores de límites para `Fr`, `Ra`, `Ac`, `Kr`, `Am`.
  - `Kr-O` (`HF/MINIS`) pasa de `E = 2.991646e-5 Eh` (lectura de H en 6-31G\*) a `E = 0.000000e+00 Eh`.
  - `Fr-O`, `Ra-O`, `Ac-O` (`r2SCAN-3c`) pasan de `E = +Infinity, grad = NaN` a energías y gradientes finitos y bien comportados (`1.557e-4`, `1.934e-4`, `3.443e-3 Eh`).

## 5. Impacto en observables (Paso 5)

- **Alto en `Fr/Ra/Ac/Am..Og` con cualquier método `-3c` (`r2SCAN-3c`, `PBEh-3c`, `HF-3c`):** produce `Infinity`/`NaN` o valores dependientes de la disposición de memoria (UB).
- **Silencioso en `Kr`, `Xe` y `Rn` con `HF/MINIS`:** introduce un término espurio de `~0.02 kcal/mol` por par (dependiente del vecino en memoria, en GCC 16.2 lee `emiss_hf_631gd(1)` de H) que contamina diferencias de energía de complejos con gases nobles pesados (`Kr`, `Xe`, `Rn`).

## 6. Reproducibilidad

Ver `repro/repro_gcp_oob.f90`, `repro/minimal.patch`, `evidence/unpatched_O2.txt`, `evidence/unpatched_bounds_check.txt` y `evidence/patched_bounds_check.txt`.
