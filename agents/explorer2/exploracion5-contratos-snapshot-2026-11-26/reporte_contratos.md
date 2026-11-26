# Auditoría guiada por catálogo — Patrón `implicit-api-contract` (`tblite`, `dftd4`, `simple-dftd3`, `mctc-lib`)

**Fecha:** 2026-11-26  
**Patrón auditado:** Contratos de API implícitos (`implicit-api-contract`, casos previos OSTRA-004, OSTRA-024, OSTRA-025, OSTRA-026)  
**Ámbito:** Barrido sistemático de combinaciones de argumentos de salida `optional` (Fortran) y punteros opcionales `NULL` (C API) en las ramas `main` actuales de `mctc-lib`, `dftd4`, `simple-dftd3` y `tblite`, excluyendo `multicharge::get_charges` (OSTRA-026) y `simple-dftd3::get_geometric_counterpoise` (OSTRA-025).  
**Inventario exhaustivo por rutina:** [rutinas_revisadas.csv](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/rutinas_revisadas.csv)

---

## 1. Resumen experimental (búsqueda guiada vs. exploración libre)

### Entorno y commits auditados

| Librería | Versión / Commit (`main`) | Compilador y banderas | Pruebas ejecutadas |
|---|---|---|---|
| `mctc-lib` | `0.7.0` (`71a21af5de3b`) | `gfortran 16.2.0` (`-Dbuildtype=debug -Dfortran_args=-fcheck=all`) | $2^2$ en 3 modelos CN + CSR + Wigner–Seitz + FD 4 pts + `meson test` (38/38) |
| `dftd4` | `4.3.1` (`4a5bbf0ffbef`) | `gfortran 16.2.0` + `Apple clang 21.0.0 -fsanitize=address` (`-Dopenmp=false -Dlapack=openblas`) | $2^2$ C API (`get_dispersion`) + $2^4$ C API (`get_properties`) + Fortran `get_dispersion`/`get_dispersion2`/`d4_model`/`d4s_model`/`damping` + FD 4 pts + `meson test` (54 OK, 1 XFAIL, 0 Fail) |
| `simple-dftd3` | `1.7.0` (`84cb7cddd8d6`) | `gfortran 16.2.0` + `Apple clang 21.0.0 -fsanitize=address` (`-Dopenmp=false`) | $2^3$ `get_dispersion` + $2^2$ C API + `d3_model` (`weight_references`, `get_atomic_c6`) + 6 amortiguamientos + ATM + `test/unit/tester` |
| `tblite` | `0.8.0` (`c6f66b93f7b3`) | `gfortran 16.2.0` + `Apple clang 21.0.0 -fsanitize=address` (`-Dopenmp=false -Dlapack=openblas -Dddx=false`) | $2^2$ `xtb_singlepoint` (GFN2, GFN2+ALPB) + FD 4 pts (coord y virial) + 4 contenedores `get_engrad` + `get_effective_qat` + `get_selfenergy` + C API + 7 suites `test/unit/tester` |

> **Aislamiento verificado:** Todos los parches se aplicaron y compilaron en copias aisladas (`repos-patched/`, `builds-patched/`). Los cuatro clones base en `repos/` permanecen limpios sin modificaciones rastreadas (`git status -s` vacío en `mctc-lib`, `dftd4`, `tblite`; solo `?? subprojects/jonquil.wrap` en `simple-dftd3` por descarga de subproyectos).

### Métricas globales del barrido guiado

| Librería | Rutinas revisadas | Filas positivas | Filas negativas | Propuestas nuevas principales (causas raíz, Nivel III) | Tiempo aprox. (min) |
|---|---:|---:|---:|---:|---:|
| `mctc-lib` | 8 | 3 | 5 | 1 (`ncoord_type%get_coordination_number`, propagado a `get_cn`/`get_cn_legacy` y a wrappers en `dftd4`/`simple-dftd3`) | 9.0 |
| `dftd4` | 18 | 13 | 5 | 2 (`get_dispersion`/`get_dispersion2`/`dftd4_get_dispersion` **SIGBUS/ASan**; `d4_model`/`d4s_model` 6 métodos) | 22.5 |
| `simple-dftd3` | 13 | 7* | 6 | 1 (`d3_model%weight_references` cero silencioso + `d3_model%get_atomic_c6` no escrito y **SIGBUS**) | 16.5 |
| `tblite` | 17 | 7 | 10 | 1 (`xtb_singlepoint` + 4 contenedores `get_engrad` + `get_effective_qat`) | 21.0 |
| **Total** | **56** | **30** | **26** | **5 causas raíz principales confirmadas (todas parcheadas y verificadas)** | **69.0** |

\*Incluye `get_counterpoise_mpi` (`src/dftd3/mpi.f90:212`), que es la contraparte MPI exacta de OSTRA-025 y **no** se cuenta como hallazgo nuevo.

### Síntesis de resultados y observaciones para el artículo

1. **Alto rendimiento de la búsqueda guiada por catálogo:** En ~69 minutos de auditoría enfocada sobre un único patrón (`implicit-api-contract`), se inspeccionaron 56 rutinas públicas (Fortran y C) en 4 librerías de la misma familia (`grimme-lab`), encontrando **5 causas raíz nuevas independientes con evidencia máxima (*defecto confirmado*)** —incluyendo **dos caídas de memoria (`SIGBUS` / `AddressSanitizer: BUS WRITE`)**, una de ellas directamente en la **API pública de C** (`dftd4_get_dispersion(err, mol, disp, param, &e, NULL, sigma)`) y otra en `d3_model%get_atomic_c6`, además de silencios `intent(out)` en las entradas principales de Fortran (`tblite::xtb_singlepoint`, `mctc-lib::get_coordination_number`, `dftd4::d4_model/d4s_model`, `simple-dftd3::d3_model`).
2. **Contraste entre librerías hermanas (valor de los negativos):** `simple-dftd3` ya corrigió el patrón en su rutina principal `get_dispersion` (Fortran y C API: las 8 máscaras de `grad/sigma/hessian` y las 4 máscaras `NULL` de C pasan con `0.000E+00` gracias a `gradient_local`/`sigma_local`), mientras que `dftd4` (código hermano de la misma familia) usó `grad = present(gradient).or.present(sigma)` sin buffers locales, provocando `SIGBUS` al omitir cualquiera de los dos. De modo análogo, `tblite_get_singlepoint` en la API de C evita el defecto de `xtb_singlepoint` porque siempre asigna ambos arreglos (`res%gradient` y `res%sigma`), mientras que un usuario directo de Fortran cae en el silencio.
3. **Verificación de identidad bit a bit (0 ULP):** Todos los parches propuestos pasan los arneses $2^k$ con error `0.00000E+00`, pasan las suites oficiales de pruebas y preservan **identidad exacta bit a bit (0 diferencias en 188 registros hexadecimales IEEE-754)** en las llamadas completas `(T, T)` y `(F, F)` ya existentes ([bit_identity_summary.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/bit_identity_summary.log)).

### Resumen de alcance por interfaz (transparencia de exposición)

| Propuesta | Librería / Rutina | Modo de fallo | Nivel de evidencia | Nivel de impacto | Fortran API | C API | Python | CLI |
|---|---|---|---|---|---|---|---|---|
| 1 | `dftd4`: `get_dispersion`, `get_dispersion2`, `dftd4_get_dispersion` | **Crash (`SIGBUS` / ASan `BUS WRITE`)** | Defecto confirmado | **III** | **Expuesta** (`grad` o `sigma` solo) | **Expuesta** (`grad=NULL, sigma!=NULL`) | No (pasa ambos cuando `grad=True`) | No (pasa ambos) |
| 2 | `tblite`: `xtb_singlepoint` (+ 4 contenedores `get_engrad`, `get_effective_qat`) | Silencio: `intent(out)` sin escribir (`999.0`) / acumulación sobre no inicializado con CDS / `0.0` | Defecto confirmado | **III** | **Expuesta** (`grad` o `sigma` solo) | No (`api/calculator.f90:684` asigna ambos) | No (envuelve C API) | No (`driver_run.f90:316` pasa ambos) |
| 3 | `mctc-lib`: `ncoord_type%get_coordination_number` (`get_cn`, `get_cn_legacy`; wrappers en `dftd4`/`s-dftd3`) | Silencio: `dcndr`/`dcndL`/`dcndrlist` sin escribir (`999.0`); con `cut>0`, UB al escalar memoria no inicializada | Defecto confirmado | **III** | **Expuesta** (también vía `dftd4`/`s-dftd3`) | No | No | No |
| 4 | `simple-dftd3`: `d3_model%weight_references` y `d3_model%get_atomic_c6` | Cero silencioso (`gwd2cn` solo); `intent(out)` sin escribir; **Crash (`SIGBUS`)** en segundas sin primeras | Defecto confirmado | **III** | **Expuesta** (`dftd3` público) | No (`get_dispersion` pasa todo) | No | No |
| 5 | `dftd4`: `d4_model` y `d4s_model` (`weight_references`, `get_atomic_c6`, `get_polarizabilities`) | Silencio: `intent(out)` sin escribir (`999.0`) al pedir derivadas solo de CN o solo de carga `q` | Defecto confirmado | **III** | **Expuesta** (`!DEC$ ATTRIBUTES DLLEXPORT`) | No | No | No |

---

## 2. Propuestas de hallazgos principales (formato `finding.md`, sin asignar número OSTRA)

---

### OSTRA-NNN (Propuesta 1): `dftd4` — `get_dispersion`, `get_dispersion2` y C API `dftd4_get_dispersion` provocan `SIGBUS` / escritura fuera de límites (ASan `BUS WRITE`) al solicitar solo `gradient` o solo `sigma`

| Field | Value |
|---|---|
| library | dftd4 (Fortran API `dftd4_disp` y C API `dftd4_api`) |
| location | `src/dftd4/disp.f90` (`get_dispersion` L89, L126–127, L135–136, L149–150; `get_dispersion2` L207, L244–245, L253–257); `src/dftd4/api.f90` (`get_dispersion_api` / `dftd4_get_dispersion` L726–744) |
| versions | `4.3.1` (`main` commit `4a5bbf0ffbef983f6b225a4a275e61ea265088f6`) |
| detected | 2026-11-26 |
| detected_by | Antigravity (búsqueda guiada `implicit-api-contract`, ronda 5); pendiente de verificación independiente |
| oracle | Contrato de API (combinaciones $2^2$ de salidas opcionales `gradient`/`sigma` y punteros `NULL` en C bajo AddressSanitizer) + consistencia energía–gradiente/virial (diferencias finitas centrales de 4 puntos) |
| libraries | dftd4 |
| pattern | implicit-api-contract |
| oracle_class | api-contract, sanitizer-code-reading, energy-potential-consistency |
| class | memory safety / implicit API contract, new |
| impact_level | III (caída inmediata por escritura sobre argumento ausente / puntero nulo en la API pública de C y de Fortran al pedir solo virial en C o solo gradiente/virial en Fortran) |
| status | candidate (evidencia local: **defecto confirmado** — anomalía reproducida, causa aislada y parche verificado) |
| report | Pending |

#### Observation
Se ejecutaron las $2^2 = 4$ combinaciones de `(gradient, sigma)` en la API de C ([probe_dftd4_c.c](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/scripts/probe_dftd4_c.c), compilado con `clang -fsanitize=address`) y en la API de Fortran ([probe_dftd4_fortran.f90](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/scripts/probe_dftd4_fortran.f90)) para una molécula periódica de H$_2$O (celda triclinica, PBE-D4, `s9 = 1.0`). Salidas crudas: [probe_dftd4_c.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/probe_dftd4_c.log), [probe_dftd4_fortran.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/probe_dftd4_fortran.log).

| Interfaz y llamada | Combinación `(grad, sigma)` | Resultado baseline | Resultado tras parche |
|---|---|---|---|
| C API `dftd4_get_dispersion` (ASan) | `(ptr, ptr)` | `EXIT=0`, `E=-8.54620162e-4`, `max|g|=3.727e-5`, `max|s|=1.213e-3` | `EXIT=0`, idéntico (`0.000e+00`) |
| C API `dftd4_get_dispersion` (ASan) | `(NULL, NULL)` | `EXIT=0`, `|dE|=0.0` | `EXIT=0`, `|dE|=0.0` |
| C API `dftd4_get_dispersion` (ASan) | `(ptr, NULL)` | `EXIT=0`, `max|dg|=0.0` (salvado por parche asimétrico en `api.f90:738`) | `EXIT=0`, `max|dg|=0.0` |
| C API `dftd4_get_dispersion` (ASan) | `(NULL, ptr)` | **`EXIT=134` — `AddressSanitizer: BUS` (`WRITE` en `disp.f90:126` desde `api.f90:744`)** | **`EXIT=0`, ASan limpio, `max|ds|=0.000e+00`** |
| Fortran `get_dispersion` | `(T, F)` solo `gradient` | **`EXIT=138` (`SIGBUS`)** | **`EXIT=0`, `max|g|=4.32252e-5` (coincide con ref)** |
| Fortran `get_dispersion` | `(F, T)` solo `sigma` | **`EXIT=138` (`SIGBUS`)** | **`EXIT=0`, `max|s|=1.21348e-3` (coincide con ref)** |
| Fortran `get_dispersion2` | `(T, F)` solo `gradient` | **`EXIT=138` (`SIGBUS`)** | **`EXIT=0`, `max|g|=4.42974e-5`** |
| Fortran `get_dispersion2` | `(F, T)` solo `sigma` | **`EXIT=138` (`SIGBUS`)** | **`EXIT=0`, `max|s|=1.25604e-3`** |

#### Artifact filter
- Compilación limpia sin optimizaciones agresivas (`-Dbuildtype=debug -Dopenmp=false -Dfortran_args=-fcheck=all`) y con `clang -g -O0 -fsanitize=address -fno-omit-frame-pointer`.
- Oráculo físico de diferencias finitas centrales de 4 puntos en el cálculo de referencia (`h = 1e-2, 1e-3, 1e-4, 1e-5`): para `h = 1e-4`, `max|g_ref - FD| = 4.121e-14` (precisión de máquina); para `sigma` con deformación de celda/coordenadas, `max|s_ref - FD| = 2.570e-6` (convergido con el corte discreto de red).
- El mismo arnés en C con `dftd4_get_properties` ($2^4 = 16$ máscaras `NULL`) y en `simple-dftd3` (`dftd3_get_dispersion`, 4 máscaras `NULL`) pasa completamente limpio bajo ASan, descartando cualquier error del arnés.

#### Documentation and paper
- En `include/dftd4.h`, `dftd4_get_dispersion` declara `double *gradient` y `double *sigma` como punteros de salida opcionales anulables (sin indicar que `sigma` no pueda pasarse con `gradient == NULL`).
- En `src/dftd4/disp.f90:57-61`, `gradient` y `sigma` se declaran `real(wp), intent(out), optional` de manera independiente, y la línea 89 (`grad = present(gradient).or.present(sigma)`) demuestra la intención explícita del diseño de aceptar cualquiera de los dos por separado. El comentario/workaround parcial en `api.f90:738-740` (que solo asignó `sigma` temporal cuando `has_grad` era verdadero) confirma que se trató de un parche incompleto para el caso `(grad, NULL)` que dejó expuesto el caso simétrico `(NULL, sigma)` y toda la API de Fortran. **Defecto inequívoco.**

#### Cause
En `src/dftd4/disp.f90`:
- Línea 89 (`get_dispersion`) y línea 207 (`get_dispersion2`): `grad = present(gradient).or.present(sigma)`.
- Sin embargo, dentro del bloque `if (grad)` (L126–127 y L244–245), el código ejecuta incondicionalmente `gradient(:, :) = 0.0_wp` y `sigma(:, :) = 0.0_wp` y los pasa a `get_dispersion2`/`get_dispersion3` y a las contracciones `gemv` (L135–136, L149–150, L253–257), desreferenciando el argumento opcional ausente (`NULL` en ABI Fortran/C).
- En `src/dftd4/api.f90:738-740`, el envoltorio de C solo creó un arreglo temporal para `sigma` cuando `has_grad` estaba activo (`else if (has_grad) then allocate(sigma(3, 3))`), omitiendo el caso simétrico `has_sigma .and. .not. has_grad`.

**Parche mínimo probado** ([dftd4-disp-capi.patch](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/patches/dftd4-disp-capi.patch)):
- En `get_dispersion` y `get_dispersion2`, asignar arreglos locales `gradient_local(3, mol%nat)` y `sigma_local(3, 3)` cuando `grad` es verdadero (exactamente el mismo patrón ya empleado por `simple-dftd3/src/dftd3/disp.f90:102-184`), y copiar al final solo si `present(gradient)` / `present(sigma)`.
- En `api.f90`, eliminar el workaround asimétrico (pasando directamente `gradient` o `sigma` solo cuando el puntero C está asociado).
- **Verificación del parche:**
  - Todas las pruebas en C bajo ASan ([patched_probe_dftd4_c.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/patched_probe_dftd4_c.log)) y Fortran ([patched_probe_dftd4_fortran.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/patched_probe_dftd4_fortran.log)) terminan con `EXIT=0` y error `0.000E+00`.
  - Identidad bit a bit en `(T, T)` y `(F, F)`: **67/67 registros hexadecimales IEEE-754 idénticos (0 diffs)** ([bit_identity_summary.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/bit_identity_summary.log)).
  - Suite oficial `meson test`: **54 OK, 1 Expected Fail, 0 Fail, 11 Skipped** ([patched_test_dftd4.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/patched_test_dftd4.log)).

#### Impact
- **Afectado:** Cualquier llamador de la **API pública de C** (`dftd4.h`) que pida solo el virial (`gradient=NULL, sigma!=NULL`, p. ej. cálculos de estrés/presión a coordenadas fijas), y cualquier llamador de la **API pública de Fortran** que llame a `get_dispersion` o `get_dispersion2` pidiendo solo `gradient` (caso muy habitual en moléculas aisladas sin condiciones periódicas) o solo `sigma`. Produce una caída inmediata (`SIGBUS` / `SIGSEGV`).
- **No afectado:** La CLI `dftd4` (`app/driver.f90`) y el binding de Python (`python/dftd4/interface.py`) siempre asignan y pasan tanto `gradient` como `sigma` cuando `grad=.true.`.

#### Verification
Pendiente de verificación independiente por un segundo agente/persona.

#### Corrections
Ninguna hasta la fecha.

---

### OSTRA-NNN (Propuesta 2): `tblite` — `xtb_singlepoint`, cuatro contenedores públicos `get_engrad` y `get_effective_qat` (CEH) omiten en silencio `gradient`, `sigma` o `dqatdr` si no se solicitan ambos simultáneamente

| Field | Value |
|---|---|
| library | tblite (`xtb_singlepoint`, contenedores `repulsion`, `halogen`, `disp/d3`, `disp/d4`, y `ceh`) |
| location | `src/tblite/xtb/singlepoint.f90` (`xtb_singlepoint` L151); `src/tblite/repulsion/effective.f90` (`get_engrad` L125); `src/tblite/classical/halogen.f90` (`get_engrad` L138); `src/tblite/disp/d3.f90` (`get_engrad` L109); `src/tblite/disp/d4.f90` (`get_engrad` L224); `src/tblite/ceh/ceh.f90` (`get_effective_qat` L1308–1309) |
| versions | `0.8.0` (`main` commit `c6f66b93f7b3ed5db77698021822eaedb8074ce6`) |
| detected | 2026-11-26 |
| detected_by | Antigravity (búsqueda guiada `implicit-api-contract`, ronda 5); pendiente de verificación independiente |
| oracle | Contrato de API (combinaciones $2^2$ con centinela `999.0_wp` y `0.0_wp`) + consistencia energía–gradiente y energía–virial (FD centrales de 4 puntos) |
| libraries | tblite |
| pattern | implicit-api-contract |
| oracle_class | api-contract, energy-potential-consistency |
| class | silent wrong result / uninitialized output, new |
| impact_level | III (un usuario de Fortran que invoque `xtb_singlepoint(..., gradient=g)` en una molécula sin pasar `sigma` recibe la memoria no inicializada o basura previa de `g` sin advertencia ni error; con solvatación CDS además acumula el término CDS sobre memoria no inicializada) |
| status | candidate (evidencia local: **defecto confirmado** — anomalía reproducida, causa aislada y parche verificado) |
| report | Pending |

#### Observation
Se probaron las combinaciones parciales de `gradient`/`sigma` en `xtb_singlepoint` (GFN2-xTB y GFN2+ALPB sobre H$_2$O no periódica), en los 4 contenedores `get_engrad` (`repulsion`, `d4_disp`, `d3_disp`, `halogen` sobre C–Br···N) y en `get_effective_qat` ([probe_tblite_fortran.f90](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/scripts/probe_tblite_fortran.f90), [probe_tblite_virial_fd.f90](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/scripts/probe_tblite_virial_fd.f90)). Salidas crudas: [probe_tblite_fortran.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/probe_tblite_fortran.log), [probe_tblite_virial_fd.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/probe_tblite_virial_fd.log).

| Rutina | Combinación probada | Valor de referencia (`max|ref|`) | Salida baseline (error) | Salida tras parche (`max|err|`) |
|---|---|---|---|---|
| `xtb_singlepoint` (GFN2) | `(T, F)` solo `gradient` (`intent(out)`, centinela `999.0`) | `2.35660E-02` | `min(g)=9.99000E+02` (sin escribir, `max|err|=999.02`) | **`0.00000E+00`** |
| `xtb_singlepoint` (GFN2) | `(F, T)` solo `sigma` (`intent(out)`, centinela `999.0`) | `3.32319E-02` | `min(s)=9.99000E+02` (sin escribir, `max|err|=999.00`) | **`0.00000E+00`** |
| `xtb_singlepoint` (GFN2+ALPB) | `(T, F)` solo `gradient` (centinela `999.0`) | `1.79126E-02` | `min(g)=9.99000E+02` (sin escribir) | **`0.00000E+00`** |
| `repulsion%get_engrad` | `(T, F)` / `(F, T)` (`intent(inout)` desde `0.0`) | `9.01147E-02` / `1.44284E-01` | `0.00000E+00` (silencio total) | **`0.00000E+00`** |
| `d4_disp%get_engrad` | `(T, F)` / `(F, T)` (`intent(inout)` desde `0.0`) | `1.28712E-10` / `1.72719E-10` | `0.00000E+00` (silencio total) | **`0.00000E+00`** |
| `d3_disp%get_engrad` | `(T, F)` / `(F, T)` (`intent(inout)` desde `0.0`) | `7.18245E-07` / `1.19873E-06` | `0.00000E+00` (silencio total) | **`0.00000E+00`** |
| `halogen%get_engrad` | `(T, F)` / `(F, T)` (`intent(inout)` desde `0.0`) | `5.04778E-03` / `2.72580E-02` | `0.00000E+00` (silencio total) | **`0.00000E+00`** |
| `get_effective_qat` (CEH) | `(dqatdr, dcn_endr)` sin `dqatdL`/`dcn_endL` | `9.58777E-02` | `0.00000E+00` (puesto a cero y saltado) | **`0.00000E+00`** |

Además, por inspección del flujo en `xtb_singlepoint` (L151–159, L209) y `solvation/cds.f90:286-293`: como `cds_solvation%get_engrad` sí comprueba `if (present(gradient))` independientemente, al llamar `xtb_singlepoint(..., gradient=g)` con solvatación CDS sin `sigma`, `xtb_singlepoint` salta `gradient(:,:) = 0.0_wp` y `cds%get_engrad` ejecuta `gemv(..., gradient, beta=1.0_wp)` acumulando el término CDS directamente sobre memoria no inicializada.

#### Artifact filter
- Oráculo físico de diferencias finitas centrales de 4 puntos para GFN2-xTB:
  - Gradiente (`probe_tblite_fortran.log`): `max|g_ref - FD| = 1.483e-7` (`h=1e-3, 1e-4`), `1.481e-7` (`h=1e-5`) — limitado únicamente por la convergencia SCF (`accuracy=1e-4`).
  - Virial (`probe_tblite_virial_fd.log`, deformación $\mathbf{r} \to (\mathbf{I} + \boldsymbol{\epsilon})\mathbf{r}$): `max|s_ref - FD| = 2.311e-11` (`h=1e-3`), `3.608e-11` (`h=1e-4`), `1.793e-10` (`h=1e-5`).
- Control positivo del mismo código (`wignerseitz_cell%get_wignerseitz_weights`): devuelve `max|err| = 0.000E+00` con `dwdr` solo.

#### Documentation and paper
Ni las firmas públicas ni los docstrings de `xtb_singlepoint` (`src/tblite/xtb/singlepoint.f90:75-104`), `container_type%get_engrad` (`src/tblite/container/type.f90:86-98`) ni `get_effective_qat` (`src/tblite/ceh/ceh.f90:1285-1304`) documentan que `gradient` y `sigma` (o `dqatdr` y `dqatdL`) deban pasarse conjuntamente. Además, en `xtb_singlepoint`, `gradient` y `sigma` son `intent(out), optional`, por lo que dejar una variable `intent(out)` presente sin asignar infringe el contrato del estándar Fortran. **Defecto inequívoco.**

#### Cause
- `src/tblite/xtb/singlepoint.f90:151`: `grad = present(gradient) .and. present(sigma)`. Si solo uno está presente, `grad` es `.false.` y se omite tanto la inicialización `= 0.0_wp` como toda la evaluación de derivadas (SCF, `h0`, dispersión, repulsión, halógeno, Coulomb y `alpb`).
- `src/tblite/repulsion/effective.f90:125`, `src/tblite/classical/halogen.f90:138`, `src/tblite/disp/d3.f90:109`, `src/tblite/disp/d4.f90:224`: mismo patrón `present(gradient) .and. present(sigma)` en `get_engrad`.
- `src/tblite/ceh/ceh.f90:1308-1309`: `if (present(dqatdr) .and. present(dqatdL) .and. present(dcn_endr) .and. present(dcn_endL))` tras poner `dqatdr = 0.0_wp` en L1305.

**Parche mínimo probado** ([tblite-singlepoint.patch](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/patches/tblite-singlepoint.patch)):
- En `xtb_singlepoint`: `grad = present(gradient) .or. present(sigma)` con `gradient_local`/`sigma_local`, copiando al final a `gradient`/`sigma` si están presentes.
- En `repulsion/effective.f90` y `classical/halogen.f90`: ramas `else if (present(gradient))` / `else if (present(sigma))` con buffer local solo para el argumento ausente.
- En `disp/d3.f90` y `disp/d4.f90`: `grad = present(gradient) .or. present(sigma)` inicializando `gradient_local` desde `gradient` (si está presente) y `sigma_local` desde `sigma` (si está presente), preservando exactamente el orden de suma en punto flotante.
- En `ceh/ceh.f90`: comprobar independientemente `present(dqatdr) .and. present(dcn_endr)` y `present(dqatdL) .and. present(dcn_endL)`.
- **Verificación del parche:**
  - Todas las combinaciones parciales dan `0.00000E+00` ([patched_probe_tblite_fortran.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/patched_probe_tblite_fortran.log), [probe_tblite_virial_fd.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/probe_tblite_virial_fd.log)).
  - Identidad bit a bit en `(T, T)` y `(F, F)`: **46/46 registros hexadecimales IEEE-754 idénticos (0 diffs)** ([bit_identity_summary.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/bit_identity_summary.log)).
  - Suites unitarias oficiales (`repulsion`, `halogen`, `dispersion`, `gfn1-xtb`, `gfn2-xtb`, `ipea1-xtb`, `ceh`): **todas PASSED (`SUITE_EXIT=0`)** ([patched_test_tblite.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/patched_test_tblite.log)).

#### Impact
- **Afectado:** La **API pública de Fortran** de `tblite` (el caso de uso más natural en moléculas aisladas es `call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, gradient)` omitiendo `sigma`, que devuelve silencio total / basura de memoria en `gradient`, rompiendo cualquier optimización de geometría o dinámica molecular que llame a `tblite` desde Fortran).
- **No afectado:** La **API pública de C** (`tblite_get_singlepoint` en `src/tblite/api/calculator.f90:684` asigna siempre ambos `res%gradient` y `res%sigma`; probado con ASan en [probe_tblite_c.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/probe_tblite_c.log)), el binding de Python y la CLI `tblite` (`app/driver_run.f90:316`).

#### Verification
Pendiente de verificación independiente por un segundo agente/persona.

#### Corrections
Ninguna hasta la fecha.

---

### OSTRA-NNN (Propuesta 3): `mctc-lib` — `ncoord_type%get_coordination_number` (`get_cn`, `get_cn_legacy` y wrappers de `dftd4`/`simple-dftd3`) omite derivadas si se pide solo `dcndr`, `dcndL` o `dcndrlist`, y lee memoria no inicializada cuando `cut > 0`

| Field | Value |
|---|---|
| library | mctc-lib (propagado a `dftd4::dftd4_ncoord` y `simple-dftd3::dftd3_ncoord`) |
| location | `src/mctc/ncoord/type.f90` (`get_coordination_number` L215–231; wrappers `get_cn` L130, `get_cn_legacy` L161) |
| versions | `mctc-lib 0.7.0` (`71a21af5de3b`); propagado en `dftd4 4.3.1` (`4a5bbf0ffbef`, `src/dftd4/ncoord.f90:72`) y `simple-dftd3 1.7.0` (`84cb7cddd8d6`, `src/dftd3/ncoord.f90`) |
| detected | 2026-11-26 |
| detected_by | Antigravity (búsqueda guiada `implicit-api-contract`, ronda 5); pendiente de verificación independiente |
| oracle | Contrato de API ($2^2$ en 3 modelos CN + versión CSR con centinela `999.0_wp`) + derivadas por diferencias finitas centrales de 4 puntos |
| libraries | mctc-lib, dftd4, simple-dftd3 |
| pattern | implicit-api-contract |
| oracle_class | api-contract, energy-potential-consistency |
| class | silent wrong result / uninitialized memory read, new |
| impact_level | III (salidas `intent(out)` sin escribir en la API pública de número de coordinación de `mctc-lib`, `dftd4` y `simple-dftd3`; además lectura de memoria no inicializada en `cut_coordination_number` cuando `cut > 0`) |
| status | candidate (evidencia local: **defecto confirmado** — anomalía reproducida, causa aislada y parche verificado) |
| report | Pending |

#### Observation
Se probaron los modelos `cn_count%exp`, `cn_count%dftd4` y `cn_count%erf` con `cut = 6.0_wp`, más la variante con lista de vecinos CSR (`dcndrlist`, `dcndL`), sobre H$_2$O periódica ([probe_mctc_ncoord.f90](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/scripts/probe_mctc_ncoord.f90)). Salida cruda: [probe_mctc_ncoord.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/probe_mctc_ncoord.log).

| Modelo / Ruta | Combinación probada | `max|ref|` | Baseline (`min..max` con centinela `999.0`) | Tras parche (`max|err|`) |
|---|---|---|---|---|
| `exp` (`get_cn`) | `(T, F)` solo `dcndr` | `9.32731E-02` | `9.99000E+02` (sin escribir) | **`0.00000E+00`** |
| `exp` (`get_cn`) | `(F, T)` solo `dcndL` | `1.70875E-01` | `9.99000E+02` (sin escribir) | **`0.00000E+00`** |
| `dftd4` (`get_cn`) | `(T, F)` / `(F, T)` | `7.02540E-02` / `1.30654E-01` | `9.99000E+02` (sin escribir) | **`0.00000E+00`** |
| `erf` con `cut = 6.0` | `(T, F)` solo `dcndr` | `4.82638E-01` | **`9.85573E+02..9.93238E+02`** (**lee y escala el centinela `999.0`!**) | **`0.00000E+00`** |
| `erf` con `cut = 6.0` | `(F, T)` solo `dcndL` | `8.04738E-01` | **`9.85573E+02..9.93238E+02`** (**lee y escala el centinela `999.0`!**) | **`0.00000E+00`** |
| CSR list (`get_coordination_number`) | `(T, F)` solo `dcndrlist` / `(F, T)` solo `dcndL` | `9.32731E-02` / `1.70875E-01` | `9.99000E+02` (sin escribir) | **`0.00000E+00`** |
| `dftd4_ncoord::get_coordination_number` | `(T, F)` / `(F, T)` | `7.02540E-02` / `1.30654E-01` | `9.99000E+02` (propaga `mctc-lib`) | Resuelto al actualizar `mctc-lib` |
| `dftd3_ncoord::get_coordination_number` | `(T, F)` / `(F, T)` | `9.32731E-02` / `1.70875E-01` | `9.99000E+02` (propaga `mctc-lib`) | Resuelto al actualizar `mctc-lib` |

#### Artifact filter
- Oráculo físico de diferencias finitas centrales de 4 puntos de `cn` respecto a coordenadas y deformación de red (`probe_mctc_ncoord.log`):
  - `exp`: `max|dcndr - FD| = 7.539e-12` (`h=1e-3`), `max|dcndL - FD| = 1.998e-11` (`h=1e-4`).
  - `dftd4`: `max|dcndr - FD| = 8.742e-13`, `max|dcndL - FD| = 1.546e-12` (`h=1e-4`).
  - `erf_cut6`: `max|dcndr - FD| = 8.833e-13`, `max|dcndL - FD| = 1.013e-11` (`h=1e-3`).
- Control interno en `mctc-lib` (`get_wignerseitz_weights`): `(T, F)` y `(F, T)` dan `0.00000E+00`.

#### Documentation and paper
- En `mctc-lib` (`src/mctc/ncoord/type.f90:130-210`), `dcndr`, `dcndL` y `dcndrlist` están declarados y documentados como argumentos `intent(out), optional` sin ninguna restricción que exija pasarlos en pares.
- **Nota honesta de código aguas abajo:** En `dftd4/src/dftd4/disp.f90:311-313` existe un comentario interno en una rutina de `dftd4` que señala *"mctc-lib currently evaluates Cartesian CN derivatives only when both dcndr and dcndL are requested"*, lo que confirma que el comportamiento es una limitación conocida en al menos un cliente interno, **no** un diseño deseado, y sigue sin documentarse ni corregirse en `mctc-lib`, produciendo además comportamiento indefinido (lectura de memoria no inicializada) en `cut_coordination_number` cuando `cut > 0`.

#### Cause
En `src/mctc/ncoord/type.f90`:
- Línea 216: `if (present(dcndrlist) .and. present(dcndL))` (rama con lista CSR).
- Línea 222: `if (present(dcndr) .and. present(dcndL))` (rama estándar).
- Si solo uno está presente, cae en la rama `else` (solo calcula `cn`), dejando el arreglo `intent(out)` solicitado sin escribir.
- Luego, en líneas 229–231, si `self%cut > 0.0_wp`, llama a `cut_coordination_number(self%cut, cn, dcndr, dcndL, dcndrlist)`, la cual **sí** comprueba `present(dcndr)` / `present(dcndL)` / `present(dcndrlist)` individualmente y multiplica `dcndr = dcndr * dcn_cut` leyendo la memoria sin inicializar.

**Parche mínimo probado** ([mctc-lib-ncoord.patch](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/patches/mctc-lib-ncoord.patch)):
- En `get_coordination_number`, usar arreglos temporales `dcndr_local`/`dcndL_local`/`dcndrlist_local` cuando se solicita solo una de las derivadas (`else if (present(dcndr))`, `else if (present(dcndL))`, `else if (present(dcndrlist))`), manteniendo la ruta original sin copias cuando ambos están presentes.
- **Verificación del parche:**
  - Todas las combinaciones (`exp`, `dftd4`, `erf_cut6`, CSR) dan `max|err| = 0.00000E+00` ([patched_probe_mctc_ncoord.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/patched_probe_mctc_ncoord.log)).
  - Identidad bit a bit en `(T, T)` y `(F, F)`: **57/57 registros hexadecimales IEEE-754 idénticos (0 diffs)** ([bit_identity_summary.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/bit_identity_summary.log)).
  - Suite oficial `meson test`: **38/38 OK** ([patched_test_mctc-lib.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/patched_test_mctc-lib.log)).

#### Impact
- **Afectado:** API pública de Fortran de `mctc-lib` (`ncoord%get_cn`, `ncoord%get_cn_legacy`, `ncoord%get_coordination_number`) y los wrappers públicos `dftd4_ncoord::get_coordination_number` y `dftd3_ncoord::get_coordination_number`. Cualquier usuario que pida solo `dcndr` (p. ej. para moléculas aisladas) recibe memoria no inicializada o basura escalada (si `cut > 0`).
- **No afectado:** APIs de C, Python y CLI de `dftd4`, `simple-dftd3` y `tblite` (que siempre pasan ambos arreglos internamente).

#### Verification
Pendiente de verificación independiente por un segundo agente/persona.

#### Corrections
Ninguna hasta la fecha.

---

### OSTRA-NNN (Propuesta 4): `simple-dftd3` — `d3_model%weight_references` devuelve ceros silenciosos para `gwd2cn` solo, y `d3_model%get_atomic_c6` deja derivadas sin escribir o provoca `SIGBUS` si se solicitan segundas derivadas sin las primeras

| Field | Value |
|---|---|
| library | simple-dftd3 (`dftd3_model`) |
| location | `src/dftd3/model.f90` (`weight_references` L245–254; `get_atomic_c6` L312–360) |
| versions | `1.7.0` (`main` commit `84cb7cddd8d6669c4ecd59b06d53e0d373656239`) |
| detected | 2026-11-26 |
| detected_by | Antigravity (búsqueda guiada `implicit-api-contract`, ronda 5); pendiente de verificación independiente |
| oracle | Contrato de API (combinaciones parciales de salidas opcionales con centinela `999.0_wp` + ejecución con `-fcheck=all`) |
| libraries | simple-dftd3 |
| pattern | implicit-api-contract |
| oracle_class | api-contract, sanitizer-code-reading |
| class | silent wrong result + memory safety (`SIGBUS`), new |
| impact_level | III (cero silencioso en `weight_references`, salida sin escribir en `get_atomic_c6` y caída `SIGBUS` en la API pública Fortran de `d3_model`) |
| status | candidate (evidencia local: **defecto confirmado** — anomalía reproducida, causa aislada y parche verificado) |
| report | Pending |

#### Observation
Se probaron combinaciones parciales de argumentos opcionales en `d3_model%weight_references` y `d3_model%get_atomic_c6` sobre H$_2$O ([probe_sdftd3_fortran.f90](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/scripts/probe_sdftd3_fortran.f90)). Salida cruda: [probe_sdftd3_fortran.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/probe_sdftd3_fortran.log).

| Rutina | Combinación probada | Valor ref | Resultado baseline | Resultado tras parche |
|---|---|---|---|---|
| `d3%weight_references` | `(F, T)` solo `gwd2cn` (sin `gwdcn`) | `max|gwd2cn| = 1.12872` | **`0.00000E+00` silencioso** (`max|err| = 1.12872`) | **`0.00000E+00`** (`max|gwd2cn|=1.12872`) |
| `d3%get_atomic_c6` | `gwdcn, dc6dcn, gwd2cn, d2c6dcn2` (sin `d2c6dcnij`) | `min(d2c6)=0.78296` | `d2c6dcn2` sin escribir (`9.99000E+02`) | **`0.00000E+00`** |
| `d3%get_atomic_c6` | `gwdcn, dc6dcn, gwd2cn, d2c6dcnij` (sin `d2c6dcn2`) | `min(d2mix)=4.0089E-03` | `d2c6dcnij` sin escribir (`9.99000E+02`) | **`0.00000E+00`** |
| `d3%get_atomic_c6` | `gwd2cn, d2c6dcn2, d2c6dcnij` **sin** `gwdcn`/`dc6dcn` | `max|d2c6|=...` | **`EXIT=138` (`SIGBUS`) en L316/L342** | **`EXIT=0`, `max|d2c6-ref| = 0.00000E+00`** |

#### Artifact filter
- Compilación debug con `-fcheck=all`, sin OpenMP (`-Dopenmp=false`).
- En el mismo arnés, el control de alto nivel `get_dispersion` ($2^3 = 8$ máscaras de `grad/sigma/hessian`) devuelve `0.000E+00` en todas las combinaciones.

#### Documentation and paper
En `src/dftd3/model.f90:200-220` y `268-303`, todos los argumentos de derivada (`gwdcn`, `gwd2cn`, `dc6dcn`, `d2c6dcn2`, `d2c6dcnij`) se declaran `optional` sin documentar acoplamiento entre sí. Además, `weight_references` inicializa `if (present(gwd2cn)) gwd2cn = 0.0_wp` (L227) fuera del `if (present(gwdcn))`, lo que evidencia intención de tratar `gwd2cn` como opcional, pero luego anida su cálculo dentro de `if (present(gwdcn))` (L245–254). En `get_atomic_c6`, el `SIGBUS` ocurre porque `second = present(gwd2cn) .and. present(d2c6dcn2) .and. present(d2c6dcnij)` (L312) se activa aunque `first` (L311) sea `.false.`, pero el bloque `if (first .or. second)` escribe incondicionalmente `dc6dcn(:,:) = 0.0_wp` (L316) y lee `gwdcn` (L342). **Defecto inequívoco.**

#### Cause
1. `src/dftd3/model.f90:245-254` (`weight_references`): `if (present(gwd2cn))` está anidado dentro de `if (present(gwdcn))`; si `gwdcn` no está presente, `gwd2cn` queda en `0.0_wp`.
2. `src/dftd3/model.f90:311-360` (`get_atomic_c6`):
   - `second` acopla con `.and.` a `d2c6dcn2` y `d2c6dcnij` (aunque `d2c6dcn2` solo depende de `gwd2cn`, y `d2c6dcnij` solo de `gwdcn`!).
   - Al evaluar `if (first .or. second)` cuando `second` es verdadero y `first` falso, L316 (`dc6dcn = 0.0_wp`) y L342 (`gwdcn(...)`) desreferencian argumentos ausentes → `SIGBUS`.

**Parche mínimo probado** ([simple-dftd3-model.patch](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/patches/simple-dftd3-model.patch)):
- En `weight_references`: calcular `dnorm` y `dgwk` cuando `present(gwdcn) .or. present(gwd2cn)`, y escribir `gwdcn` y `gwd2cn` con sus guardas independientes.
- En `get_atomic_c6`: separar las condiciones exactas de dependencia matemática: `has_dc6 = present(gwdcn) .and. present(dc6dcn)`, `has_d2cn2 = present(gwd2cn) .and. present(d2c6dcn2)`, `has_d2cnij = present(gwdcn) .and. present(d2c6dcnij)` (y poner a cero `d2c6dcnij` solo si corresponde; si `d2c6dcnij` se pide sin `gwdcn`, no desreferencia `gwdcn` ni `dc6dcn`).
- **Verificación del parche:**
  - Todas las pruebas parciales y el antiguo caso `SIGBUS` (`c6_second_no_first`) terminan con `EXIT=0` y `max|err| = 0.00000E+00` ([patched_probe_sdftd3_fortran.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/patched_probe_sdftd3_fortran.log)).
  - Identidad bit a bit en llamadas completas: **18/18 registros hexadecimales IEEE-754 idénticos (0 diffs)** ([bit_identity_summary.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/bit_identity_summary.log)).
  - Suite unitaria oficial (`test/unit/tester`): **UNIT_EXIT=0**, todas las suites PASSED ([patched_test_simple-dftd3.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/patched_test_simple-dftd3.log)).

#### Impact
- **Afectado:** Usuarios externos de la API Fortran pública `dftd3` que invoquen `d3_model%weight_references` o `d3_model%get_atomic_c6` directamente (p. ej. al acoplar derivadas D3 o hessianas parciales en otro código).
- **No afectado:** `get_dispersion` de alto nivel (Fortran y C API), Python y CLI de `simple-dftd3`.

#### Verification
Pendiente de verificación independiente por un segundo agente/persona.

#### Corrections
Ninguna hasta la fecha.

---

### OSTRA-NNN (Propuesta 5): `dftd4` — `d4_model` y `d4s_model` (`weight_references`, `get_atomic_c6`, `get_polarizabilities`) acoplan con `.and.` las derivadas respecto a CN y respecto a carga `q`, dejando salidas `intent(out)` sin escribir

| Field | Value |
|---|---|
| library | dftd4 (`dftd4_model_d4` y `dftd4_model_d4s`, exportados públicamente con `!DEC$ ATTRIBUTES DLLEXPORT`) |
| location | `src/dftd4/model/d4.f90` (`weight_references` L270, `get_atomic_c6` L512–513, `get_polarizabilities` L671–672); `src/dftd4/model/d4s.f90` (`weight_references` L267, `get_atomic_c6` L520–521, `get_polarizabilities` L679–680) |
| versions | `4.3.1` (`main` commit `4a5bbf0ffbef983f6b225a4a275e61ea265088f6`) |
| detected | 2026-11-26 |
| detected_by | Antigravity (búsqueda guiada `implicit-api-contract`, ronda 5); pendiente de verificación independiente |
| oracle | Contrato de API (combinaciones parciales CN-solo y carga-solo con centinela `999.0_wp` frente a la llamada completa) |
| libraries | dftd4 |
| pattern | implicit-api-contract |
| oracle_class | api-contract |
| class | silent wrong result / uninitialized output, new |
| impact_level | III (en 6 métodos públicos de la API Fortran de `d4_model` y `d4s_model`, solicitar únicamente derivadas respecto a CN —p. ej. con cargas fijas— o únicamente respecto a `q` deja el arreglo `intent(out)` completamente sin escribir en silencio) |
| status | candidate (evidencia local: **defecto confirmado** — anomalía reproducida, causa aislada y parche verificado) |
| report | Pending |

#### Observation
Se probaron las llamadas parciales (solo CN o solo carga `q`) en los 6 métodos de `d4_model` y `d4s_model` sobre H$_2$O ([probe_dftd4_fortran.f90](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/scripts/probe_dftd4_fortran.f90)). Salidas crudas: [probe_dftd4_fortran.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/probe_dftd4_fortran.log) y [patched_probe_dftd4_fortran.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/patched_probe_dftd4_fortran.log).

| Método (`d4_model` / `d4s_model`) | Combinación | Baseline (`min` con centinela `999.0`) | Tras parche (`max|err|` vs ref) |
|---|---|---|---|
| `d4%weight_references` | `gwdcn` solo / `gwdq` solo | `9.99000E+02` / `9.99000E+02` | **`0.00000E+00`** / **`0.00000E+00`** |
| `d4%get_atomic_c6` | `dc6dcn` solo / `dc6dq` solo | `9.99000E+02` / `9.99000E+02` | **`0.00000E+00`** / **`0.00000E+00`** |
| `d4%get_polarizabilities` | `dadcn` solo / `dadq` solo | `9.99000E+02` / `9.99000E+02` | **`0.00000E+00`** / **`0.00000E+00`** |
| `d4s%weight_references` | `gwdcn` solo / `gwdq` solo | `9.99000E+02` / `9.99000E+02` | **`0.00000E+00`** / **`0.00000E+00`** |
| `d4s%get_atomic_c6` | `dc6dcn` solo / `dc6dq` solo | `9.99000E+02` / `9.99000E+02` | **`0.00000E+00`** / **`0.00000E+00`** |
| `d4s%get_polarizabilities` | `dadcn` solo / `dadq` solo | `9.99000E+02` / `9.99000E+02` | **`0.00000E+00`** / **`0.00000E+00`** |

#### Artifact filter
Compilación debug `-fcheck=all -Dopenmp=false`. Cada variante (`d4` y `d4s`) se compara contra su propia referencia de llamada completa.

#### Documentation and paper
Las seis rutinas forman parte de la interfaz pública exportada (`!DEC$ ATTRIBUTES DLLEXPORT`) de `d4_model` y `d4s_model`, y declaran las parejas `(gwdcn, dc6dcn/dadcn)` y `(gwdq, dc6dq/dadq)` como argumentos `optional` separados. Matemáticamente son derivadas parciales independientes ($\partial/\partial \mathrm{CN}$ y $\partial/\partial q$); en esquemas con cargas fijas o en acoplamientos autoconsistentes donde solo se necesita $\partial C_6 / \partial q$ o $\partial C_6 / \partial \mathrm{CN}$, solicitar una sola es un uso legítimo. Dejar `intent(out)` sin escribir viola el contrato. **Defecto inequívoco.**

#### Cause
En ambos archivos (`d4.f90` y `d4s.f90`):
- `weight_references`: `if (present(gwdcn) .and. present(gwdq))` en lugar de `.or.` con escrituras guardadas por `present(...)`.
- `get_atomic_c6`: `if (present(gwdcn) .and. present(dc6dcn) .and. present(gwdq) .and. present(dc6dq))` en lugar de `(present(gwdcn).and.present(dc6dcn)) .or. (present(gwdq).and.present(dc6dq))`.
- `get_polarizabilities`: `if (present(gwdcn) .and. present(dadcn) .and. present(gwdq) .and. present(dadq))` en lugar de `.or.` por pareja.

**Parche mínimo probado** ([dftd4-model.patch](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/patches/dftd4-model.patch)):
- Desacoplar las comprobaciones CN y `q` con `.or.` en los 6 métodos de `d4.f90` y `d4s.f90`.
- **Verificación del parche:**
  - Los 12 casos parciales dan `max|err| = 0.00000E+00` ([patched_probe_dftd4_fortran.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/patched_probe_dftd4_fortran.log)).
  - Identidad bit a bit en `d4s` completo: **incluida en los 67/67 registros hexadecimales IEEE-754 idénticos (0 diffs)** ([bit_identity_summary.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/bit_identity_summary.log)).
  - Suite oficial `meson test`: **54 OK, 1 Expected Fail, 0 Fail** ([patched_test_dftd4.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/patched_test_dftd4.log)).

#### Impact
- **Afectado:** Usuarios externos de la API Fortran pública `dftd4` (`d4_model` / `d4s_model`).
- **No afectado:** `get_dispersion` de alto nivel, C API, Python y CLI de `dftd4` (que pasan siempre ambos conjuntos de derivadas cuando calculan gradientes).

#### Verification
Pendiente de verificación independiente por un segundo agente/persona.

#### Corrections
Ninguna hasta la fecha.

---

## 3. Hallazgos secundarios / de capa baja y exclusiones explícitas

Se separan explícitamente de las 5 propuestas principales para mantener una evaluación honesta del impacto y evitar inflar el conteo:

1. **Rutinas de amortiguamiento y suma de Fourier en `simple-dftd3` y `dftd4` (capa interna-pública, acoplamiento 4-vía / 6-vía `.and.`):**
   - **`simple-dftd3`:** `get_dispersion2` en los 6 tipos de amortiguamiento (`rational.f90:220`, `zero.f90:243`, `mzero.f90:194`, `optimizedpower.f90:179`, `cso.f90:177`, `z.f90:176`), `get_atm_dispersion` (`atm.f90:433`) y las sumas de Fourier (`ewald.f90:47`, `spme.f90:78`) exigen `present(dc6dcn) .and. present(dEdcn) .and. present(gradient) .and. present(sigma)`. Probado dinámicamente en los 6 amortiguamientos + ATM ([probe_sdftd3_fortran.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/probe_sdftd3_fortran.log)): omitir `sigma` deja `gradient` sin modificar (`999.0`).
   - **`dftd4`:** `rational_damping_param%get_dispersion2` (`rational.f90:66`, L113–115), `get_dispersion3` (`rational.f90:362`) y `get_atm_dispersion` (`atm.f90:445`, L517–519) exigen los 6 opcionales a la vez (`dc6dcn .and. dc6dq .and. dEdcn .and. dEdq .and. gradient .and. sigma`). Probado dinámicamente ([probe_dftd4_fortran.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/probe_dftd4_fortran.log)): omitir `sigma` o `gradient` deja las salidas en `999.0`.
   - *Evaluación honesta:* Son bloques de construcción de bajo nivel (no expuestos en C/Python/CLI y llamados con argumentos completos por `get_dispersion`), aunque son módulos públicos con `optional`. Se documentan como hallazgos secundarios de capa baja (Nivel III/IV para llamadores Fortran de bajo nivel) y no se parchearon por separado dado que `get_dispersion` de alto nivel ya los envuelve con buffers locales.
2. **`tblite::get_selfenergy` (`src/tblite/xtb/h0.f90:114`, L141–143) — Nivel IV (menor):**
   - `dsedcn_en` (`intent(out), optional`) solo se escribe dentro de `if (present(cn_en))`. Probado dinámicamente ([probe_tblite_fortran.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/probe_tblite_fortran.log)): pasar `dsedcn_en` sin `cn_en` lo deja en `999.0`. A diferencia de los casos anteriores, aquí `cn_en` es una entrada lógicamente necesaria para la derivada dependiente de electronegatividad, aunque dejar un `intent(out)` sin inicializar sigue siendo un borde frágil (debería inicializarse a `0.0_wp` o documentarse explícitamente).
3. **`simple-dftd3::get_counterpoise_mpi` (`src/dftd3/mpi.f90:212`, L251) — Excluido del conteo nuevo:**
   - Contiene exactamente `grad = present(gradient) .and. present(sigma)`, el mismo defecto que `get_geometric_counterpoise` (**OSTRA-025**). Se anota únicamente como contraparte MPI de OSTRA-025 para completitud, **sin** doble contabilización.

---

## 4. Resultados negativos confirmados (rutinas limpias)

Documentar los negativos es parte esencial del experimento del catálogo:

| Librería | Rutina(s) | Verificación realizada | Por qué está limpia |
|---|---|---|---|
| `simple-dftd3` | `get_dispersion` (Fortran: `scalar`, `atomic`, `v2`, `mpi`) | Dinámica $2^3 = 8$ máscaras `(grad, sigma, hessian)` ([probe_sdftd3_fortran.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/probe_sdftd3_fortran.log)) | Usa `gradient_local` y `sigma_local` (`disp.f90:102, 153-154, 183-184`); error `0.000E+00` en las 8 máscaras. |
| `simple-dftd3` | `dftd3_get_dispersion` (C API) | Dinámica $2^2 = 4$ máscaras `NULL` con ASan ([probe_sdftd3_c.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/probe_sdftd3_c.log)) | ASan limpio y `0.000e+00` en todas las combinaciones. |
| `dftd4` | `dftd4_get_properties` (C API) y `get_properties` (Fortran) | Dinámica $2^4 = 16$ máscaras `NULL` con ASan ([probe_dftd4_c.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/probe_dftd4_c.log)) | Comprueba `present(cn)`, `present(q)`, `present(c6)`, `present(alpha)` de forma independiente; 16/16 limpias (`0.000e+00`). |
| `tblite` | `tblite_get_singlepoint` y `tblite_get_result_*` (C API) | Dinámica con ASan ([probe_tblite_c.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/probe_tblite_c.log)) + lectura | `api/calculator.f90:684` asigna siempre ambos `res%gradient` y `res%sigma`; cada getter consulta un único campo validando `NULL`. |
| `mctc-lib` & `tblite` | `get_wignerseitz_weights` | Dinámica $2^2$ con centinela `999.0` ([probe_mctc_ncoord.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/probe_mctc_ncoord.log), [probe_tblite_fortran.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/probe_tblite_fortran.log)) | Comprueba `present(dwdr)` y `present(dwdL)` independientemente (`0.000E+00`). |
| `tblite` | `ceh_singlepoint`, `cds_solvation%get_engrad`, `shift_solvation%get_engrad`, `container_list%get_engrad`, `alpb/cds%get_gradient`, `get_hamiltonian*` | Lectura estática exhaustiva | Sin salidas opcionales acopladas (`cds%get_engrad` usa `if (present(gradient))` solo; `container_list` despacha transparentemente; resto con salidas obligatorias). |
| Todas (4 librerías) | `get_*_version`, `get_pairwise_dispersion`, `dftd*_get_*_hessian`, constructores `new_*_ncoord`, `add_coordination_number_derivs` | Lectura estática | Comprobaciones independientes de `present()` o salidas obligatorias sin opcionales acoplados. |

---

## 5. Inventario de entregables y evidencia reproducible

- **Tabla completa de rutinas revisadas:** [rutinas_revisadas.csv](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/rutinas_revisadas.csv)
- **Parches mínimos (`git diff`, probados y compilados en `repos-patched/`):**
  - [mctc-lib-ncoord.patch](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/patches/mctc-lib-ncoord.patch)
  - [dftd4-disp-capi.patch](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/patches/dftd4-disp-capi.patch)
  - [dftd4-model.patch](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/patches/dftd4-model.patch)
  - [simple-dftd3-model.patch](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/patches/simple-dftd3-model.patch)
  - [tblite-singlepoint.patch](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/patches/tblite-singlepoint.patch)
- **Reproductores y arneses (`scripts/`):**
  - [probe_mctc_ncoord.f90](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/scripts/probe_mctc_ncoord.f90), [bitid_mctc.f90](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/scripts/bitid_mctc.f90)
  - [probe_dftd4_c.c](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/scripts/probe_dftd4_c.c), [probe_dftd4_fortran.f90](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/scripts/probe_dftd4_fortran.f90), [bitid_dftd4.f90](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/scripts/bitid_dftd4.f90)
  - [probe_sdftd3_c.c](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/scripts/probe_sdftd3_c.c), [probe_sdftd3_fortran.f90](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/scripts/probe_sdftd3_fortran.f90), [bitid_sdftd3.f90](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/scripts/bitid_sdftd3.f90)
  - [probe_tblite_c.c](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/scripts/probe_tblite_c.c), [probe_tblite_fortran.f90](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/scripts/probe_tblite_fortran.f90), [probe_tblite_virial_fd.f90](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/scripts/probe_tblite_virial_fd.f90), [bitid_tblite.f90](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/scripts/bitid_tblite.f90)
- **Salidas crudas (`outputs/`):**
  - Baseline: [probe_mctc_ncoord.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/probe_mctc_ncoord.log), [probe_dftd4_c.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/probe_dftd4_c.log), [probe_dftd4_fortran.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/probe_dftd4_fortran.log), [probe_sdftd3_c.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/probe_sdftd3_c.log), [probe_sdftd3_fortran.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/probe_sdftd3_fortran.log), [probe_tblite_c.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/probe_tblite_c.log), [probe_tblite_fortran.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/probe_tblite_fortran.log), [probe_tblite_virial_fd.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/probe_tblite_virial_fd.log)
  - Parcheado y suites oficiales: [patched_probe_mctc_ncoord.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/patched_probe_mctc_ncoord.log), [patched_probe_dftd4_c.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/patched_probe_dftd4_c.log), [patched_probe_dftd4_fortran.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/patched_probe_dftd4_fortran.log), [patched_probe_sdftd3_fortran.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/patched_probe_sdftd3_fortran.log), [patched_probe_tblite_c.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/patched_probe_tblite_c.log), [patched_probe_tblite_fortran.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/patched_probe_tblite_fortran.log), [bit_identity_summary.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/bit_identity_summary.log), [patched_test_mctc-lib.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/patched_test_mctc-lib.log), [patched_test_dftd4.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/patched_test_dftd4.log), [patched_test_simple-dftd3.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/patched_test_simple-dftd3.log), [patched_test_tblite.log](file:///Users/napolesd/Documents/Antigravity/PhysAud/exploracion5-contratos/outputs/patched_test_tblite.log)
