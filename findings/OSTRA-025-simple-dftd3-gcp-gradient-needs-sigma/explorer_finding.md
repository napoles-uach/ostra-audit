# Cand-07: Contrato implícito en `simple-dftd3` `get_geometric_counterpoise` y API C `dftd3_get_counterpoise`: omisión de `sigma` deja `gradient` sin escribir en silencio

- **Librería:** `simple-dftd3` (`v1.7.0`, commit `84cb7cddd8d6669c4ecd59b06d53e0d373656239`, 2026-10-28)
- **Archivos fuente:** `src/dftd3/gcp.f90` (línea 104), `src/dftd3/api.f90` (líneas 1367–1428); contraste interno con `src/dftd3/disp.f90` (línea 102)
- **Patrón del catálogo:** **Patrón 3 — Contratos de API implícitos (salidas opcionales que no se escriben sin aviso)** (tercer caso independiente junto a `OSTRA-004` en `PySCF` y `Cand-06` en `libecpint`)
- **Oráculo:** Consistencia del gradiente analítico con diferencias finitas (FD) + postcondición de salida de argumento opcional `gradient`
- **Clase:** Contrato de API / error silencioso (error del 100 % en el gradiente)
- **Nivel de evidencia:** **Defecto confirmado** (reproducido en `simple-dftd3` original; corregido con parche mínimo en copia aislada, coincidiendo con FD a `1.24e-11`)
- **Estado upstream:** Presente en `v1.7.0` / `HEAD` (`84cb7cdd`); no consta en `ostra-audit` ni en issues conocidos de `dftd3/simple-dftd3`. Sin reportar (Paso 6).

---

## 1. Detección (Paso 1)

En `simple-dftd3`, la subrutina pública `get_geometric_counterpoise` (y `get_geometric_counterpoise_atomic`, `src/dftd3/gcp.f90:104`) declara `gradient` y `sigma` como argumentos `optional` independientes, pero decide si calcular derivadas mediante un **`.and.`** en lugar de **`.or.`**:
```fortran
! src/dftd3/gcp.f90:104
grad = present(gradient) .and. present(sigma)
```
En contraste directo, dentro de la **misma librería**, la rutina de dispersión D3 (`src/dftd3/disp.f90:102`) usa **`.or.`** (`grad = present(gradient) .or. present(sigma)`), permitiendo solicitar `gradient` sin `sigma` (el caso estándar para toda molécula no periódica!).

Lo mismo se propaga a la **API C** (`src/dftd3/api.f90:1367-1428`, `dftd3_get_counterpoise`): si el usuario en C/Python pasa un arreglo `gradient` válido y `sigma = NULL` (documentado como puntero opcional para sistemas moleculares), `c_associated(c_sigma)` es falso, se omite `sigma` al llamar a Fortran, `grad` evalúa a `.false.`, se llama a `gcp_energy` en lugar de `gcp_deriv`, y **`gradient` se devuelve lleno de ceros (en C) o completamente sin escribir (en Fortran) sin ningún error ni advertencia**.

### Evidencia observada (sin parchear)

En el dímero molecular `C–N` ($R = 2.5\ \text{Bohr}$, `r2SCAN-3c`, no periódico, `evidence/unpatched.txt`):

| Llamada | `Energy` (Eh) | `Grad_z` (Eh/Bohr) | Error absoluto vs FD |
|---|---:|---:|---:|
| Referencia FD ($h = 10^{-5}$) | — | `-2.3043983216e-03` | — |
| `get_geometric_counterpoise(..., gradient=g, sigma=s)` | `+4.7649130194e-03` | `-2.3043983340e-03` | `1.24e-11` |
| **`get_geometric_counterpoise(..., gradient=g)` (`g` inicializado a `-999`)** | `+4.7649130194e-03` | **`-9.9900000000e+02`** | **Sin escribir!** |
| **`get_geometric_counterpoise(..., gradient=g)` (`g = 0`, equivalente API C con `sigma=NULL`)** | `+4.7649130194e-03` | **`0.0000000000e+00`** | **`2.3044e-03` (100 %!)** |

## 2. Filtrado de artefactos (Paso 2)

- Reproducido directamente compilando el código original `src/dftd3/gcp.f90`.
- La energía es idéntica (`4.7649130194e-3 Eh`) con o sin `sigma`, lo que hace el fallo aún más silencioso: el llamador recibe una energía plausible y exacta junto a un gradiente exactamente nulo (o basura no escrita).

## 3. Cruce con documentación y código fuente (Paso 3)

- Ni la firma ni los comentarios de `get_geometric_counterpoise` (`gcp.f90:52-97`) ni de `dftd3_get_counterpoise` (`api.f90:1367-1428`) documentan que `gradient` sea ignorado si `sigma` se omite.
- La inconsistencia directa con `disp.f90:102` (`grad = present(gradient) .or. present(sigma)`) confirma que se trata de un contrato implícito no intencionado (probablemente introducido porque `gcp_deriv` requiere ambos argumentos no opcionales internamente).

## 4. Aislamiento de la causa y parche mínimo (Paso 4)

- **Línea exacta:** `src/dftd3/gcp.f90:104` (`grad = present(gradient) .and. present(sigma)`).
- **Parche mínimo verificado en copia aislada (`repro/minimal.patch`):** cambiar `.and.` por `.or.`, acumular en buffers locales `grad_loc(3, mol%nat)` y `sigma_loc(3, 3)`, y sumar selectivamente con `if (present(gradient))` y `if (present(sigma))`.
- **Resultado con parche (`evidence/patched.txt`):** `get_geometric_counterpoise(..., gradient=g)` sin `sigma` devuelve `Grad_z = -2.3043983340e-03`, coincidiendo con la referencia FD a **`1.2426e-11`** (error relativo `5.4e-9`).

## 5. Impacto en observables (Paso 5)

- **Crítico (100 % de error en la fuerza gCP)** en cualquier interfaz molecular (optimización de geometría, dinámica molecular, cálculo de frecuencias por FD del gradiente) que invoque `get_geometric_counterpoise` o `dftd3_get_counterpoise(..., gradient, NULL)` omitiendo el virial `sigma` por tratarse de un sistema no periódico: la contribución gCP a las fuerzas desaparece silenciosamente mientras la energía sí incluye gCP, rompiendo la consistencia energía–fuerza.

## 6. Reproducibilidad

Ver `repro/repro_gcp_grad_contract.f90`, `repro/minimal.patch`, `evidence/unpatched.txt` y `evidence/patched.txt`.
