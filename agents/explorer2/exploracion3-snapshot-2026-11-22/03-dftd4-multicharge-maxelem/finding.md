# Cand-03: Dominio de elementos desacoplado entre `dftd4` (`max_elem = 118`) y `multicharge` (`max_elem = 103`): parámetros centinela `-1.0` aceptados en silencio

- **Librerías:** `dftd4` (commit `4a5bbf0ffbef983f6b225a4a275e61ea265088f6`, 2026-11-16) × `multicharge` (commit `4bb6c9637508bd1d21a28f72fa514fcf4e968fa0`, 2026-11-17)
- **Archivos fuente:**
  - `dftd4`: `src/dftd4/model/d4.f90` (líneas 106–112, 196, 205)
  - `multicharge`: `src/multicharge/param/eeq2019.f90` (línea 56, getters líneas 400+), `src/multicharge/param/eeqbc2025.f90` (línea 86), `src/multicharge/model/eeq.f90` (líneas 66–99)
- **Patrón del catálogo:** **Patrón 1 — Tablas de tamaño fijo desacopladas del índice máximo** (acoplamiento entre librerías: `max_elem = 118` en el llamador vs `max_elem = 103` en la tabla de parámetros EEQ)
- **Oráculo:** Positividad física de radios atómicos ($r_A > 0$) y dureza diagonal del modelo EEQ + propagación de `error_type` para elementos sin parametrizar
- **Clase:** Error numérico nuevo / contrato de rango de tabla silencioso
- **Nivel de evidencia:** **Defecto confirmado** (reproducido en `multicharge` compilado sin modificar; eliminado con parche mínimo en copia aislada)
- **Estado upstream:** Presente en `HEAD` de ambos repositorios (`2026-11-16`/`17`); no consta en `ostra-audit` ni en los issues consultados. Sin reportar (Paso 6).

---

## 1. Detección (Paso 1)

1. En `dftd4/src/dftd4/model/d4.f90:106-112`, la validación de elementos en `new_d4_model` rechaza explícitamente solo `Z = 104..111` (Rf–Rg) y **permite explícitamente `Z = 112..118` (Cn–Og)** porque `dftd4` dispone de tablas hasta `max_elem = 118`:
   ```fortran
   ! check for unsupported elements (104 (Rf) - 111 (Rg))
   do isp = 1, mol%nid
      if (mol%num(isp) > 103 .and. mol%num(isp) < 112) then
         call fatal_error(error, "Structure contains unsupported element ...")
         return
      end if
   end do
   ```
2. A continuación (`d4.f90:196` y `:205`), para los modelos de carga por defecto `d4_qmod%eeq` y `d4_qmod%eeqbc`, `dftd4` llama a `new_eeq2019_model(mol, d4%mchrg, error)` o `new_eeqbc2025_model(mol, d4%mchrg, error)` de `multicharge`.
3. Sin embargo, las tablas de `multicharge` (`param/eeq2019.f90:56` y `param/eeqbc2025.f90:86`) tienen **`max_elem = 103` (`H`–`Lr`)** y sus funciones `get_eeq_chi`, `get_eeq_rad`, `get_eeq_eta` y `get_eeq_kcnchi` devuelven el valor centinela **`-1.0_wp`** cuando `Z > 103`.
4. `new_eeq_model` (`multicharge/src/multicharge/model/eeq.f90:66-99`) copia `chi`, `rad`, `eta` y `kcnchi` directamente a `self` **sin verificar que `rad > 0` ni que `Z <= max_elem`**, y nunca asigna `error`.

### Evidencia observada (sin parchear)

Ejecutando los módulos originales sin modificar de `multicharge` sobre dímeros `X–O` (`R = 4.0 Bohr`, `evidence/unpatched.txt`):

| Dímero | `Z` | ¿`allocated(err)`? | `(chi, rad, eta, kcn)` de `X` | Cargas EEQ `Q(X), Q(O)` |
|---|---:|:---:|---|---|
| `Lr–O` | 103 | `F` | `(1.14, 1.34, -0.00, 0.01)` | `+0.55046, -0.55046` |
| **`Cn–O`** | **112** | **`F` (silencioso!)** | **`(-1.00, -1.00, -1.00, -1.00)`** | **`-1.55287, +1.55287`** |
| **`Og–O`** | **118** | **`F` (silencioso!)** | **`(-1.00, -1.00, -1.00, -1.00)`** | **`-1.55287, +1.55287`** |

Dado que las interacciones fuera de la diagonal usan `rad(izp)**2 = (-1.0)**2 = +1.0` (evitando `NaN` en `sqrt`), el sistema lineal se resuelve sin fallar, pero el término diagonal `self%eta(izp) + sqrt2pi / self%rad(izp)` (`eeq.f90:297, 328`) se vuelve **negativo** (`-1.0 - 0.7979 = -1.7979`), invirtiendo la física del modelo y asignando una carga altamente positiva (`+1.55 e`) al oxígeno frente al elemento superpesado (`-1.55 e`).

## 2. Filtrado de artefactos (Paso 2)

- Reproducido directamente en `multicharge` sin parchear. El fallo ocurre en la inicialización de parámetros (`new_eeq2019_model` -> `new_eeq_model`) y en el acoplamiento real `d4.f90:108, 196`. No depende de tolerancias ni del shim de `mctc_ncoord`.
- Nótese que el defecto afecta tanto a `multicharge` autónomo (para todo `Z > 103`, es decir, `104..118`) como a `dftd4` con `qmod = eeq` o `eeqbc` para `Z = 112..118` (los elementos que `d4.f90:108` deja pasar expresamente).

## 3. Cruce con documentación y código fuente (Paso 3)

- `param/eeq2019.f90` y `param/eeqbc2025.f90` definen deliberadamente `-1.0_wp` como centinela de «dato inexistente» cuando `Z > max_elem` (`103`), pero el constructor `new_eeq_model` nunca inspecciona ese centinela para emitir `fatal_error`.
- En `dftd4`, el modelo `gfn2` tiene cargas de referencia hasta `Z = 118`, pero `d4.f90:108` aplicó la misma lista de elementos permitidos (`1..103, 112..118`) a los tres modelos (`eeq`, `eeqbc`, `gfn2`) sin comprobar si `multicharge` soporta `112..118`.

## 4. Aislamiento de la causa y parche mínimo (Paso 4)

- **Líneas exactas:** `multicharge/src/multicharge/model/eeq.f90:91-94` (ausencia de validación de `rad <= 0.0_wp`) y `dftd4/src/dftd4/model/d4.f90:108, 196, 205`.
- **Parche mínimo verificado en copia aislada (`repro/minimal.patch`):** comprobar `if (any(rad <= 0.0_wp))` al inicio de `new_eeq_model` (y análogamente en `new_eeqbc_model`) y llamar a `fatal_error(error, ...)`. Como `dftd4/src/dftd4/model/d4.f90:197, 206` ya hace `if (allocated(error)) return`, este parche en `multicharge` protege automáticamente tanto a `multicharge` como a `dftd4`.
- **Resultado con parche (`evidence/patched.txt`):**
  - `Lr–O` (`Z = 103`): sigue funcionando idéntico (`Q = +0.55046, -0.55046`).
  - `Cn–O` (`Z = 112`) y `Og–O` (`Z = 118`): devuelven limpiamente `ERROR -> Missing EEQ parameters (non-positive atomic radius; unsupported element)`.

## 5. Impacto en observables (Paso 5)

- Cualquier cálculo EEQ (`multicharge`) con elementos trans-lawrencio (`Z > 103`) o cálculo DFT-D4 con elementos superpesados `Cn–Og` (`Z = 112..118`, explícitamente permitidos en `d4.f90`) produce silenciosamente cargas no físicas (p. ej. `Q(O) = +1.55 e` en `CnO`) y, en consecuencia, coeficientes $C_6(q)$ y energías de dispersión D4 corruptos sin aviso alguno.

## 6. Reproducibilidad

Ver `repro/repro_maxelem.f90`, `repro/minimal.patch`, `evidence/unpatched.txt` y `evidence/patched.txt`.
