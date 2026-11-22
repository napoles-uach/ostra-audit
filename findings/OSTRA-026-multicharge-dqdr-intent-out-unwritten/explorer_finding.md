# Cand-08: Contrato implícito en `multicharge` `get_charges` y `mchrg_model_type%solve`: omisión de `dqdL`/`sigma` deja `intent(out) :: dqdr` (y `gradient`) sin escribir

- **Librería:** `multicharge` (commit `4bb6c9637508bd1d21a28f72fa514fcf4e968fa0`, 2026-11-17)
- **Archivos fuente:** `src/multicharge/charge.f90` (líneas 53, 67, 125, 153) y `src/multicharge/model/type.F90` (líneas 219–221, 310–331)
- **Patrón del catálogo:** **Patrón 3 — Contratos de API implícitos (salidas no escritas / `intent(out)` sin inicializar)** (cuarto caso independiente, junto a `OSTRA-004` en `PySCF`, `Cand-06` en `libecpint` y `Cand-07` en `simple-dftd3`)
- **Oráculo:** Contrato Fortran de `intent(out)` + consistencia de $\partial q_i / \partial \mathbf{R}_j$ con diferencias finitas (FD)
- **Clase:** Contrato de API / memoria sin inicializar (`intent(out)`)
- **Nivel de evidencia:** **Defecto confirmado** (reproducido en `multicharge` sin modificar; corregido con parche mínimo en copia aislada, coincidiendo con FD a `1.46e-11`)
- **Estado upstream:** Presente en `HEAD` (`4bb6c963`); no consta en `ostra-audit` ni en issues conocidos de `grimme-lab/multicharge`. Sin reportar (Paso 6).

---

## 1. Detección (Paso 1)

En `multicharge/src/multicharge/charge.f90:38-67`, las rutinas públicas `get_charges`, `get_eeq_charges` y `get_eeqbc_charges` declaran `dqdr` y `dqdL` como argumentos `optional` con `intent(out)`:
```fortran
real(wp), intent(out), optional :: dqdr(:, :, :)
real(wp), intent(out), optional :: dqdL(:, :, :)
...
grad = (present(dqdr) .and. present(dqdL)) .or. present(cache)
```
Exactamente el mismo acoplamiento aparece en `mchrg_model_type%solve` (`src/multicharge/model/type.F90:219-221`):
```fortran
dcn = present(dcndr) .and. present(dcndL)
grad = present(gradient) .and. present(sigma) .and. dcn
cpq = ((present(dqdr) .and. present(dqdL)) .or. present(response)) .and. dcn
```
Si un llamador en un sistema molecular (no periódico, donde la derivada respecto a la red `dqdL` o el virial `sigma` no tiene uso) solicita únicamente `dqdr` (p. ej. `call get_charges(model, mol, err, q, dqdr=dqdr)` o `call get_eeq_charges(mol, err, q, dqdr=dqdr)`), `grad` evalúa a `.false.`, se omite todo el cálculo de derivadas y **`dqdr` — declarado `intent(out)`! — queda completamente sin escribir** (conteniendo memoria sin inicializar o valores residuales) sin emitir ningún `error`. Análogamente, llamar a `solve(..., gradient=g)` sin `sigma` omite silenciosamente `gradient`.

### Evidencia observada (sin parchear)

En el dímero `O–H` ($R = 1.8\ \text{Bohr}$, modelo `EEQ2019`, no periódico, `evidence/unpatched.txt`):

| Llamada | `Q(O), Q(H)` | $\partial q_{\text{O}} / \partial z_{\text{H}}$ | Error vs FD |
|---|---:|---:|---:|
| Referencia FD ($h = 10^{-5}$) | `-0.36150, +0.36150` | `+2.0865449037e-01` | — |
| `get_charges(..., q, dqdr=dqdr, dqdL=dqdL)` | `-0.36150, +0.36150` | `+2.0865449039e-01` | `2.0e-10` |
| **`get_charges(..., q, dqdr=dqdr)` (sin `dqdL`, valor previo `-999`)** | `-0.36150, +0.36150` | **`-9.9900000000e+02`** | **Sin escribir (`intent(out)`)!** |

## 2. Filtrado de artefactos (Paso 2)

- Reproducido directamente en los 12 archivos fuente originales de `multicharge` (`charge.f90`, `model/type.F90`, `model/eeq.f90`, etc.) sin modificar.
- La referencia FD y la llamada con `dqdL` concuerdan a 10 dígitos (`2.0865449037e-01` vs `2.0865449039e-01`), descartando cualquier artefacto del shim de `mctc_ncoord`.

## 3. Cruce con documentación y código fuente (Paso 3)

- Ni la firma ni la documentación de `get_charges` / `get_eeq_charges` / `get_eeqbc_charges` indican que `dqdr` requiera pasar también `dqdL` para ser escrito.
- De hecho, en `type.F90:322` los autores ya escribieron `if (present(dqdL)) call gemm(...)` (y en `326`: `if (present(dqdr)) ...`), pero en la línea `330` faltó `if (present(dqdr))` y en `charge.f90:67` y `type.F90:220-221` usaron `.and.` en lugar de `.or.`.
- Es notable que exactamente el mismo idioma frágil (`present(grad) .and. present(virial)`) aparece tanto en `simple-dftd3` (`Cand-07`) como en `multicharge` (`Cand-08`), mostrando un patrón recurrente en librerías Fortran del ecosistema Grimme-lab.

## 4. Aislamiento de la causa y parche mínimo (Paso 4)

- **Líneas exactas:** `src/multicharge/charge.f90:67` y `src/multicharge/model/type.F90:220-221, 310-313, 330`.
- **Parche mínimo verificado en copia aislada (`repro/minimal.patch`):** cambiar `.and.` por `.or.` para `grad` y `cpq`, y proteger las escrituras opcionales con `if (present(gradient))`, `if (present(sigma))`, `if (present(dqdr))` y `if (present(dqdL))`.
- **Resultado con parche (`evidence/patched.txt`):** `get_charges(..., q, dqdr=dqdr)` sin `dqdL` devuelve `+2.0865449039e-01`, coincidiendo con la referencia FD a **`1.4550e-11`**.

## 5. Impacto en observables (Paso 5)

- **Crítico** en cualquier código molecular que acople cargas EEQ a gradientes (p. ej. gradientes de dispersión dependiente de carga o campos electrostáticos) pidiendo solo `dqdr` sin `dqdL`: en Fortran un argumento `intent(out)` no asignado queda indefinido (comportamiento indefinido / basura arbitraria en memoria), corrompiendo por completo las derivadas de carga y las fuerzas dependientes de $\partial q / \partial \mathbf{R}$ sin error ni aviso.

## 6. Reproducibilidad

Ver `repro/repro_dqdr_contract.f90`, `repro/minimal.patch`, `evidence/unpatched.txt` y `evidence/patched.txt`.
