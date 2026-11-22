# Exploración 3 (Explorador 2) — Segundos casos para el catálogo OSTRA

Esta carpeta contiene los resultados verificados de la fase 3 de auditoría física enfocada en encontrar **segundos (y terceros/cuartos) casos independientes**, idealmente en librerías distintas, para los tres patrones del catálogo de `~/Documents/ostra-audit/catalog/README.md` que hasta ahora contaban con un único caso:

1. **Patrón 1:** Tablas de tamaño fijo desacopladas del índice máximo que el código puede solicitar (caso existente: `OSTRA-008` en `libcint`).
2. **Patrón 2:** Comentarios o documentación que afirman un hecho numérico no verificado (caso existente: `OSTRA-010` en `libxc`).
3. **Patrón 3:** Contratos de API implícitos — salidas opcionales o arreglos que no siempre se escriben sin aviso (caso existente: `OSTRA-004` en `PySCF libcgto`).

Ninguno de los 8 candidatos siguientes estaba registrado en `~/Documents/ostra-audit` (`findings/`, `candidates/`, `discarded/`), todos fueron verificados contra las ramas `HEAD`/`master` actuales (noviembre 2026) y **todos cuentan con reproductor mínimo y parche mínimo aplicado sobre copia aislada** que elimina la anomalía observada. **Ningún hallazgo ha sido reportado upstream (Paso 6).**

---

## Tabla resumen de candidatos

| # | Candidato (carpeta) | Librería (versión / commit) | Patrón del catálogo | Nivel de evidencia | Estado |
|---|---|---|---|---|---|
| 01 | [01-simple-dftd3-gcp-zeff-oob](01-simple-dftd3-gcp-zeff-oob/finding.md) | `simple-dftd3` (`v1.7.0`, `84cb7cdd`) | **Patrón 1** (tablas fijas 35/36 vs `zeff = 36..39, 95..118`) | **Defecto confirmado** | Nuevo; sin reportar |
| 02 | [02-pyscf-ecp-fixed-tables](02-pyscf-ecp-fixed-tables/finding.md) | `PySCF libcgto` (`2.14.0` y `master`) | **Patrón 1** (`z==16.0` `entry=400` OOB + `omega_nuc[CART_CUM]` stack smash; + nota Patrón 2 `order=32`) | **Defecto confirmado** | Nuevo en `master` (`_factorial[24]` de 2.14.0 ya corregido en `master`); sin reportar |
| 03 | [03-dftd4-multicharge-maxelem](03-dftd4-multicharge-maxelem/finding.md) | `dftd4` (`4a5bbf0f`) × `multicharge` (`4bb6c963`) | **Patrón 1** (`max_elem = 118` en `dftd4` vs `103` en `multicharge`, centinela `-1.0` silencioso) | **Defecto confirmado** | Nuevo; sin reportar |
| 04 | [04-libecpint-bessel-upper-bound](04-libecpint-bessel-upper-bound/finding.md) | `libecpint` (`v1.1`, `d7966f19`) | **Patrón 2** (`upper_bound(z, L)` es cota *inferior* para $L \ge 1$ en $(0, z_{\max}(L))$, subestima hasta `1,958×`) | **Defecto confirmado** | Nuevo; sin reportar |
| 05 | [05-simple-dftd3-gcp-bint-baux](05-simple-dftd3-gcp-bint-baux/finding.md) | `simple-dftd3` (`v1.7.0`, `84cb7cdd`) | **Patrón 2** (`"12 terms seem accurate enough"` + cruce en `|za-zb| < 0.1` $\to$ cancelación en `baux_all`/`g3s3s`) | **Defecto confirmado** | Nuevo (clase: precisión numérica); sin reportar |
| 06 | [06-libecpint-bessel-contract](06-libecpint-bessel-contract/finding.md) | `libecpint` (`v1.1`, `d7966f19`) | **Patrón 3** (`calculate(z<=0)` deja `values[1..maxL]` sin escribir y la versión escalar devuelve `1.0` para $L \ge 1$) | **Defecto confirmado** | Nuevo; sin reportar |
| 07 | [07-simple-dftd3-gcp-grad-contract](07-simple-dftd3-gcp-grad-contract/finding.md) | `simple-dftd3` (`v1.7.0`, `84cb7cdd`) | **Patrón 3** (`get_geometric_counterpoise` y API C dejan `gradient` sin escribir / en ceros si se omite `sigma`) | **Defecto confirmado** | Nuevo; sin reportar |
| 08 | [08-multicharge-dqdr-contract](08-multicharge-dqdr-contract/finding.md) | `multicharge` (`4bb6c963`) | **Patrón 3** (`get_charges` deja `intent(out) :: dqdr` sin escribir si se omite `dqdL`; ídem `solve` con `gradient` sin `sigma`) | **Defecto confirmado** | Nuevo; sin reportar |

---

## Cobertura lograda de los tres patrones prioritarios

- **Patrón 1 (Tablas fijas desacopladas):** `OSTRA-008` (`libcint`) + **3 casos nuevos independientes** en `simple-dftd3`, `PySCF libcgto master` y el acoplamiento `dftd4 × multicharge`.
- **Patrón 2 (Comentarios numéricos no verificados):** `OSTRA-010` (`libxc`) + **2 casos nuevos independientes** en `libecpint` y `simple-dftd3` (más evidencia adicional en el nuevo comentario `order=32` de `PySCF master`).
- **Patrón 3 (Contratos de API implícitos / salidas no escritas):** `OSTRA-004` (`PySCF libcgto`) + **3 casos nuevos independientes** en `libecpint`, `simple-dftd3` y `multicharge`.

## Otros documentos en esta carpeta

- [searched.md](searched.md): registro completo de todas las búsquedas (patrón × librería × método), incluidos los resultados negativos, falsos positivos descartados por el Paso 2 y observaciones secundarias en `spglib`.
- [lda_c_pk09_response.md](lda_c_pk09_response.md): respuesta exacta al pendiente de la Tanda 1 (`lda_c_pk09` requiere `XC_POLARIZED` con $0 < |\zeta| < 1$, p. ej. $\zeta = 0.4$, alrededor de $r_{s,\text{polo}} = 22.739412132480812$).
- `repos/`: clones superficiales de solo lectura (`--depth 1`, sin modificar) de `libecpint`, `simple-dftd3`, `dftd4`, `multicharge` y `spglib`.
