# OSTRA-020 (parte B): Aislamiento completo de las caídas de integrales ECP con derivadas en PySCF (`master 484fd7c`)

- **Librería:** `PySCF` (`libcgto`: `pyscf/lib/gto/nr_ecp.h`, `nr_ecp.c`, `nr_ecp_deriv.c`)
- **Versión/Commit auditado:** `master` `484fd7c114ce51ca007140190f4adc9f6a50f1c9` (2026-11-22); verificado también sobre el release instalado `2.14.0` (`nr_ecp.h` y `nr_ecp_deriv.c` son idénticos byte a byte a `master`; `nr_ecp.c` solo difiere en `_factorial[65]`).
- **Oráculo:** AddressSanitizer (`clang -O1 -g -fsanitize=address`) en compilación aislada de `master 484fd7c` + consistencia de gradiente nuclear analítico (`mf.nuc_grad_method()`) contra diferencias finitas de 4 puntos (FD4) + regresión bit a bit para $l_b \le 5$.
- **Clase:** Seguridad de memoria (desbordamientos de pila y de *heap*), nuevo.
- **Estado upstream:** Verificado en `master 484fd7c`; no encontrado en los issues de `pyscf/pyscf` (los issues históricos `#866` y `#1065` tratan fallos distintos). Sin reportar (Paso 6).

---

## 1. Causa exacta de la caída con $l_b = 6$ y cualquier $l_c$ (y hallazgos adicionales revelados por ASan)

### Por qué $l_b = 6$ falla para **cualquier** $l_c$ (incluso $l_c = 0$ o ECP puramente local)

La verificación observó que `mol.intor("ECPscalar_ipnuc")` fallaba con $l_b = 6$ y $l_c \in \{0, 2, 3, 4\}$ aunque $l_b + 1 + l_c \le 11 < 12$. La razón confirmada por AddressSanitizer (`evidence/unpatched_asan_traces.txt`, Traza 1) es que **`omega_nuc[CART_CUM]` (`CART_CUM = 456`, válido solo hasta $l_{\max} = 12$) existe en DOS funciones independientes**:

1. **Canal local (`ECPtype1_cart` $\to$ `type1_rad_ang`, `nr_ecp.c:5253`):**
   El canal local ($l = -1$) **no depende de $l_c$**, sino de la suma de los momentos angulares de los dos orbitales de base:
   $$l_{\max}^{\text{local}} = l_i^{\text{eff}} + l_j^{\text{eff}}.$$
   En `ecp_cartcum.py`, `mol.intor("ECPscalar_ipnuc")` evalúa todos los pares de capas `(ish, jsh)`, incluido el bloque diagonal `(0, 0)` con $l_i = l_b = 6$ y $l_j = l_b = 6$. `_deriv1_cart` (`nr_ecp_deriv.c:254-256`) eleva $l_i \to l_i + 1 = 7$ y llama a `ECPtype1_cart` con $l_{\max}^{\text{local}} = 7 + 6 = 13 > 12$.
   - **Primer acceso inválido (ASan, Traza 1):** escritura de 8 bytes en `ang_nuc_in_cart` (`nr_ecp.c:5023`) llamada desde `type1_rad_ang` (`nr_ecp.c:5258`), desbordando `double omega_nuc[CART_CUM]` (`nr_ecp.c:5253`, tamaño `456 * 8 = 3648 B`, offset `3680`) por **`105` doubles (`840` bytes)** — **independientemente de $l_c$** (ocurre incluso con un ECP que solo tiene canal local $l = -1$, caso `lb6_loconly`!).
   - En cambio, el bloque fuera de la diagonal `(0, 1)` ($l_i = 6, l_j = 0 \Rightarrow l_{\max}^{\text{local}} = 7 \le 12$) **no** falla (`OK shls=(0,1)`).
   - Por la misma razón, `ECPscalar` (sin derivada) con $l_b = 6$ alcanza exactamente $l_{\max}^{\text{local}} = 6 + 6 = 12 \le 12$ y funciona, pero falla a partir de $l_b = 7$ ($7 + 7 = 14 > 12$).

2. **Canal no local (`ECPtype2_cart` $\to$ `type2_facs_ang`, `nr_ecp.c:5078`):**
   Depende de $l_{\max}^{\text{no-local}} = \max(l_i^{\text{eff}}, l_j^{\text{eff}}) + l_c$. En `lb = 5, lc = 7`, el canal local cabe ($6 + 5 = 11 \le 12$), pero `type2_facs_ang` (`nr_ecp.c:5083`) desborda `omega_nuc[456]` (`nr_ecp.c:5078`) con $l_i^{\text{eff}} + l_c = 6 + 7 = 13 > 12$ (ASan, Traza 2).

### Tres defectos adicionales de memoria descubiertos al barrer derivadas segundas y $l_c > 5$ con ASan

Al extender el barrido ASan a todas las integrales con derivada en `nr_ecp_deriv.c`, afloraron **tres defectos adicionales** en los mismos archivos de `master 484fd7c`:

3. **Subestimación del búfer de *heap* en `ECPscalar_cache_size` para derivadas (`nr_ecp.c:6048`, `nr_ecp_deriv.c:302, 391`):**
   `_cart_factory` y `_sph_factory` calculan `cache_size` con los momentos angulares **sin incrementar** `(li, lj)`, aunque luego `_deriv1_cart` / `_igv_cart` elevan `li -> li + 1`, `_ipipv_cart` eleva **`li -> li + 2`**, y `_ipvip_cart` eleva **`li -> li + 1, lj -> lj + 1`**.
   - **ASan, Traza 3 (`lb = 5, lc = 5`, `ECPscalar_ipipnuc` y `ECPscalar_ipiprinv`):** **`SEGV` en `cache_3dfac` (`nr_ecp.c:5054`)** $\leftarrow$ `type2_facs_ang` (`nr_ecp.c:5145`) $\leftarrow$ `ECPtype2_cart` (`nr_ecp.c:5535`) — **ya con funciones `h` ($l_b = 5$, p. ej. `cc-pV5Z`) y $l_c = 5$**, donde $l_{\max} = 12 \le 12$ cabía en `CART_CUM`, pero `cache` se quedó corto al crecer `li` de `5` a `7`!
4. **`ECP_LMAX = 5` fijado en duro al reservar `cache` en `ECPtype2_cart` (`nr_ecp.c:5407-5419`) y `ECPscalar_cache_size` (`nr_ecp.c:6062-6089`):**
   Cuando el ECP contiene proyectores con $l_c > 5$, `ECPtype2_cart` reserva `rad_all`, `angi`, `angj`, `buf`, `radi`, `radj` usando `ECP_LMAX = 5`, pero itera con el $l_c$ real.
   - **ASan, Traza 4 (`lb = 5, lc = 6`, `ECPscalar_ipnucip` y `ECPscalar_iprinvip`):** **`heap-buffer-overflow` en `ECPsph_ine` (`nr_ecp.c:4675`)** $\leftarrow$ `type2_facs_rad` (`nr_ecp.c:5199`).
5. **Tablas auxiliares acopladas e inicialización sobredimensionada en `nr_ecp.c`:**
   - `_offset_cart[15]` (`nr_ecp.c:4567`) solo llega hasta índice `14` ($l \le 13$; acceder `_offset_cart[lmax+1]` en `type2_facs_ang:5085` desborda con $l_{\max} \ge 14$, p. ej. `lb = 6` + derivadas segundas).
   - `_cart_pow_y[120]` y `_cart_pow_z[120]` (`nr_ecp.c:4541-4556`) solo cubren hasta $l = 14$ (`120` entradas), mientras `cint.h` define `ANG_MAX = 15`, `CART_MAX = 136` y `ang_nuc_in_cart` usa `xx[16]` (hasta $l = 15$).
   - En `ECPtype1_cart` (`nr_ecp.c:5930`), `rad_all` se reserva para `npi*npj*d2`, pero la línea de puesta a cero escribe `d3 = d1*d2` elementos.
   - En `ECPtype2_cart` (`nr_ecp.c:5408-5413, 5535`), `RESTORE_STACK` puede dejar `cache` por debajo del final de `angi/angj/buf` cuando el bloque angular supera al radial (p. ej. `li' = 8, lj = 6, lc = 5`), haciendo que `type2_facs_ang` sobrescriba `angj/buf`.

---

## 2. Límite real por combinación de $(l_b, l_c)$ y orden de derivada

Con una capa de base de momento $l_b$ y un ECP con canal local ($l = -1$) y proyectores hasta $l_c$, las condiciones de seguridad de memoria en `master 484fd7c` **sin parchear** son:

| Integral ECP | Desplazamiento $(\Delta l_i, \Delta l_j)$ | Falla `omega_nuc[CART_CUM]` (local, si hay $l=-1$) | Falla `omega_nuc[CART_CUM]` (no local) | Otros fallos ASan en `master 484fd7c` |
|---|---|---|---|---|
| `ECPscalar` | $(0, 0)$ | $2 l_b > 12 \Rightarrow \mathbf{l_b \ge 7}$ | $l_b + l_c > 12$ | $l_c > 5$ solapa sub-búferes en `ECPtype2_cart` |
| `ECPscalar_ipnuc`, `ECPscalar_iprinv`, `ECPscalar_ignuc` | $(+1, 0)$ | $2 l_b + 1 > 12 \Rightarrow \mathbf{l_b \ge 6}$ (cualquier $l_c$) | $l_b + 1 + l_c > 12$ (p. ej. $l_b=5, l_c \ge 7$) | Subestimación de `cache` para $\Delta l_i = +1$ |
| `ECPscalar_ipipnuc`, `ECPscalar_ipiprinv` | $(+2, 0)$ | $2 l_b + 2 > 12 \Rightarrow \mathbf{l_b \ge 6}$ (cualquier $l_c$) | $l_b + 2 + l_c > 12$ (p. ej. $l_b=5, l_c \ge 6$) | **`SEGV` en `cache_3dfac` ya en $l_b = 5, l_c = 5$** ($\Delta l_i = +2$ ignorado en `cache_size`) |
| `ECPscalar_ipnucip`, `ECPscalar_iprinvip` | $(+1, +1)$ | $2 l_b + 2 > 12 \Rightarrow \mathbf{l_b \ge 6}$ (cualquier $l_c$) | $l_b + 1 + l_c > 12$ (p. ej. $l_b=5, l_c \ge 7$) | **`heap-buffer-overflow` en `ECPsph_ine` en $l_b = 5, l_c = 6$** (`ECP_LMAX = 5` + $\Delta l_i = \Delta l_j = +1$) |

El resultado completo de los **96 casos** probados bajo ASan (`8` integrales $\times$ `12–13` combinaciones, incluyendo capas contraídas `nprim = 3` y ECP solo local `lb6_loconly`) está en `evidence/asan_matrix_comparison.txt`: **58 fallos ASan sin parchear $\to$ 0 fallos ASan (96/96 `PASS`) con el parche mínimo**.

El techo duro posterior al parche queda determinado por `libcint` (`ANG_MAX = 15` en `CINTc2s_bra_sph`/`CINTs2c_bra_sph`): soporta de forma limpia hasta $l_{\max} \le 15$, es decir:
- `ECPscalar`: hasta $l_b \le 7$ ($2 l_b = 14 \le 15$).
- Derivadas primeras (`ipnuc`, `iprinv`, `ignuc`): hasta $l_b \le 7$ ($2 l_b + 1 = 15 \le 15$).
- Derivadas segundas (`ipipnuc`, `ipiprinv`, `ipnucip`, `iprinvip`): hasta $l_b \le 6$ ($2 l_b + 2 = 14 \le 15$).

---

## 3. Corrección mínima y verificación (ASan, regresión y gradiente nuclear FD4)

El parche unificado contra `master 484fd7c` está en `repro/minimal.patch` (generado por `repro/make_patch.py`):
1. `nr_ecp.h`: `CART_CUM` `456` $\to$ `816` ($\sum_{k=0}^{15} \binom{k+2}{2} = 816$, alineado con `ANG_MAX = 15` y `CART_MAX = 136` de `cint.h`).
2. `nr_ecp.c`:
   - Extender `_offset_cart[]` a 17 entradas (`0..16`, terminando en `680, 816`) y `_cart_pow_y[]`/`_cart_pow_z[]` a 136 entradas ($l = 0..15$).
   - Calcular `ecp_lmax = MAX(ECP_LMAX, max_lc)` dinámicamente desde `ecpbas` en `ECPtype2_cart` y en el nuevo `ECPscalar_cache_size_deriv(comp, dli, dlj, ...)`, manteniendo sin cambios la ABI pública de `ECPscalar_cache_size` (que delega con `dli = 0, dlj = 0`).
   - Proteger el puntero `ang_work = max(cache, cache_ang)` antes de llamar a `type2_facs_ang`, y corregir la inicialización de `rad_all` en `ECPtype1_cart` (`d3` $\to$ `npi*npj*d2`).
3. `nr_ecp_deriv.c`: pasar `(dli, dlj)` a las funciones `static` `_cart_factory` y `_sph_factory` (`(1, 0)` para `ipnuc`/`iprinv`/`ignuc`, `(2, 0)` para `ipipnuc`/`ipiprinv`, `(1, 1)` para `ipnucip`/`iprinvip`).

### Resultados de verificación

1. **ASan limpio (`evidence/asan_matrix_comparison.txt`):** **96/96 `PASS`** (0 errores de AddressSanitizer frente a 58 fallos sin parchear).
2. **Regresión bit a bit para $l_b \le 5$ (`evidence/regression_lb_le_5.txt`):** en `scalar`, `ipnuc`, `iprinv`, `ignuc`, `ipipnuc` e `ipnucip` para `(lb=4, lc=4)` y `(lb=5, lc=4)`, **`max |unpatched - patched| = 0.000e+00`** en todos los casos (idéntico bit a bit).
3. **Gradiente nuclear analítico (`mf.nuc_grad_method().kernel()`) vs FD4 con funciones `i` ($l_b = 6$) (`evidence/nuc_grad_fd_validation.txt`):**
   - **`NaH` con capa `i` ($l_b = 6$) + ECP local y no local (`[-1, 0, 1]`):**
     - Sin parchear: **aborta con `exit code -6` (`SIGABRT`)** dentro de `mf.nuc_grad_method().kernel()`.
     - Parcheado (`libcgto_patched.dylib` compilado desde `master 484fd7c + minimal.patch`): `grad_z = +5.994325447891e-02 Eh/Bohr`, suma traslacional `1.11e-16`, diferencia frente a FD4:
       - `h = 1e-2 Bohr`: `|anal - FD4| = 6.414e-11` (`rel = 1.07e-09`)
       - `h = 1e-3 Bohr`: **`|anal - FD4| = 5.168e-14` (`rel = 8.62e-13`)**
       - `h = 1e-4 Bohr`: `|anal - FD4| = 1.567e-12` (`rel = 2.62e-11`)
   - **`AgH` con bases/ECP integrados reales (`Ag: cc-pv5z-pp` + ECP `def2-svp`, `H: sto-3g`, `nao = 149`, $R = 3.06\ \text{Bohr}$):**
     - Sin parchear: RHF converge (`E = -146.647710898546 Eh`), pero **`mf.nuc_grad_method().kernel()` aborta con `exit code -6` (`SIGABRT`)**.
     - Parcheado: `grad_z = -1.626501750585e-02 Eh/Bohr`, suma traslacional `8.88e-16`, frente a FD4 (`h = 1e-3 Bohr`): **`|anal - FD4| = 1.219e-10` (`rel = 7.50e-09`)**.

---

## 4. Relevancia práctica en química computacional

La caída **afecta directamente a gradientes nucleares reales en PySCF usando únicamente bases y ECPs estándar incluidos de serie en PySCF** (sin bases exóticas ni `basis-set-exchange`):

- **Bases integradas en PySCF con funciones `i` ($l = 6$):** un escaneo de los 186 archivos `.dat` de `pyscf/gto/basis/` encontró **27 conjuntos de bases integrados** con funciones `i`, incluidos:
  - **`cc-pv5z-pp`** (21 metales pesados: `Ag, Au, Cd, Cu, Hf, Hg, Ir, Mo, Nb, Os, Pd, Pt, Re, Rh, Ru, Ta, Tc, W, Y, Zr, Zn`) y **`cc-pwCV5Z-PP`** (`At, Bi, Pb, Po, Rn, Tl`).
  - **`cc-pv5z`**, **`aug-cc-pv5z`**, **`pc-4`**, **`pcseg-4`** en todos los metales de transición 3d (`Sc`–`Zn`), y **`cc-pvqz-dk`** en 30 lantánidos/actínidos.
  - Además, fuera de los `.dat` integrados, `cc-pV6Z` / `aug-cc-pV6Z` (vía `basis-set-exchange`) introducen funciones `i` en elementos ligeros (`C, N, O, F, ...`).
- **Condición exacta para que se dispare en cálculos reales:** basta que el sistema combine **cualquier base con funciones `i` ($l = 6$)** con **cualquier ECP que posea canal local no nulo (`ANG_OF == -1`)** — como **`def2-ECP`** (`def2-svp`/`tzvp`/`qzvp`), **`lanl2dz`**, **`sbkjc`**, **`crenbl`**, **`bfd`** o **`ccECP`** — y que exista solapamiento no despreciable entre la capa `i` y el centro ECP (`check_3c_overlap != 0`).
  - Ejemplos reales probados que abortan con `SIGABRT` (`rc = -6`) en PySCF sin parchear:
    1. **`CuI`** ($R = 2.33\ \text{\AA}$): `Cu` con `cc-pv5z` (todos los electrones, $l_b = 6$) + `I` con `lanl2dz` (base y ECP integrados).
    2. **`AgH`** ($R = 3.06\ \text{Bohr}$) y **`AgI`** ($R = 2.55\ \text{\AA}$): `Ag` con `cc-pv5z-pp` ($l_b = 6$) + ECP `def2-svp`.
  - *(Nota de precisión: si el ECP no contiene canal local $l = -1$, como el ECP `cc-pv5z-pp` Stuttgart-MDF de `Au` cuyos canales son exclusivamente `[0, 1, 2, 3, 4]`, `ECPtype1_cart` no entra en el lazo y `ipnuc` no cae para $l_b = 6, l_c = 4$, aunque `ipipnuc` con $l_b = 5, l_c = 5$ sí cae por el defecto de `ECPscalar_cache_size`.)*

---

## 5. Archivos de evidencia y reproducción

- `environment.txt`: entorno exacto, compilador y SHA de PySCF `master`.
- `repro/make_patch.py` y `repro/minimal.patch`: parche mínimo contra `master 484fd7c`.
- `repro/run_ecp_asan.c`: arnés C compilado con `-fsanitize=address -g`.
- `repro/test_nuc_grad_fd.py`: validación de `mf.nuc_grad_method()` vs FD4 con `libcgto_patched.dylib`.
- `evidence/unpatched_asan_traces.txt`: trazas completas de AddressSanitizer para los 4 modos de fallo.
- `evidence/asan_matrix_comparison.txt`: tabla completa de los 96 casos (sin parchear vs parcheado).
- `evidence/regression_lb_le_5.txt`: verificación bit a bit (`0.000e+00`) para $l_b \le 5$.
- `evidence/nuc_grad_fd_validation.txt`: `SIGABRT` sin parchear vs concordancia a `5.17e-14` (`NaH`) y `1.22e-10` (`AgH`) con FD4 tras aplicar el parche.
