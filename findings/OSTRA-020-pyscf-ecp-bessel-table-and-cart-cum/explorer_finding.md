# Cand-02: Tablas de tamaño fijo desacopladas en `PySCF` `libcgto` (`nr_ecp.c` / `nr_ecp.h` en `master`): OOB en `z == 16.0` y desbordamiento de pila en `omega_nuc[CART_CUM]`

- **Librería:** `PySCF` (`libcgto`: `pyscf/lib/gto/nr_ecp.c`, `pyscf/lib/gto/nr_ecp.h`)
- **Versiones verificadas:** Release instalado `pyscf 2.14.0` **y** rama `master` de GitHub (descargada el 2026-11-22)
- **Patrón del catálogo:** **Patrón 1 — Tablas de tamaño fijo desacopladas del índice máximo** (tercer caso independiente junto a `OSTRA-008` y `Cand-01`; incluye nota de Patrón 2 sobre el nuevo comentario de `master`)
- **Oráculo:** AddressSanitizer (`-fsanitize=address`) + referencia analítica `mpmath` (50 dígitos) para $K_l(z) = e^{-z} i_l(z)$ + evaluación directa `mol.intor("ECPscalar_ipnuc")`
- **Clase:** Error numérico nuevo / acceso fuera de límites (OOB global y de pila)
- **Nivel de evidencia:** **Defecto confirmado** (OOB en `z == 16.0` probado con ASan en código textual de `master` y eliminado con parche mínimo; desbordamiento de pila en `omega_nuc[CART_CUM]` probado en `pyscf 2.14.0` con `nr_ecp.h` idéntico a `master`)
- **Estado upstream:**
  - **Importante (verificación contra `master`):** `_factorial[24]` (que corrompía `ECPsph_ine` desde `order = 12` en `2.14.0`, inyectando un error relativo de `8.54e-5` en `ECPscalar_cart` para `l = 6`) **ya fue ampliado a 65 entradas en `master`**. Por ello NO lo reportamos como defecto vigente en `master`.
  - **Dos defectos de Patrón 1 siguen sin corregir en `master`:** (A) el error por uno `z > 16` en `ECPsph_ine_opt` (`entry = 400` fuera de tablas de 400 filas) y (B) el tamaño fijo `CART_CUM = 456` (`// upto l = 12`) / `_offset_cart[15]` frente a derivadas de capas `l >= 6`. Además, el nuevo comentario en `master` que declara soporte hasta `order = 32` mediante la suma alternante $R_l(-z)/(2z)$ sufre cancelación catastrófica para `order >= 20` cerca de `z = 16..20`.
  - Sin reportar upstream (Paso 6).

---

## 1. Detección (Paso 1)

### Defecto A (vigente en `2.14.0` y `master`): OOB en `z == 16.0` dentro de `ECPsph_ine_opt`
En `pyscf/lib/gto/nr_ecp.h:21-22` se definen `K_TAB_ENTRIES = 400` y `K_TAB_INTERVAL = 16.0 / 400 = 0.04` (cubriendo `[0, 16)`). Sin embargo, `ECPsph_ine_opt` (`nr_ecp.c:4717` en `master`, `:4694` en `2.14.0`) comprueba:
```c
if (z < 1e-7 || z > 16 || order > K_TAB_COL - K_TAYLOR_MAX - 1)
```
Cuando `z == 16.0`, la condición `z > 16` es falsa y entra en la rama tabulada con:
```c
int entry = (int)floor(z / K_TAB_INTERVAL); // entry = 400 (índice válido: 0..399)
```
Esto lee la fila `400` (un bloque completo fuera de los arreglos globales `_sph_ine_tab_order7[400]` para `order <= 7` y `_sph_ine_tab[400]` para `8 <= order <= 16`), cayendo en la memoria estática contigua (`_factorial[]`):

- **Con AddressSanitizer en código textual de `master`** (`evidence/unpatched_asan.txt`):
  - `order = 7, z = 16.0`: `AddressSanitizer: global-buffer-overflow — READ of size 64 ... 0 bytes after global variable '_sph_ine_tab_order7' of size 204800`.
  - `order = 16, z = 16.0`: `AddressSanitizer: global-buffer-overflow — READ of size 144 ... 0 bytes after global variable '_sph_ine_tab' of size 76800`.
- **Sin sanitizador (`-O2`, silencioso)** (`evidence/unpatched_values.txt` vs `mpmath` 50 dígitos):
  - `order = 7, z = 16.0`: `out[0] = -2.500000e-03` (signo invertido! exacto: `+3.125000e-02`), `out[7] = +6.526275e-01` (exacto: `+5.310895e-03`, error `122×`).
  - `order = 16, z = 16.0`: `out[0] = -6.747435e-03` (signo invertido!), `out[16] = +4.848490e-01` (exacto: `+9.474534e-06`, error `51,170×`).

### Defecto B (vigente en `2.14.0` y `master`): Desbordamiento de pila en `omega_nuc[CART_CUM]`
`nr_ecp.h` es **idéntico byte a byte** entre `2.14.0` y `master` (`diff` vacío) y define:
```c
#define CART_CUM (455+1) // upto l = 12
```
Mientras que `type1_rad_ang` (`nr_ecp.c:5253` en `master`) y `type2_facs_ang` (`nr_ecp.c:5078`) declaran `double omega_nuc[CART_CUM]` y escriben hasta `lmax = li + lj` (o `li + 1 + lj` en `ECPscalar_ipnuc`, `li + 2 + lj` en derivadas segundas). Para una base con funciones `i` (`l = 6`, como `cc-pV6Z`) y `mol.intor("ECPscalar_ipnuc")`, `lmax = 6 + 1 + 6 = 13 > 12`, escribiendo `105` doubles (`840` bytes) más allá de `omega_nuc[456]` en la pila.

- **Evidencia observada** (`evidence/unpatched_values.txt`):
  - `mol.intor("ECPscalar_ipnuc")` con `Au` (ECP local) y `H` (`l = 6`) aborta inmediatamente con **`exit code -6` (`SIGABRT`, stack smashing)**!

### Nota adicional (Patrón 2 en el nuevo comentario de `master`)
El commit reciente de `master` que amplió `_factorial` añadió el comentario (`nr_ecp.c:4464`):
```c
// Bessel order can reach 2*ANG_MAX+2 = 32 for second derivatives of a pair of
// l=15 shells; the expansion needs (2*order)!.
```
Sin embargo, evaluar directamente la serie asintótica alternante $R_l(-z)/(2z)$ para `z > 16` hasta `order = 32` sufre **cancelación catastrófica** (evidencia en `evidence/unpatched_values.txt` vs `mpmath`):
- `z = 16.01`: `l = 24` da `+5.595e-09` vs exacto `+1.701e-09` (**`229%` de error**); `l = 28` da `+6.115e-07` vs `+9.720e-12` (**`62,900×`**); `l = 32` da `+1.563e-04` vs `+3.416e-14` (**`4.58e9×`**).
- `z = 20.0`: `l = 32` invierte el signo (`-3.328e-09` vs exacto `+2.136e-12`, error `1,559×`).

## 2. Filtrado de artefactos (Paso 2)

- Verificado contra la rama `master` actual (descartando `_factorial[24]` como defecto de `master`, pues ya se corrigió allí, aunque sigue afectando al release `2.14.0`).
- Verificado con ASan (sin falsos positivos posibles) y referencia analítica de 50 dígitos (`mpmath`).

## 3. Cruce con documentación y código fuente (Paso 3)

- `nr_ecp.h:22` documenta explícitamente `[0, 16]` pero usa `floor(z / 0.04)` en un arreglo de 400 filas (`0..399`), incompatible con incluir `z = 16.0`.
- `nr_ecp.c` en `master` pretende soportar hasta `l = 15` (`order = 32`), pero `CART_CUM = 456` (`// upto l = 12`) y `_offset_cart[15]` (`0..14`) no se actualizaron en paralelo.

## 4. Aislamiento de la causa y parche mínimo (Paso 4)

- **Parche mínimo verificado en copia aislada (`repro/minimal.patch`):**
  1. Cambiar `z > 16` a `z >= 16` en `ECPsph_ine_opt` y `ECPsph_ine`.
  2. Sustituir la suma alternante $R_l(-z)/(2z)$ para `z >= 16` por recurrencia ascendente mientras $2l + 1 \le z$ y recurrencia descendente de Miller por encima (análoga a `libecpint`).
  *(Nota: para el Defecto B en `type1_rad_ang`/`type2_facs_ang`, el arreglo `omega_nuc` debe dimensionarse dinámicamente en `cache` o ampliar `CART_CUM` y `_offset_cart` hasta el `lmax` máximo permitido.)*
- **Resultado con parche (`evidence/patched_asan.txt`):**
  - **Cero errores de AddressSanitizer** en `z = 16.0` (`order = 7` y `order = 16`).
  - Concordancia con `mpmath` de 50 dígitos a `< 2e-15` relativo en todos los órdenes `l = 0..32` para `z = 16.00, 16.01, 17.00, 20.00` (p. ej. `z = 16.01, l = 32`: `+3.4157155895e-14` exacto vs `+1.563e-04` sin parchear).

## 5. Impacto en observables (Paso 5)

- **Defecto A (`z == 16.0`):** ocurre en un punto exacto de medida cero, pero cualquier cuadrícula radial donde $2\alpha R r_i = 16.0$ (p. ej. pruebas analíticas o mallas específicas) lee memoria fuera de límites e invierte el signo de $K_0..K_{10}$.
- **Defecto B (`omega_nuc[CART_CUM]`):** impide calcular gradientes ECP (`ECPscalar_ipnuc`) con funciones `i` (`l = 6`, p. ej. `cc-pV6Z`) — provoca `SIGABRT` inmediato.
- **Nota `2.14.0` (`_factorial[24]`, ya corregido en `master`):** en el release `2.14.0` todavía publicado en PyPI, `ECPscalar_cart` con `l = 6` (`a = 0.7`, `R = 2.5`) sufre un error silencioso de `8.54e-5` relativo (`1.36e-4 Eh`).

## 6. Reproducibilidad

Ver `repro/repro_master_ecp_bessel.c`, `repro/repro_pyscf_stack_and_mpmath.py`, `repro/minimal.patch`, `evidence/unpatched_asan.txt`, `evidence/unpatched_values.txt` y `evidence/patched_asan.txt`.
