# Cand-05: `simple-dftd3` `gcp.f90` — comentario no verificado (`"12 terms seem accurate enough"`) y cruce de rama sobre `|za - zb| < 0.1` en lugar de `|bx|`

- **Librería:** `simple-dftd3` (`v1.7.0`, commit `84cb7cddd8d6669c4ecd59b06d53e0d373656239`, 2026-10-28)
- **Archivo fuente:** `src/dftd3/gcp.f90` (líneas 806, 1106–1132, 1209, 1303, 1419–1473)
- **Patrón del catálogo:** **Patrón 2 — Comentarios que afirman un hecho numérico no verificado** (+ cambio de rama mal parametrizado)
- **Oráculo:** Consistencia analítica con diferencias finitas de 4 puntos (FD4 con barrido de paso $h$) para $\partial g_z / \partial R$ vs $H_{zz}$, más referencia `mpmath` de 50 dígitos para las integrales auxiliares $B_k(x) = \int_{-1}^{1} t^k e^{-xt}\,dt$
- **Clase:** Precisión numérica (cancelación catastrófica + truncamiento de serie + discontinuidad en $|x| = 10^{-6}$)
- **Nivel de evidencia:** **Defecto confirmado** (reproducido en `simple-dftd3` sin modificar; reducido en 6 órdenes de magnitud — de `4.68e-3` a `1.75e-11` en FD4 — con parche mínimo en copia aislada)
- **Estado upstream:** Presente en `v1.7.0` / `HEAD` (`84cb7cdd`); no consta en `ostra-audit` ni en issues conocidos. Sin reportar (Paso 6).

---

## 1. Detección (Paso 1)

En las rutinas de solapamiento de orbitales de Slater de gCP (`ssovl`, `gsovl`, `dsovl` en `src/dftd3/gcp.f90`), las integrales auxiliares $B_k(x)$ con $x = b_x = (\zeta_B - \zeta_A) R / 2$ se evalúan por dos ramas:
1. **Rama de serie (`bint`, `:1106-1132`, y `bint_all`, `:1445-1473`)**, con el comentario explícito:
   ```fortran
   ! 'infinite' sum is numerically unstable. 12 terms seem
   ! accurate enough
   ```
   que además introduce un corte abrupto `if (abs(x) < 1e-6)` donde $B_1, B_3, B_5, B_7$ se fuerzan exactamente a `0.0` (salto relativo del `100%` en $|x| = 10^{-6}$).
2. **Rama cerrada (`B0..B6`, `g3s3s`, `baux_all` en `:1419-1440`)**, que calcula $e^x S_+ - e^{-x} S_-$ con términos $O(k!/x^{k+1})$.

El criterio de selección (`gcp.f90:806, 1209, 1303`) decide sobre `abs(za - zb) < 0.1_wp` en vez de sobre el argumento real $|b_x| = |\zeta_B - \zeta_A| R / 2$. Esto produce dos regímenes problemáticos:

- **Régimen A (`|za - zb| >= 0.1` a distancias de enlace $R \sim 0.8$–$2.0\ \text{Bohr}$, donde $|b_x| \sim 0.05$–$0.15$):**
  Se invoca la fórmula cerrada `baux_all` / `g3s3s` en zona de cancelación catastrófica (`evidence/unpatched.txt` vs `mpmath` 50 dígitos):
  - `x = 0.05`: $B_6$ tiene error relativo `6.51e-5`, mientras **$B_7$ y $B_8$ cancelan por completo a `0.0000` (error `100%`)**.
  - `x = 0.10`: $B_7$ tiene error `1.61e-3`, **$B_8$ tiene error `1.54%` (`1.54e-2`)**.
  - `x = 0.1218` (par real **`Cl–Ge`** en `r2SCAN-3c` a $R = 2.0\ \text{Bohr}$, $|\Delta\zeta| = 0.1218$): $B_7$ y $B_8$ pierden ~4 dígitos (`2.43e-4` y `2.71e-4`), y `g3s3s` presenta ruido de cancelación en $1/F_b^7$ que contamina las derivadas numéricas del gradiente (`FD4` de $g_z$ vs $H_{zz}$ diverge a `4.682e-3` con $h = 10^{-5}$).
- **Régimen B (`|za - zb| < 0.1` con `bint` / `bint_all` fija a 12 términos):**
  - En $|x| = 9.9 \times 10^{-7}$: salto del `100%` en todos los $B_k$ impares ($B_1, B_3, B_5, B_7$).
  - A distancias grandes (o cutoffs mayores): para `x = 1.49` (p. ej. `Fe–Br` en `r2SCAN-3c` a $R = 30\ \text{Bohr}$, el cutoff por defecto) el error de truncamiento es `5.44e-9`; para `x = 3.0`, sube a `1.22e-5`; para `x = 5.0`, llega a **`1.64e-3`** (refutando la afirmación general de que 12 términos bastan sin acotar $|b_x|$).

## 2. Filtrado de artefactos (Paso 2)

- Se ejecutó un barrido de paso FD4 (`h = 1e-2, 1e-3, 1e-4, 1e-5`) en `Cl–Ge` (`r2SCAN-3c`, $R = 2.0\ \text{Bohr}$) para separar el error analítico de $H_{zz}$ del ruido de cancelación en `gsovl` (`g3s3s`):
  - **Sin parchear (`evidence/unpatched.txt`):** el error FD4 **crece** al reducir $h$ (`1.88e-6` en `h=1e-2` $\to$ `1.36e-5` en `1e-3` $\to$ `2.95e-5` en `1e-4` $\to$ `4.68e-3` en `1e-5`), firma inequívoca de ruido de cancelación de $\sim 10^{-8}$ en `g3s3s`. Al comparar directamente con la versión parcheada (sin FD), el error intrínseco en `Cl–Ge` es `6.67e-6` relativo en $H_{zz}$ a $R = 1.0\ \text{Bohr}$ (`1.72e-6` en $E$ a $R = 0.8\ \text{Bohr}$, `2.55e-7` en $H_{zz}$ a $R = 2.0\ \text{Bohr}$).
  - **Con parche (`evidence/patched.txt`):** en el **mismo** barrido de paso, el error FD4 cae a **`6.73e-12` (`h=1e-3`)**, **`8.62e-12` (`h=1e-4`)** y **`1.75e-11` (`h=1e-5`)** — una mejora de **6 órdenes de magnitud**.

## 3. Cruce con documentación y código fuente (Paso 3)

- `gcp.f90:1207-1208` ya reconoce que las fórmulas cerradas son singulares cuando $\zeta_B \to \zeta_A$ (`! Nearly equal exponents: the closed-form expressions below are singular for zb -> za`), pero el umbral `abs(za - zb) < 0.1` ignora el factor $R/2$ del argumento verdadero $b_x = (\zeta_B - \zeta_A) R / 2$.
- En el rango por defecto `cutoff%gcp = 30 Bohr` y pares con $|\Delta\zeta| < 0.1$, $|b_x| \le 1.5$ mantiene el truncamiento de 12 términos en $\sim 5 \times 10^{-9}$, pero para cutoffs más amplios (o llamadas directas a `bint`/`bint_all`) el truncamiento fijo falla rápidamente, y el corte `abs(x) < 1e-6` nunca fue necesario en la serie de potencias.

## 4. Aislamiento de la causa y parche mínimo (Paso 4)

- **Líneas exactas:** `src/dftd3/gcp.f90:806, 1118-1129, 1209, 1303, 1448-1472`.
- **Parche mínimo verificado en copia aislada (`repro/minimal.patch`):**
  1. Conmutar `ssovl`, `gsovl` y `dsovl` sobre `abs(bx) < 1.5_wp` en vez de `abs(za - zb) < 0.1_wp`.
  2. Eliminar la rama discontinua `abs(x) < 1e-6_wp` y sumar la serie de `bint` / `bint_all` adaptativamente hasta convergencia de máquina (`abs(pw(i)) < epsilon(1.0_wp) * 1e-2_wp`, máx. 60 términos).
- **Resultado con parche (`evidence/patched.txt`):** consistencia FD4 entre gradiente y Hessiano analítico en `Cl–Ge` pasa de `4.68e-3` a **`6.73e-12`–`1.75e-11`** (6 órdenes de magnitud).

## 5. Impacto en observables (Paso 5)

- **Bajo en energías de equilibrio** (`~10^-8`–`10^-6` relativo en pares como `Cl–Ge`, `Mn–Se`, `O–S` en `r2SCAN-3c`).
- **Moderado en frecuencias vibracionales / Hessianos numéricos** calculados diferenciando el gradiente gCP (ruido `~10^-5`–`4.7e-3` en `Cl–Ge` para pasos estándar `h = 1e-4`–`1e-5 Bohr`) y en distancias comprimidas ($R \lesssim 1.0\ \text{Bohr}$, donde $B_7, B_8$ colapsan a `0.0`).

## 6. Reproducibilidad

Ver `repro/repro_bint_baux.f90`, `repro/check_bk_mpmath.py`, `repro/minimal.patch`, `evidence/unpatched.txt` y `evidence/patched.txt`.
