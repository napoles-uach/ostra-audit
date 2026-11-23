# Verificación independiente — OSTRA-008 (`libcint` v6.2.1, `src/fmt.c` `TURNOVER_POINT`)

- **Repositorio y commit auditados:** `https://github.com/sunqm/libcint`, tag `v6.2.1`, commit `ab88a545c5799f45a2bcc720551eda2c43039f0f` (clon intacto en [`libcint/`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/libcint)).
- **Entorno de prueba:** macOS arm64 (`Apple clang 21.0.0`, `GCC 16.2.0` con `libquadmath` habilitado, Python 3.13.3 con `mpmath 1.3.0`, `numpy 2.5.1` y `PySCF 2.14.0` con `libcint.6.dylib` v6.1.3).
- **Disciplina de verificación:** Fase 1 ejecutada completamente a ciegas (sin consultar `ostra-audit` ni issues externos); Fase 2 ejecutada solo tras cerrar por completo todos los resultados de Fase 1.

---

## Veredicto y separación de evidencias

| Nivel de evidencia | Estado | Alcance exacto |
|---|---|---|
| **Anomalía reproducida** | **Confirmada** | En las 4 compilaciones probadas (`clang -O2`, `clang -O0`, `gcc-16 -O2 + libquadmath`, y `libcint.6.dylib` v6.1.3 de PySCF), tanto `gamma_inc_like`/`fmt_erfc_like` para $m \ge 40$ como integrales SR reales (`int2e_sph`, `int3c2e_sph`) con $\mathrm{nroots} \in \{21, 22, 23, 24\}$ devuelven errores relativos desde $3\times 10^{-7}$ hasta $10^{+50}$ (y hasta $10^{+85}$ en $F_m(T)$), **en silencio** (sin error de retorno ni aviso en `stderr`) en amplios regímenes de $R$ incluso con `libquadmath` de 128 bits. |
| **Causa identificada** | **Confirmada** | Lectura fuera de límites (`global-buffer-overflow`, `READ of size 8`) sobre `static double TURNOVER_POINT[40]` en [`src/fmt.c:42-83`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/libcint/src/fmt.c#L42-L83) en 6 funciones sin comprobación de rango (`fmt.c:214, 256, 363, 404, 482, 532`). Cuando el valor basura leído en memoria es menor que $T$ (p. ej. subnormales $\sim 10^{-314}$ o `0.0` de `_JACOBI_COEF_ORDER`), se elige la recurrencia ascendente en lugar de la serie descendente para $T < t_0(m)$, provocando cancelación catastrófica. |
| **Defecto confirmado** | **Confirmado (con matices)** | El parche de Fase 2 (extender `TURNOVER_POINT` hasta $m=71$) elimina por completo el `global-buffer-overflow` en ASan para todos los caminos internos de `libcint`, restaura $F_m(T)$ a precisión de máquina ($\le 4.6\times 10^{-16}$) para $T < t_0(m)$, reduce el error de integrales SR con `libquadmath` de $10^{+25}$ a $5\times 10^{-12}\text{–}7\times 10^{-11}$, y preserva **identidad bit a bit (`0 ULP`)** para $m \le 39$ y $\ell$ bajos. |

### Resumen honesto: qué se confirma, qué se matiza y qué contradice o excede la afirmación

- **Confirmado sin reservas:**
  1. `TURNOVER_POINT` tiene exactamente 40 entradas (`m = 0..39`, 320 bytes); toda invocación con $m \ge 40$ (y $t \ne 0$, y $t\cdot\mathrm{lower}^2 \le 200$ en `erfc_like`) realiza una lectura OOB de 8 bytes.
  2. Afecta de forma directa a integrales SR (`omega < 0`) de momento angular alto (`int2e_sph` y `int3c2e_sph`) cuando $\mathrm{nroots} \in \{21, 22, 23, 24\}$ ($m = 2\cdot\mathrm{nroots} - 1 \in \{41, 43, 45, 47\}$), y produce resultados **numéricamente erróneos en silencio** (verificado con `stderr = ''` tanto en `gcc-16` con `libquadmath` —p. ej. 3.77 % de error en `(11,11|10,10)` a $R=0.8$— como en `libcint.6.dylib` v6.1.3 de PySCF).
  3. El resultado depende fuertemente del binario (de qué símbolo coloque el compilador/enlazador tras `TURNOVER_POINT`).

- **Matices y puntos no anticipados en la afirmación original:**
  1. **Alcance en $\ell$ para química cuántica práctica:** en `int2e_sph` sin derivadas se requiere $L_{\mathrm{tot}} = \ell_i+\ell_j+\ell_k+\ell_l \ge 40$ (p. ej. `(10,10|10,10)`, orbitales tipo $n$); incluso con derivadas de 4.º orden (`sum(ng)=4`) se requiere $L_{\mathrm{tot}} \ge 36$ (`(9,9|9,9)`). En `int3c2e_sph` se requiere $\ell_i+\ell_j+\ell_l \ge 40$ (p. ej. `(14,14|12)`). Está muy por encima de bases estándar ($g, h, i$ con $\ell \le 6$), aunque dentro de `ANG_MAX = 15`.
  2. **Techo duro en `CINTsr_rys_roots`:** `src/rys_roots.c:250-255` solo soporta hasta `nroots = 24` ($m = 47$); para `nroots >= 25` imprime error y llama `exit(1)`. Por tanto, en integrales SR los únicos índices OOB alcanzables son los **impares** $m \in \{41, 43, 45, 47\}$.
  3. **También puede afectar a Coulomb completo / LR en el *fallback* de Schmidt:** en `CINTrys_roots` (`src/rys_roots.c:51-53, 113`), si `CINTqrys_jacobi` o `CINTqrys_laguerre` devuelve error para `nroots >= 20`, `segment_solve` recurre a `CINTqrys_schmidt`, que llama a `qgamma_inc_like(..., nroots * 2)` con $m \in \{40, 42, \dots, 62\}$. Observamos dinámicamente este *fallback* en `pyscf_v613` (arm64 sin `quadmath`) con `omega = +0.5` en `l = (10,10|10,10)`, $R = 0.2$ (`rys_roots fails: nroots=21 x=0.00799599`).
  4. **Inestabilidad ortogonal en arm64 sin `libquadmath`:** en macOS arm64 (`sizeof(long double) == 8`, sin `HAVE_QUADMATH_H`), `CINTqrys_jacobi` se mapea a `CINTlrys_jacobi` (64 bits), y `lrys_wheeler` pierde precisión a $\mathrm{nroots} \ge 20$ independientemente del OOB (p. ej. en el control `m=39` a $R=1.0$ da $1.07\times 10^{-5}$ frente a $1.95\times 10^{-12}$ en `gcc-16` con `libquadmath`). Esta inestabilidad persiste tras el parche en builds sin `quadmath` a $R=1.0$, mientras que en `gcc-16 + libquadmath` el parche elimina 100 % del error ($10^{-12}\text{–}10^{-11}$).
  5. **Sobre el parche propuesto (Fase 2):** `m = 71` **no** es el máximo exacto alcanzable (es un margen holgado por encima de 47 en SR, 62 bajo `assert(rys_order < MXRYSROOTS)`, y 63 por el búfer `fmt[MXRYSROOTS*2]`). Además, el parche sin comprobación de límites sigue leyendo OOB si un usuario externo llama a la API exportada (`gamma_inc_like`, etc.) con $m \ge 72$, y la fórmula $t_0(m)$ en precisión `double` aún produce pérdida de dígitos justo por encima de $t_0(m)$ para $m$ muy grande ($1.46\times 10^{-3}$ en $m=71$, $T=27.31$, frente a $<10^{-15}$ en `__float128`).

---

## Fase 1 — Verificación a ciegas

### 1.1 Análisis estático de `src/fmt.c`

- **Tamaño de `TURNOVER_POINT`:** [`src/fmt.c:42-83`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/libcint/src/fmt.c#L42-L83) define `static double TURNOVER_POINT[]` con exactamente **40 elementos** (índices `0..39`, 320 bytes): `m=0, 1` valen `0.`, `m=2` vale `0.866025403784`, ..., `m=39` vale `14.93687515212`.
- **Todas las funciones que indexan `TURNOVER_POINT[m]`** (6 en total, todas indexan directamente con el argumento entero `m`, sin ningún control de rango):

| Función | Archivo y línea | Expresión exacta | Condición de lectura OOB |
|---|---|---|---|
| `gamma_inc_like` | [`src/fmt.c:214`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/libcint/src/fmt.c#L214) | `else if (t < TURNOVER_POINT[m])` | `t != 0` y **`m >= 40`** |
| `lgamma_inc_like` | [`src/fmt.c:256`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/libcint/src/fmt.c#L256) | `else if (t < TURNOVER_POINT[m])` | `t != 0` y **`m >= 40`** |
| `fmt_erfc_like` | [`src/fmt.c:363`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/libcint/src/fmt.c#L363) | `else if (t < TURNOVER_POINT[m])` | `lower != 0`, `t != 0`, `t*lower^2 <= 200` y **`m >= 40`** (si `lower == 0` delega a `gamma_inc_like` en L343) |
| `fmt_lerfc_like` | [`src/fmt.c:404`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/libcint/src/fmt.c#L404) | `else if (t < TURNOVER_POINT[m])` | `lower != 0`, `t != 0`, `t*lower^2 <= 200` y **`m >= 40`** (si `lower == 0` delega a `lgamma_inc_like` en L384) |
| `qgamma_inc_like` (`HAVE_QUADMATH_H`) | [`src/fmt.c:482`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/libcint/src/fmt.c#L482) | `else if (t < TURNOVER_POINT[m])` | `t != 0` y **`m >= 40`** |
| `fmt_qerfc_like` (`HAVE_QUADMATH_H`) | [`src/fmt.c:532`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/libcint/src/fmt.c#L532) | `else if (t < TURNOVER_POINT[m])` | `lower != 0`, `t != 0`, `t*lower^2 <= 200` y **`m >= 40`** (si `lower == 0` delega a `qgamma_inc_like` en L512) |

### 1.2 Relación entre $m$, `nroots` y momento angular $\ell$

1. **Sitios de llamada a `fmt.c` en todo `libcint/src/`:**
   - **Ruta Wheeler–Jacobi (camino principal para SR con $\mathrm{nroots} \ge 7$):** [`src/rys_wheeler.c:3678, 3703, 6334`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/libcint/src/rys_wheeler.c#L3678-L3715) (`CINT*rys_jacobi`). Para `lower != 0` llama siempre a `*naive_jacobi_moments(nroots * 2, x, lower, moments)`, que a su vez invoca `fmt_*erfc_like(fmt, t, lower, n - 1)` ([`rys_wheeler.c:3346, 3538, 6196`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/libcint/src/rys_wheeler.c#L3346)):
     $$m = 2\cdot\mathrm{nroots} - 1.$$
     Por tanto $m \ge 40 \iff \mathrm{nroots} \ge 21$. En [`src/rys_roots.c:244-255`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/libcint/src/rys_roots.c#L244-L255), `CINTsr_rys_roots` soporta `nroots = 21..24` usando `segment_solve1(..., lower_bp1=0.25, breakpoint=60, CINTqrys_jacobi, CINTqrys_laguerre, CINTqrys_jacobi)` (nota: `lower_bp2` no se usa en `segment_solve1`, [`rys_roots.c:128-143`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/libcint/src/rys_roots.c#L128-L143)), por lo que **siempre** entra en `CINTqrys_jacobi` cuando `lower >= 0.25` o cuando `lower < 0.25` y `x <= 60`. Para `nroots >= 25` entra en `default:` y aborta con `exit(1)`.
     - **Índices concretos leídos OOB por SR:** `TURNOVER_POINT[41]` (`nroots=21`), `[43]` (`nroots=22`), `[45]` (`nroots=23`), `[47]` (`nroots=24`).
   - **Ruta Schmidt (*fallback* de `CINTrys_roots` para Coulomb/LR y API directa):** [`src/rys_roots.c:1762, 1764, 1863, 1865, 1983, 1985`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/libcint/src/rys_roots.c#L1758-L1986) llama a `*gamma_inc_like` / `fmt_*erfc_like` con
     $$m = 2\cdot\mathrm{nroots}.$$
     Aquí el OOB comienza en `nroots = 20` (`TURNOVER_POINT[40, 42, ..., 62]`). En `CINTrys_roots` (`lower == 0`), el camino principal usa `flocke_jacobi_moments` (algoritmo de Miller, no llama a `fmt.c`, pues `x <= SMALLX_LIMIT = 3e-7` se ataja antes en [`rys_roots.c:59`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/libcint/src/rys_roots.c#L59)), pero si `fn1`/`fn2` falla, `segment_solve` cae en `CINTqrys_schmidt` ([`rys_roots.c:51-53`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/libcint/src/rys_roots.c#L51-L53)).

2. **Combinaciones de $\ell$ que disparan $m \ge 40$ en SR:**
   En [`src/g2e.c:69-81`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/libcint/src/g2e.c#L69-L81), $\mathrm{nroots} = \mathrm{rys\_order} = \lfloor L_{\mathrm{ceil}}/2 \rfloor + 1$ (para $\mathrm{rys\_order} > 3$), con $L_{\mathrm{ceil}} = \sum (\ell_a + \mathrm{ng}_a)$. Por ello $\mathrm{nroots} \ge 21 \iff L_{\mathrm{ceil}} \ge 40$:

| Integral | Condición en $L_{\mathrm{tot}}$ | Ejemplos concretos de $\ell$ | `nroots` | Índice OOB `TURNOVER_POINT[m]` |
|---|---|---|---|---|
| `int2e` (control) | $L_{\mathrm{tot}} = 38, 39$ | `(10,10|10,9)` | 20 | `m = 39` (en rango) |
| `int2e` (`ng=0`) | $L_{\mathrm{tot}} = 40, 41$ | `(10,10|10,10)`, `(14,14|6,6)` | 21 | **`m = 41`** |
| `int2e` (`ng=0`) | $L_{\mathrm{tot}} = 42, 43$ | `(11,11|10,10)` | 22 | **`m = 43`** |
| `int2e` (`ng=0`) | $L_{\mathrm{tot}} = 44, 45$ | `(12,11|11,10)` | 23 | **`m = 45`** |
| `int2e` (`ng=0`) | $L_{\mathrm{tot}} = 46, 47$ | `(12,12|11,11)` | 24 | **`m = 47`** |
| `int2e_ip1` (`sum(ng)=1`) | $L_{\mathrm{tot}} \ge 39$ | `(10,10|10,9)` | 21 | **`m = 41`** |
| `int2e` deriv. 4.º orden (`sum(ng)=4`) | $L_{\mathrm{tot}} \ge 36$ | `(9,9|9,9)` | 21 | **`m = 41`** |
| `int3c2e` ([`g3c2e.c:70`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/libcint/src/g3c2e.c#L70)) | $\ell_i+\ell_j+\ell_l \ge 40$ | `(14,14|12)` | 21 | **`m = 41`** |
| `int2c2e`, `int1e_grids` | máx. $2\cdot 14 + \mathrm{ng} \le 33$ | — | $\le 17$ | No alcanza $m \ge 40$ |

### 1.3 Confirmación dinámica con AddressSanitizer

Compilamos `libcint v6.2.1` sin modificar con `clang -O1 -g -fsanitize=address -fno-omit-frame-pointer` (ver [`scripts/build_libcint.py`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/scripts/build_libcint.py) y [`scripts/test_asan.c`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/scripts/test_asan.c)). Log completo en [`outputs/asan_phase1.log`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/outputs/asan_phase1.log):

- **Casos en rango / sin `fmt.c` (`EXIT=0`, ASan limpio):** `boys_ok` ($m=39$), `eri_sr_l9` (`int2e_sph`, `(10,10|10,9)`, `omega=-0.5`, $m=39$), `eri_full_l10` (`(10,10|10,10)`, `omega=0`), `eri_lr_l10` (`(10,10|10,10)`, `omega=+0.5`).
- **Casos OOB (`EXIT=134`, `AddressSanitizer: global-buffer-overflow`):** `boys_gamma_40` (`fmt.c:214`), `boys_lgamma_40` (`fmt.c:256`), `boys_erfc_40` (`fmt.c:363`), `boys_lerfc_41` (`fmt.c:404`), `eri_sr_l10` (`int2e_sph` `(10,10|10,10)` `omega=-0.5`), y `eri_sr_3c2e` (`int3c2e_sph` `(14,14|12)` `omega=-0.5`).

Mensaje exacto de ASan para `gamma_inc_like(..., m=40)` ([`outputs/asan_phase1.log:6-13`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/outputs/asan_phase1.log#L6-L13)):
```text
==24605==ERROR: AddressSanitizer: global-buffer-overflow on address 0x000104b52b40 at pc 0x0001049754e0 bp 0x00016baf1b60 sp 0x00016baf1b58
READ of size 8 at 0x000104b52b40 thread T0
    #0 0x0001049754dc in gamma_inc_like fmt.c:214
    #1 0x00010430cb6c in main test_asan.c:135

0x000104b52b40 is located 0 bytes after global variable 'TURNOVER_POINT' defined in '.../libcint/src/fmt.c' (0x000104b52a00) of size 320
SUMMARY: AddressSanitizer: global-buffer-overflow fmt.c:214 in gamma_inc_like
```

Mensaje exacto de ASan para la integral SR real `int2e_sph` con `l=(10,10|10,10)`, `omega=-0.5` ([`outputs/asan_phase1.log:182-196`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/outputs/asan_phase1.log#L182-L196)):
```text
==24635==ERROR: AddressSanitizer: global-buffer-overflow on address 0x000104e6ab48 at pc 0x000104c8e05c bp 0x00016b77b600 sp 0x00016b77b5f8
READ of size 8 at 0x000104e6ab48 thread T0
    #0 0x000104c8e058 in fmt_lerfc_like fmt.c:404
    #1 0x000104c903b8 in CINTlrys_jacobi rys_wheeler.c:3712
    #2 0x000104ca37a0 in CINTsr_rys_roots rys_roots.c
    #3 0x000104c7fe60 in CINTg0_2e g2e.c:4465
    #4 0x000104c34654 in CINT2e_loop_nopt cint2e.c:248
    #5 0x000104c3d890 in CINT2e_drv cint2e.c:841
    #6 0x000104c41ff0 in int2e_sph cint2e.c:1196

0x000104e6ab48 is located 8 bytes after global variable 'TURNOVER_POINT' defined in '.../libcint/src/fmt.c' (0x000104e6aa00) of size 320
SUMMARY: AddressSanitizer: global-buffer-overflow fmt.c:404 in fmt_lerfc_like
```
(8 bytes después de `TURNOVER_POINT[0..39]` corresponde exactamente a `TURNOVER_POINT[41]`, ya que $2\times 21 - 1 = 41$.)

### 1.4 Oráculo físico (`mpmath` a 60 dígitos e identidad exacta $\mathrm{SR} = \mathrm{full} - \mathrm{LR}$) y 1.5 Dependencia del binario

Probamos **cuatro compilaciones distintas** ([`scripts/phase1_oracle_and_builds.py`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/scripts/phase1_oracle_and_builds.py), [`scripts/phase1_deep_check.py`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/scripts/phase1_deep_check.py)):
1. `clang_O2_v621`: `Apple clang 21.0.0 -O2` (`libcint v6.2.1`, sin `quadmath`).
2. `clang_O0_v621`: `Apple clang 21.0.0 -O0` (`libcint v6.2.1`, sin `quadmath`).
3. `gcc16_O2_quad_v621`: `GCC 16.2.0 -O2` (`libcint v6.2.1`, **con `libquadmath` de 128 bits**).
4. `pyscf_v613`: binario empaquetado de PySCF 2.14.0 (`libcint.6.dylib`, `v6.1.3`).

#### Qué hay en memoria justo después de `TURNOVER_POINT` en cada binario

Inspeccionando el mapa de símbolos (`nm -n`) y leyendo los bytes exactos IEEE-754 de cada Mach-O ([`outputs/phase1_deep_check.log`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/outputs/phase1_deep_check.log)):

| `m` | Símbolo/valor en `clang_O2_v621` y `pyscf_v613` (`_JACOBI_COEF_ORDER`) | `clang_O0_v621` (`_JACOBI_ALPHA`) | `gcc16_O2_quad_v621` (constantes + `_qJACOBI_SN`) | Valor verdadero $t_0(m)$ |
|---:|---:|---:|---:|---:|
| 39 (en rango) | `14.93688` | `14.93688` | `14.93688` | `14.93688` |
| **40 (OOB)** | `0.000000e+00` | `3.333333e-01` | `8.862269e-01` | `15.30898` |
| **41 (OOB)** | `4.940656e-324` | `5.238095e-01` | `2.000000e-20` | `15.68098` |
| **42 (OOB)** | `4.243992e-314` | `5.064935e-01` | `-3.105036e+231` | `16.05289` |
| **43 (OOB)** | `2.121996e-314` | `5.030303e-01` | `1.301833e-02` | `16.42470` |
| **45 (OOB)** | `2.121996e-314` | `5.011442e-01` | `1.923278e+00` | `17.16808` |
| **47 (OOB)** | `1.482197e-323` | `5.005974e-01` | `1.937500e+00` | `17.91115` |

Cada binario lee una constante diferente dependiendo del compilador y del nivel de optimización.

#### Oráculo de Boys $F_m(T)$ (`gamma_inc_like`) frente a `mpmath` (error relativo máx. en $k=0\dots m$)

Datos completos en [`outputs/phase1_oracle.log`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/outputs/phase1_oracle.log) (y comportamiento análogo en `fmt_erfc_like`):

| $m$ | $T$ | `clang_O2_v621` | `clang_O0_v621` | `gcc16_O2_quad_v621` | `pyscf_v613` | Causa de divergencia entre binarios |
|---:|---:|---:|---:|---:|---:|---|
| 39 (ctrl) | 0.5 | `2.17e-16` | `2.17e-16` | `2.17e-16` | `2.17e-16` | En rango ($T < 14.94$) |
| 39 (ctrl) | 14.0 | `2.39e-16` | `2.39e-16` | `3.44e-16` | `2.39e-16` | En rango ($T < 14.94$) |
| **40** | **0.5** | **`9.51e+44`** | **`9.51e+44`** | **`2.17e-16`** | **`9.51e+44`** | En `gcc16` lee `0.886 > 0.5` (rama desc.); en los otros `< 0.5` |
| **41** | **0.5** | **`7.90e+46`** | **`2.17e-16`** | **`7.90e+46`** | **`7.90e+46`** | En `clang_O0` lee `0.524 > 0.5` (rama desc.); en los otros `< 0.5` |
| **41** | **2.0** | `4.60e+22` | `4.60e+22` | `4.12e+22` | `4.60e+22` | Todos `< 2.0` → recurrencia ascendente |
| **41** | **10.0** | `4.18e-03` | `4.18e-03` | `4.16e-03` | `4.18e-03` | Todos `< 10.0` (verdadero $t_0=15.68$) |
| **45** | **0.5** | **`4.73e+54`** | **`2.17e-16`** | **`2.17e-16`** | **`4.73e+54`** | `clang_O0` (`0.501`) y `gcc16` (`1.923`) son $> 0.5$ |
| **47** | **5.0** | `6.11e+13` | `6.11e+13` | `5.49e+13` | `6.11e+13` | Todos $< 5.0$ (verdadero $t_0=17.91$) |
| **47** | **10.0** | `3.59e+01` | `3.59e+01` | `3.57e+01` | `3.59e+01` | Error relativo del 3587 % |

#### Oráculo de integrales SR de alto $\ell$ (`int2e_sph`, $\omega = -0.5$) frente a la identidad exacta $\mathrm{full}(\omega=0) - \mathrm{LR}(\omega=+0.5)$

Parámetros: dos centros separados por $R$ con exponentes `(1.2, 0.9 | 1.1, 0.8)` ($a_0 = 0.9975$, $x = a_0 R^2$, $\mathrm{lower} = \sqrt{\theta} = 0.448$). Error relativo en norma $\ell_2$ respecto a `full - LR` de `gcc16` ([`outputs/phase1_deep_check.log`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/outputs/phase1_deep_check.log)):

| Caso ($\ell$, `nroots`, $m$) | $R$ | $x$ | `clang_O2_v621` | `clang_O0_v621` | `gcc16_O2_quad_v621` | `pyscf_v613` |
|---|---:|---:|---:|---:|---:|---:|
| Control `(10,10\|10,9)` (`n=20, m=39`) | 0.2 | 0.040 | `3.00e-12` | `3.04e-12` | `2.02e-12` | `3.00e-12` |
| Control `(10,10\|10,9)` (`n=20, m=39`) | 0.6 | 0.359 | `1.10e-11` | `1.10e-11` | `2.46e-12` | `1.10e-11` |
| Control `(10,10\|10,9)` (`n=20, m=39`) | 0.8 | 0.638 | `6.25e-11` | `6.25e-11` | `1.85e-12` | `6.25e-11` |
| **OOB `(10,10\|10,10)` (`n=21, m=41`)** | **0.2** | **0.040** | **`1.73e+40`** | **`4.83e-12`** | **`2.40e+22`** *(silencioso)* | **`1.73e+40`** |
| **OOB `(10,10\|10,10)` (`n=21, m=41`)** | **0.6** | **0.359** | **`1.51e+21`** | **`6.57e-12`** | **`9.46e+01`** | **`1.51e+21`** |
| **OOB `(10,10\|10,10)` (`n=21, m=41`)** | **0.8** | **0.638** | **`3.24e+16`** | **`3.24e+16`** | **`3.30e-04`** *(silencioso)* | **`3.24e+16`** *(silencioso)* |
| **OOB `(10,10\|10,10)` (`n=21, m=41`)** | **1.0** | **0.997** | **`4.37e+11`** | **`4.37e+11`** | **`3.02e-07`** *(silencioso)* | **`4.37e+11`** *(silencioso)* |
| **OOB `(11,11\|10,10)` (`n=22, m=43`)** | **0.8** | **0.638** | **`3.69e+18`** | **`3.69e+18`** | **`3.77e-02`** *(silencioso)* | **`3.69e+18`** *(silencioso)* |
| **OOB `(11,11\|10,10)` (`n=22, m=43`)** | **1.0** | **0.997** | **`3.19e+13`** | **`3.19e+13`** | **`2.21e-05`** *(silencioso)* | **`3.19e+13`** *(silencioso)* |
| **OOB `(12,12\|11,11)` (`n=24, m=47`)** | **0.2** | **0.040** | **`1.42e+50`** | **`3.67e-11`** | **`3.68e-11`** | **`1.42e+50`** |
| **OOB `(12,12\|11,11)` (`n=24, m=47`)** | **0.8** | **0.638** | **`6.50e+22`** | **`6.50e+22`** | **`3.98e-11`** | **`6.50e+22`** |

Nótese la demostración directa de **dependencia del binario en la propia integral SR**:
- En `(10,10|10,10)` a $R = 0.6$ ($x = 0.359$), `clang_O0` lee `TURNOVER_POINT[41] = 0.5238 > 0.359` y acierta (`6.57e-12`), mientras `clang_O2`/`pyscf_v613` (`4.94e-324`) y `gcc16` (`2.0e-20`) fallan catastróficamente (`1.51e+21` y `9.46e+01`); al pasar a $R=0.8$ ($x=0.638 > 0.5238$), `clang_O0` también salta abruptamente a `3.24e+16`.
- En `(12,12|11,11)`, `gcc16` lee accidentalmente `1.9375` (de `_qJACOBI_SN`) y acierta hasta $x < 1.94$, mientras `clang_O2`/`pyscf` fallan por $10^{+50}$ y `clang_O0` falla a partir de $x > 0.5006$.
- La verificación de `stderr` a nivel de descriptor de archivo C (`fd 2`) confirmó que en `gcc16` (`R=0.2, 0.8, 1.0`) y `pyscf_v613` (`R=0.8, 1.0`), **no se emite ningún aviso ni error** (`stderr == ''`), devolviendo números plausibles erróneos en completo silencio (p. ej. 3.77 % de error en `gcc16` para `(11,11|10,10)` a $R=0.8$).

---

## Fase 2 — Evaluación del parche propuesto (`m = 0..71`)

### (a) ¿Es `m = 71` el máximo alcanzable con `ANG_MAX` de `libcint`?

**No es el máximo exacto alcanzable; es una cota superior conservadora (suficiente para todos los caminos internos, pero mayor que los techos reales del código).**

Desglose exacto desde el código de `libcint v6.2.1`:
1. [`include/cint.h.in:123-126`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/libcint/include/cint.h.in#L123-L126): `ANG_MAX = 15` (comentario `l = 0..15`, `LMAX1 = 16`, `CART_MAX = 136`) y `MXRYSROOTS = 32`.
2. **Integrales SR (`CINTsr_rys_roots`, [`src/rys_roots.c:247-255`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/libcint/src/rys_roots.c#L247-L255)):** techo duro `nroots <= 24` (`nroots >= 25` ejecuta `exit(1)`). Además `JACOBI_ALPHA`/`lJACOBI_ALPHA` tienen 48 entradas ($2\times 24$). Por tanto:
   $$m_{\max}^{\mathrm{SR}} = 2 \times 24 - 1 = \mathbf{47}.$$
3. **Integrales respetando las aserciones de `g2e.c` ([`src/g2e.c:87-99`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/libcint/src/g2e.c#L87-L99)):** `assert(i_l < ANG_MAX)` ($\ell \le 14$) y `assert(rys_order < MXRYSROOTS)` (`nroots <= 31`). Con $\ell = 14$ en 4 centros y el máximo incremento de derivadas 4c2e en `src/autocode/` (`sum(ng[:4]) = 5` en `breit1.c`/`intor4.c`), $\mathrm{rys\_order} = \lfloor (56 + 5)/2 \rfloor + 1 = 31$. Por tanto:
   - vía `CINTqrys_jacobi`: $m_{\max} = 2\times 31 - 1 = \mathbf{61}$;
   - vía `CINTqrys_schmidt`: $m_{\max} = 2\times 31 = \mathbf{62}$.
4. **Techo de los búferes de pila (`MXRYSROOTS = 32`):** en `*naive_jacobi_moments` y `CINTrys_schmidt`, el búfer de Boys tiene tamaño `MXRYSROOTS * 2 = 64` ([`rys_wheeler.c:3342, 3534, 6192`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/libcint/src/rys_wheeler.c#L3342), [`rys_roots.c:1760`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/libcint/src/rys_roots.c#L1760)), por lo que $m = 63$ es el último índice que cabe sin desbordar la pila (`m = 64` ya desborda `fmt_ints[64]`).
5. **Incluso ignorando `MXRYSROOTS` y `-DNDEBUG`** con $\ell = 15$ en los 4 centros y `sum(ng[:4]) = 5`: $\mathrm{rys\_order} = \lfloor 65/2 \rfloor + 1 = 33 \implies m = 65$ (Jacobi) o $66$ (Schmidt), aún inferior a 71 (y ya desbordaría `u[MXRYSROOTS]` antes de llegar a `fmt.c`).
6. **Observación de robustez:** aunque 72 entradas (`0..71`) cubren con holgura todos los caminos internos, `gamma_inc_like`/`fmt_erfc_like` son símbolos públicos exportados y **siguen careciendo de verificación de rango**; en nuestras pruebas ASan con el parche, `gamma_inc_like(f, 10.0, 72)` sigue provocando `global-buffer-overflow` (`EXIT=-6`). Es recomendable añadir además un guardia o fórmula asintótica para `m >= sizeof(TURNOVER_POINT)/sizeof(double)`.

### (b) ¿Reproduce la fórmula las 40 entradas existentes?

Evaluando $t_0(m) = \frac{1}{2}\left(\frac{(2m-1)!!}{\sqrt{2m-1}}\right)^{1/(m-1)}$ a 60 dígitos en `mpmath` ([`scripts/phase2_eval_patch.py`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/scripts/phase2_eval_patch.py), [`outputs/phase2_eval.log`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/outputs/phase2_eval.log)):
- **`m = 0, 1`:** la fórmula no está definida en $m=0$ ni en $m=1$ (división por $m-1=0$); la tabla los fija en `0.` por convención.
- **`m = 2..39`:** **sí, coincide con todas las entradas existentes hasta el truncamiento/redondeo decimal de la tabla** (12–13 cifras significativas en `src/fmt.c`):
  - **Desviación absoluta máxima:** **$8.245 \times 10^{-12}$** en $m = 27$ (tabla `10.46165248270` vs. fórmula `10.4616524827082`).
  - **Desviación relativa máxima:** **$7.881 \times 10^{-13}$** en $m = 27$.
  (Todas las entradas de $m=2\dots 39$ coinciden con la fórmula truncada a la última cifra impresa en `src/fmt.c`.)

### (c) Resultados con el parche aplicado en copia aislada

Aplicamos el parche en [`patched/fmt.c`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/patched/fmt.c) (72 entradas, `m=0..71`; clon original sin tocar) y compilamos `build-patched-asan`, `build-patched-o2` (`clang -O2`) y `build-patched-gcc` (`gcc-16 -O2 + libquadmath`). Resultados completos en [`outputs/phase2_eval.log`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/outputs/phase2_eval.log) y [`outputs/asan_phase2_patched.log`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/outputs/asan_phase2_patched.log):

1. **AddressSanitizer limpio:**
   - Todos los casos previamente fallidos pasan con **`EXIT=0`**: `boys_gamma_40`, `boys_lgamma_40`, `boys_erfc_40`, `boys_lerfc_41`, `boys_gamma_71`, `eri_sr_l10` (`(10,10|10,10)`), y `eri_sr_3c2e` (`(14,14|12)`).
   - `boys_gamma_72` dispara ASan (`EXIT=-6`), confirmando que falta comprobación de límites para llamadas externas con $m \ge 72$.

2. **Función de Boys $F_m(T)$ (`gamma_inc_like`) parcheada frente a `mpmath`:**
   - Para todo $T < t_0(m)$ en $m \in \{40, 41, 45, 47, 50, 60, 71\}$, el error cae de $10^{+21}\text{–}10^{+85}$ a **precisión de máquina ($\le 4.62 \times 10^{-16}$)**.
   - **Matiz relevante en `double` (`fp64`) justo por encima de $t_0(m)$:** la fórmula de `fmt.c:40` aún permite pérdida creciente de dígitos en la recurrencia ascendente de 64 bits inmediatamente arriba de $t_0(m)$: a $T = t_0 + 0.5$, el error en `double` es $4.4\times 10^{-9}$ ($m=40$), $8.9\times 10^{-8}$ ($m=47$), $7.5\times 10^{-6}$ ($m=60$) y **$1.46\times 10^{-3}$** ($m=71$). En `__float128` (`qgamma_inc_like`/`fmt_qerfc_like`, usado por `CINTqrys_jacobi` en integrales SR), el margen de 34 dígitos absorbe esta cancelación sin afectar la precisión `double` final.

3. **Integrales SR de alto $\ell$ (`int2e_sph`, $\omega = -0.5$) parcheadas frente a `full - LR`:**

| Caso | $R$ | $x$ | `clang` orig | `clang` parche | `gcc16+quad` orig | `gcc16+quad` parche |
|---|---:|---:|---:|---:|---:|---:|
| `(10,10\|10,10)` (`m=41`) | 0.2 | 0.040 | `1.73e+40` | **`4.81e-12`** | `2.40e+22` | **`4.96e-12`** |
| `(10,10\|10,10)` (`m=41`) | 0.6 | 0.359 | `1.51e+21` | **`6.62e-12`** | `9.46e+01` | **`5.22e-12`** |
| `(10,10\|10,10)` (`m=41`) | 0.8 | 0.638 | `3.24e+16` | **`2.55e-11`** | `3.30e-04` | **`5.88e-12`** |
| `(10,10\|10,10)` (`m=41`) | 1.0 | 0.997 | `4.37e+11` | `2.84e-05`\* | `3.02e-07` | **`6.67e-12`** |
| `(11,11\|10,10)` (`m=43`) | 0.6 | 0.359 | `3.06e+23` | **`1.15e-11`** | `1.91e+04` | **`9.63e-12`** |
| `(11,11\|10,10)` (`m=43`) | 0.8 | 0.638 | `3.69e+18` | **`4.17e-11`** | `3.77e-02` | **`1.18e-11`** |
| `(11,11\|10,10)` (`m=43`) | 1.0 | 0.997 | `3.19e+13` | `1.60e-03`\* | `2.21e-05` | **`9.60e-12`** |
| `(12,11\|11,10)` (`m=45`) | 0.2 | 0.040 | `8.98e+46` | **`1.02e-11`** | `1.15e-11` | **`1.15e-11`** |
| `(12,12\|11,11)` (`m=47`) | 0.2 | 0.040 | `1.42e+50` | **`3.75e-11`** | `3.68e-11` | **`3.68e-11`** |

*\*Nota: el error residual en `clang` (arm64 sin `libquadmath`) a $R=1.0$ no proviene de `TURNOVER_POINT`, sino de la singularidad en `lrys_wheeler` al ejecutarse en 64 bits (`sizeof(long double) == 8`), fenómeno ya presente en el control en rango `(10,10|10,9)` (`1.07e-5`). Con `libquadmath` (`gcc16+quad`), todas las integrales SR parcheadas coinciden con `full - LR` a $5\times 10^{-12}\text{–}7\times 10^{-11}$.*

4. **Identidad bit a bit para $\ell$ bajo y $m \le 39$:**
   - `gamma_inc_like` y `fmt_erfc_like` en $m = 0\dots 39$ (10 valores de $T \in [0, 50]$): **`bit_identical = True`, `max_abs = 0.0`** en tanto `clang -O2` como `gcc16 + quadmath`.
   - `int2e_sph` para $\ell \in \{(0,0|0,0), (1,1|1,1), (2,2|2,2), (3,3|3,3), (4,4|4,4), (5,5|5,5), (10,10|10,9)\}$, $\omega \in \{0, +0.5, -0.5\}$, $R \in \{0.5, 1.5, 3.5\}$: **`bit_identical = True`, `max_abs = 0.0`, 0 discrepancias** (comparación exacta por bits `uint64`) en ambos compiladores.

---

## Inventario de scripts y salidas en `verif-OSTRA-008/`

- Scripts: [`build_libcint.py`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/scripts/build_libcint.py), [`test_asan.c`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/scripts/test_asan.c), [`phase1_oracle_and_builds.py`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/scripts/phase1_oracle_and_builds.py), [`phase1_deep_check.py`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/scripts/phase1_deep_check.py), [`phase2_eval_patch.py`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/scripts/phase2_eval_patch.py).
- Parche aislado: [`patched/fmt.c`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/patched/fmt.c).
- Salidas crudas: [`asan_phase1.log`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/outputs/asan_phase1.log), [`phase1_oracle.log`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/outputs/phase1_oracle.log), [`phase1_deep_check.log`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/outputs/phase1_deep_check.log), [`phase2_eval.log`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/outputs/phase2_eval.log), [`asan_phase2_patched.log`](file:///Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008/outputs/asan_phase2_patched.log), más JSON estructurados en `outputs/`.
