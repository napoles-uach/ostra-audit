# Auditoría de umbrales y cambios de rama en `libxc` (`devel`, commit `3c7794c50`)

## Entorno de medición y versiones comparadas
- **Versión auditada (`devel`):** clon local `/Users/napolesd/Documents/libxc-devel` (espejo de solo lectura en `libxc-ref/`, commit `3c7794c50`, cadena `xc_version_string() = "7.1.2"`, 711 IDs únicos / 727 nombres), compilada fuera del clon en `build-libxc/libxc.dylib` con `-std=c99 -O2 -fPIC -ffp-contract=off -DHAVE_CBRT -DXC_MAXORDER=2 -DXC_ENFORCE_FERMI_HOLE_CURVATURE`.
- **Versión de comparación (PySCF 2.14.0):** `.venv/lib/python3.13/site-packages/pyscf/lib/deps/lib/libxc.15.dylib` (nota medida: la librería empaquetada en este entorno reporta `xc_version_string() = "7.0.0"`, con 679 IDs únicos).
- **Grupo de control ("la cara"):** Los cambios de rama en **PBE / PBE0** (`m_saturate` en `y = 1`), **SCAN** (`scan_f_alpha` en $\alpha = 0.9818, 1, 1.0217$ y `scan_gx` en $s = 0.01884$), **r²SCAN** (`r2scan_f_alpha` en $\alpha = 0$ y $2.5$), **B3LYP / B88 / PW91** (`x_asinh_cx` en $p = 1/16$) y **ωB97X** (`attenuation_erf` en $a = 1/16, 1/2, 2^{24}$) son continuos hasta precisión de máquina ($10^{-12}$–$10^{-16}$ en valor y primera derivada; solo `r2scan_f_alpha` en $\alpha = 2.5$ muestra un residuo de $\sim 3 \times 10^{-9}$ en `vtau` por el redondeo publicado de los coeficientes polinómicos).

---

## Lista corta de candidatos que superan los pasos 1 a 3 del protocolo

### Candidato 1 — Regresión numérica nueva en `xc_E1_scaled_jet` (`z = 18.0`)
- **Clase tentativa:** Error numérico nuevo (regresión en `devel` / serie `7.1.x` respecto a `7.0.0`).
- **Ubicación exacta:** `src/expint_e1.c:140` (`if (z > 18.0)`), invocado desde `src/maple2c/gga_exc/gga_x_wpbeh.c:3383,3394`.
- **Funcionales afectados (prioridad alta):** `hyb_gga_xc_hse06` (428), `hyb_gga_xc_hse03` (427), `hyb_gga_xc_hse12` (479), `hyb_gga_xc_hse12s` (480), `hyb_gga_xc_lc_wpbe` (478), `gga_x_wpbeh` (524); también invocado en `gga_c_ft97` (114) y en derivadas de orden $\ge 3$ de `lda_x_1d_exponential` (600).
- **Causa aislada (Paso 4):** Para $z > 18$, `xc_E1_scaled_jet` sustituye la recurrencia exacta a partir del ajuste Chebyshev de 16 dígitos `xc_expint_e1_impl(z, 1)` por el truncamiento óptimo de la serie asintótica divergente de $e^z E_1(z)$. A diferencia de `xc_erfcx_jet` en $z = 8$ (resto $\sim e^{-64} \approx 1.6 \times 10^{-28}$), el resto superasintótico mínimo de $e^z E_1(z)$ vale $\sim \sqrt{2\pi/z}\,e^{-z} \approx 1.52 \times 10^{-8}$ en $z = 18$. Como `exc` evalúa `xc_E1_scaled(aux5)` pero `vxc` evalúa `xc_E1_scaled_d1(aux5)`, la energía y el potencial usan aproximaciones incompatibles en $18 < \text{aux5} < 37$.
- **Evidencia medida (Pasos 1–3 y 5):**
  - `xc_E1_scaled_d1(18 ± δ)` salta de `-2.7906111529e-3` a `-2.7906065711e-3` (**salto relativo `1.642e-6`** idéntico en $\delta = 10^{-9}$ y $10^{-14}$).
  - En `hyb_gga_xc_lc_wpbe` ($\omega = 0.40, s = 0.2, \rho = 1.638 \times 10^{-4}\text{ Bohr}^{-3}$, con `dens_threshold = 1e-30`, `sigma_threshold = 1e-40`): `vrho` salta **`6.076e-6`** relativo; lejos del corte ($\text{aux5} = 18.5, \rho = 1.572 \times 10^{-4}$), la diferencia finita centrada de 4 puntos mantiene un error relativo constante de **`3.90e-6`** para $h_{\text{rel}} \in \{10^{-4}, 10^{-5}, 10^{-6}\}$. En `libxc 7.0.0` (PySCF) el salto era `1.18e-13` y el error FD4 era `7.47e-9`.
  - En `hyb_gga_xc_hse06` ($\omega = 0.11, s = 0.2, \rho = 3.407 \times 10^{-6}\text{ Bohr}^{-3}$): salto relativo de `vrho` = **`4.135e-7`** y error FD4 = **`2.66e-7`** (`6.97e-14` y `1.94e-10` en `7.0.0`).
  - Impacto molecular (PySCF, `def2-TZVP`, malla nivel 4): entre 996 y 2,292 puntos por molécula ($\text{H}_2\text{O}, \text{N}_2, \text{CO}$) caen en $18 < \text{aux5} < 40$; al ocurrir en densidades de cola ($\rho \sim 10^{-4}$–$10^{-6}$), la perturbación integrada sobre elementos de matriz $V_{xc}^{\text{KS}}$ es pequeña ($\max|\Delta V_{\mu\nu}| = 2.63 \times 10^{-10}\text{ E}_h$).
  - Cruce con `ChangeLog.md`: el issue `#629` / `!811` corrigió un `SIGFPE` por desborde en `xc_E1_scaled_jet`, pero conservó el corte prematuro `z > 18.0` (elevarlo a `z > 38.0` elimina el salto hasta precisión de máquina).
- **Reproductor:** `reproducers/repro_01_expint_e1_jet_hse_lcwpbe.py`

---

### Candidato 2 — Tres discontinuidades en `gga_x_wpbeh` / `s_scaling_2` ($s = 1$, $s = 15$ y $w = 14$)
- **Clase tentativa:** Diseño / artefacto de implementación numérica (presente en `7.0.0` y `devel`).
- **Ubicación exacta:** `maple/scaling.mpl:16-19` (`s_scaling_2`) y `maple/gga_exc/gga_x_wpbeh.mpl:66,294` (`eb1` y `term1` con `wcutoff = 14`).
- **Funcionales afectados (prioridad alta):** `hyb_gga_xc_hse06` (428), `hyb_gga_xc_hse03` (427), `hyb_gga_xc_hse12` (479), `hyb_gga_xc_hse12s` (480), `hyb_gga_xc_lc_wpbe` (478), `gga_x_wpbeh` (524).
- **Evidencia medida y hallazgos nuevos:**
  1. **Corte inferior `s = 1.0` (`maple/scaling.mpl:17`):** $\bar{s}$ salta en $-5.140 \times 10^{-4}$. Medido en $\rho = 0.08$: salto relativo persistente ($\delta = 10^{-9}$ y $10^{-14}$) de `9.84e-5` en `zk` y `6.80e-4` en `vsigma` para `HSE06`, y de `1.51e-4` (`zk`) y `1.53e-3` (`vsigma`) para `LC-ωPBE`.
  2. **Corte superior `s = 15.0` (`maple/scaling.mpl:18`, no reportado previamente):** $\bar{s}$ pasa de $15 - \ln(1 + e^{15 - 8.572844}) = 8.5712282$ a `smax = 8.572844`, un salto de **`+1.616e-3` (3.1 veces mayor que en $s = 1$)** mientras $d\bar{s}/ds$ cae de `1.614e-3` a `0` (**100%**). Medido en $\rho = 10^{-4}$: `vsigma` salta **`157.9%` (`1.58`)** en `HSE06`, **`15.1%`** en `LC-ωPBE` y **`100%`** en `gga_x_wpbeh`.
  3. **Corte `w = \omega/k_F = 14.0` (`gga_x_wpbeh.mpl:66,294`, nuevo):** Además del salto entre el ajuste de 8 gaussianas y el de 1 gaussiana (`8.69e-8` en `zk` de `LC-ωPBE`, `7.39e-9` en `HSE06`), `eb1` usa `w < 14` mientras `term1` usa `w > 14`. Exactamente en el flotante `w == 14.0`, `term1` evalúa la rama `w <= 14` pero con el parámetro equivocado `eb1 = 2` en vez de `1.455915...`, generando un **pico puntual adicional** (`5.82e-6` relativo en `devel` y `5.81e-3` en `7.0.0` sobre `HSE06` en $s = 0.5$).
- **Verificación bibliográfica (Regla 6) e impacto molecular (Paso 5):**
  - El PDF tipografiado de Henderson, Izmaylov, Scalmani & Scuseria, *J. Chem. Phys.* **131**, 044108 (2009) está tras muro de pago y **no pudo verificarse directamente** en este entorno. Sin embargo, el comentario en `maple/scaling.mpl:13` atribuye `s_scaling_2` a una primera versión de código de Henderson/Gaussian, mientras que la implementación abierta independiente de **CP2K** (`src/xc/xc_xwpbe.F`) emplea un escalamiento continuo en $s = 8.3$ ($\bar{s} = s$ para $s \le 8.3$ y $\bar{s} = s_{\max} - 18.79622316/s^2$ para $s > 8.3$, sin saltos en $s = 1$ ni $s = 15$, y con `<` / `>=` consistentes en `wcutoff = 14`).
  - Impacto molecular medido (PySCF, `def2-TZVP`, nivel 4): entre **25.4% y 26.8% de la carga electrónica total** de $\text{H}_2\text{O}$, $\text{N}_2$ y $\text{CO}$ reside en $s \in [1, 15]$. Sustituir `s_scaling_2` de `libxc` por la forma continua de CP2K desplaza la energía total de intercambio en **`+0.164` a `+0.275 mE_h` (`+0.10` a `+0.17 kcal/mol`) en HSE06** y **`+0.184` a `+0.324 mE_h` (`+0.12` a `+0.20 kcal/mol`) en LC-ωPBE**.
- **Reproductor:** `reproducers/repro_02_wpbeh_s_scaling_and_w14.py`

---

### Candidato 3 — Corte prematuro `pc07_thr = 1/40` en **r²SCAN-L y SCAN-L** (`mgga_k_pc07.mpl:24-27`)
- **Clase tentativa:** Error numérico nuevo (presente en `7.0.0` y `devel`).
- **Ubicación exacta:** `maple/mgga_exc/mgga_k_pc07.mpl:24-27` (`pc07_thr := 1/40; pc07_zlo := pc07_thr*params_a_pc07_a`).
- **Funcionales afectados (prioridad alta en sólidos/catálisis):** `mgga_x_r2scanl` (718), `mgga_c_r2scanl` (719), `mgga_x_scanl` (700), `mgga_c_scanl` (702), `mgga_x_revscanl` (701), `mgga_c_revscanl` (703), `mgga_k_pc07_opt` (634).
- **Causa aislada (Paso 4):** El comentario de `mgga_k_pc07.mpl:22` afirma que `pc07_fab0(a/40) ~ 1.6e-17 < eps`. Eso era cierto para `mgga_k_pc07` original ($b = 3$, donde vale $4.17 \times 10^{-51}$), pero al reoptimizar **PC07-opt** para construir **SCAN-L y r²SCAN-L** (Mejía-Rodríguez & Trickey 2017–2020), el exponente bajó 11.6 veces hasta **$b = 0.258304$** ($a = 1.784720$) sin escalar `pc07_thr` con $1/b$. En $z_{\text{lo}} = a/40 = 0.044618$, $\text{pc07\_fab0}(z_{\text{lo}}) = \mathbf{4.594 \times 10^{-5}}$, nueve órdenes de magnitud por encima de $\epsilon_{\text{mach}}$, truncándose abruptamente a cero.
- **Evidencia medida (Pasos 1–5):**
  - Control `mgga_k_pc07` ($b = 3$): salto relativo en $z_{\text{lo}} = 0.013473$ de `3.48e-23` en `zk` y `0` en `vlapl`.
  - En `mgga_c_r2scanl` y `mgga_c_scanl` ($\rho = 0.1, \sigma = 0, z_{\text{lo}} = 0.044618 \pm \delta$): salto persistente de **`1.794e-6`** en `zk` y del **`100%`** en `vlapl` (`0` a `-2.655e-6`).
  - En `mgga_x_r2scanl` y `mgga_x_scanl`: salto persistente de **`2.027e-7` / `1.969e-7`** en `zk` y del **`100%`** en `vlapl` (`0` a `+5.365e-6`). En `mgga_k_pc07_opt` con $\sigma = 0$, `zk` y `vlapl` saltan ambos **`100%`**.
  - Impacto molecular (PySCF, `def2-TZVP`): entre 3.0% ($\text{N}_2$) y 5.0% ($\text{H}_2\text{O}$) de la carga electrónica tiene $z_{\text{PC07}} \le 0.044618$ (110 puntos en $|z - z_{\text{lo}}| < 0.005$ en $\text{H}_2\text{O}$). El desplazamiento en energía total integrada de moléculas cerradas en equilibrio es pequeño ($\sim 10^{-10}\text{ E}_h$), pero el potencial generalizado $v_{\nabla^2\rho}$ presenta una discontinuidad de escalón del 100% en la superficie iso-orbital $z_{\text{PC07}} = 0.044618$.
- **Reproductor:** `reproducers/repro_03_pc07_r2scanl_cutoff.py`

---

### Candidato 4 — Divergencia $(1 - |\zeta|)^{-1/3}$ y corte duro `1 - |\zeta| = 1e-12` en **TPSS / TPSSh**
- **Clase tentativa:** Limitación analítica del funcional + discontinuidad por parche duro en `libxc`.
- **Ubicación exacta:** `maple/tpss_c.mpl:34-36` (`my_piecewise3(1 - m_abs(z) <= 1e-12, add(cc[i], i=1..4), ...)`).
- **Funcionales afectados (13 funcionales, prioridad alta):** `mgga_c_tpss` (231), `hyb_mgga_xc_tpssh` (457), `hyb_mgga_xc_tpss0` (248), `mgga_c_revtpss` (241), `hyb_mgga_xc_revtpssh` (458), `mgga_c_tm` (251), etc.
- **Causa aislada y evidencia medida:** Incluso cuando $\nabla\zeta = 0$ exactamente ($\xi^2 = 0$), la derivada parcial $\partial\xi^2/\partial\sigma_{\beta\beta} = 4\rho_\alpha^2/(\rho^4 c_\xi^2)$ es $O(1) \neq 0$. Multiplicada por $(1 - \zeta)^{-4/3}$ en `tpss_C0_den` y por `tpss_par_def` $= O(1 - \zeta)$, produce una divergencia $(1 - \zeta)^{-1/3}$ en $v_{\sigma_{\beta\beta}}$ (reconocida en el comentario de `maple/tpss_c.mpl:28-31`). Al cruzar $1 - |\zeta| = 10^{-12}$ ($\rho_{\text{tot}} = 0.2$), `tpss_C0` pasa a ser constante y **`vsigma_bb` cae abruptamente de `+873.3329 Eh` a `+0.059977 Eh` en `mgga_c_tpss` (y de `+849.35 Eh` a `-23.92 Eh` en `TPSSh`)**, independientemente de `dens_threshold`, `sigma_threshold` y `zeta_threshold`.
- **Reproductor:** `reproducers/repro_04_tpss_c0_ferromagnetic_cut.py`

---

### Candidato 5 — Violación de contrato de API en `xc_func_set_fhc_enforcement` (38 meta-GGAs mixtos)
- **Clase tentativa:** Contrato de API (error nuevo de implementación C).
- **Ubicación exacta:** `src/functionals.c:511-518`.
- **Funcionales afectados (38 meta-GGAs mixtos/híbridos, prioridad alta):** `hyb_mgga_xc_tpssh` (457), `hyb_mgga_xc_r2scan0` (660), `hyb_mgga_x_scan0` (264), `hyb_mgga_xc_tpss0` (248), `hyb_mgga_x_ms2h` (224), `mgga_c_scan_rvv10` (292), `hyb_mgga_xc_pw6b95` (456), `hyb_mgga_xc_mpw1b95` (443), entre otros (los 38 de `n_func_aux > 0` con `XC_FLAGS_ENFORCE_FHC`).
- **Causa aislada y evidencia medida:** Mientras `xc_func_set_{dens,zeta,sigma,tau}_threshold` (`src/functionals.c:463-509`) recorren `for(ii=0; ii<p->n_func_aux; ii++)`, `xc_func_set_fhc_enforcement` modifica únicamente `p->info->flags` del contenedor padre. Como `xc_mix_func` delega el cálculo en `p->func_aux[ii]`, llamar `xc_func_set_fhc_enforcement(p, 0)` (introducido expresamente para códigos PAW donde el seudotensor $\tilde{\tau} < \tilde{\tau}_W$, `ChangeLog.md:189`) produce exactamente **`|dzk| = 0.0` y `|dvtau| = 0.0` en los 38 meta-GGAs mixtos** (frente a `|dvtau| = 4.05e-2` en `mgga_x_r2scan` y `5.01e-1` en `mgga_x_tpss` puros en $\rho=0.12, \sigma=0.10, \tau=0.02$).
- **Reproductor:** `reproducers/repro_05_fhc_enforcement_mixed_mgga.py`

---

### Candidatos secundarios (funcionales especializados que superan los pasos 1 a 3)
- **Candidato 6 — `mgga_c_lak` (`maple/mgga_exc/mgga_c_lak.mpl:61,64`, Error numérico nuevo con FHC apagado):** `lak_fc` define explícitamente una rama para $\alpha < 0$, pero `lak_at` en la línea 64 escribe `(a - 1)/m_max(rs*a, XC_EPSILON)` en vez de `m_max(rs*m_abs(a), XC_EPSILON)`. Para todo $\alpha < 0$, `m_max(rs*a, eps)` evalúa a `XC_EPSILON`, congelando `zk = -0.0264748758` constante y provocando un salto del **100%** en `vtau` en $\alpha = 0$ (de `+1.15e-18` a `-6.90e-3`).
- **Candidato 7 — `lda_c_pk09` (`maple/lda_exc/lda_c_pk09.mpl:25,44`, Desborde a `NaN` en cola abierta):** El polinomio denominador de la Ec. (24) tiene una raíz real en $r_s = 22.7394121325$ ($\rho = 2.0304 \times 10^{-5}\text{ Bohr}^{-3}$). La guarda `|f_r_den(rs)| <= 1e-12` cubre solo $|r_s - r_{s,\text{polo}}| \lesssim 10^{-16}$; fuera de ella, $(1 \pm \zeta)^{ss}$ desborda a `inf` para cualquier $\zeta \neq 0$, devolviendo **`NaN` en `zk` y `vrho`** en toda la franja $|r_s - r_{s,\text{polo}}| \lesssim 3 \times 10^{-5}$ (y $|v_\rho| > 10^7\text{ E}_h$ en $|\Delta r_s| = 10^{-4}$).
- **Candidato 8 — `gga_x_lb` (`maple/gga_vxc/gga_x_lb.mpl:31`) y `mgga_x_br89*` (`maple/mgga_exc/mgga_x_br89.mpl:36`):**
  - En **LB94** (`gga_x_lb`, `gga_x_lbm`), el cambio a la asíntota en $x_s = 300$ descarta prematuramente el `1 +` de $1 + 3\beta x_s \operatorname{asinh}(\gamma x_s) \approx 306$, produciendo un salto persistente de **`0.322%` (`3.218e-3`)** en `vrho`.
  - En **Becke–Roussel** (`mgga_x_br89`, `mgga_x_br89_explicit`, `mgga_x_b00`, `mgga_c_b94`), `br89_cQ` recorta $|Q| < 5 \times 10^{-13}$ a $\pm 5 \times 10^{-13}$, anulando `vlapl` y `vtau` en esa franja (y creando un salto de `3.10e-7` en $Q = 0$ en `br89_explicit` por ligera asimetría del aproximante de Padé).
- **Reproductor agrupado:** `reproducers/repro_06_secondary_lak_pk09_lb94_br89.py`
